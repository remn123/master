#include <iostream>
#include <filesystem>
#include <functional>
#include <string>
#include <map>
#include <memory>

#include <Dummies.h>
#include <Mesh.h>
#include <SD.h>
#include <Time.h>

namespace fs = std::filesystem;

int main()
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "vortex" / "vortex.msh").string());

  int p = 3;
  int order = p+1;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  auto gamma = 1.4;
  
  Ghost::Mach = 0.05;
  Ghost::U = Ghost::Mach;
  Ghost::V = 0.0;
  Ghost::T = 300.0;
  //Ghost::p = 10E+5;
  Ghost::R = 287.15;

  Ghost::rho = 1.0;
  Ghost::p = 1.0/gamma;
  
  sd->setup(
    mesh, 
    FIELDS::VORTEX_FIELD_MAPPING
  );

  double min_dx = 0.1;
  //double min_dx = sd->get_min_dx(mesh);
  // std::cout << "Minimum dx = " << min_dx << "\n";

  auto filename_msh = (cur_path.parent_path() / "results" / "output" / "vortex" / "3" / "mesh" / "vortex.vtk"  ).string();
  mesh->to_vtk(filename_msh);


  std::cout << "Saving Initial Condition ...\n";
  auto filename = (cur_path.parent_path() / "results" / "output" / "vortex" / "3" / "iterations" / "pp_vortex_").string();
  std::string tstamp = std::to_string(0);
  tstamp.insert(tstamp.begin(), 5 - tstamp.length(), '0');
  sd->to_vtk(mesh, filename + tstamp + std::string{".vtk"});

  /*
    2) Solver Loop (for each element in Mesh)
       2.1) Apply Boundary Conditions on interfaces (flux points)
       2.2) Interpolate Q and dQ from SPs to FPs
       2.3) Calculate Fluxes at internal FPs
       2.4) Calculate Fluxes at faces FPs (riemann solver and BR2)
       2.5) Interpolate Fluxes derivatives from FPs to SPs
       2.6) Calculate Residue
  */
  std::cout << "Running solver first step\n";
  sd->solve(mesh);

  /*
    3) Time Marching Loop
      3.1) Calculate residue norm
      3.2) Check if it's already converged
      3.3) (if not) Apply time iteration then go to (2)
  */
  double CFL = 0.000001;
  long MAX_ITER = 1000;
  int rk_order = 3;
  int stages = 3;
  int size = mesh->Nel * (order * order)*4; // overall number of solution points

  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>(
    CFL, MAX_ITER, stages, rk_order, size, p, Ghost::U, min_dx
  );
  
  std::cout << "Time integration\n";
  time->loop(
    mesh, 
    [&sd](std::shared_ptr<Mesh> & m){sd->solve(m);},
    filename,
    [&sd](const std::shared_ptr<Mesh> & m, const std::string & f){sd->to_vtk(m, f);}
  );

  
  time->save(
    mesh, 
    filename, 
    [&sd](const std::shared_ptr<Mesh> & m, const std::string & f){sd->to_vtk(m, f);}
  );

  return 0;
}
