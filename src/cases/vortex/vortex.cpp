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

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "overset" / "vortex" / "msh" / "vortex_64_bkg.msh").string());

  int p = 3;
  int order = p+1;
  fs::path result_folder_path = "results/output/single/vortex_64/p=" + std::to_string(order);
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  // auto gamma = 1.4;
  
  // Ghost::Mach = 0.05;  
  // Ghost::T = 300.0;
  // Ghost::p = 1.0E+5;
  // Ghost::R = 287.15;
  // Ghost::U = Ghost::Mach*std::sqrt(gamma*Ghost::R*Ghost::T);
  // Ghost::V = 0.0;

  // Ghost::rho = Ghost::p/(Ghost::R*Ghost::T);

  Ghost::T   = 1.0;
  Ghost::R   = 1.0;
  Ghost::rho = 1.0;
  Ghost::p   = 1.0;
  Ghost::U   = 1.0;
  Ghost::V   = 1.0;
  Ghost::Mach = 1.0/std::sqrt(1.4);
  Ghost::time = 0.0;
  double min_dx = 5.0;
  double Umag = std::sqrt(std::pow(Ghost::U, 2.0)+std::pow(Ghost::V, 2.0));
  Ghost::analytical_solution = FIELDS::VORTEX_T_FIELD_MAPPING;

  sd->setup(
    mesh, 
    FIELDS::VORTEX_FIELD_MAPPING
  );

  fs::path results_root_path = cur_path.parent_path() / result_folder_path;

  if (!fs::is_directory(results_root_path) || !fs::exists(results_root_path)) 
  {
    fs::create_directories((results_root_path/"mesh"));
    fs::create_directories((results_root_path/"iterations"));
    fs::create_directories((results_root_path/"log"));
  }

  //double min_dx = 0.1;
  //double min_dx = sd->get_min_dx(mesh);
  std::cout << "Minimum dx = " << min_dx << "\n";

  auto filename_msh = (cur_path.parent_path() / results_root_path / "mesh" / "vortex.vtk"  ).string();
  mesh->to_vtk(filename_msh);


  std::cout << "Saving Initial Condition ...\n";
  auto filename = (cur_path.parent_path() / results_root_path / "iterations" / "pp_vortex_").string();
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
  double CFL = 1.0/1000.0;
  long MAX_ITER = 600;
  int rk_order = 3;
  int stages = 3;
  int size = mesh->Nel * (order * order)*4; // overall number of solution points

  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>(
    CFL, MAX_ITER, stages, rk_order, 
    size, p, Umag, min_dx, results_root_path.string()
  );
  
  std::cout << "Time integration\n";
  time->loop(
    mesh, 
    [&sd](std::shared_ptr<Mesh> & m){sd->solve(m);},
    [&sd](std::shared_ptr<Mesh> & m, double t)->std::vector<double>{return sd->get_property_error_at_time(m, t);},
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
