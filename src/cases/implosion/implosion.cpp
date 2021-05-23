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

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "mesh_implosion.msh").string());
  

  int order = 2;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  //Ghost::analytical_solution = FIELDS::IMPLOSION_FIELD_MAPPING;

  sd->setup(
    mesh, 
    FIELDS::IMPLOSION_FIELD_MAPPING
  );

  auto filename = (cur_path.parent_path() / "results" / "implosion" / "pp_mesh_implosion_").string();
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

  sd->solve(mesh);

  /*
    3) Time Marching Loop
      3.1) Calculate residue norm
      3.2) Check if it's already converged
      3.3) (if not) Apply time iteration then go to (2)
  */
  double CFL = 0.01;
  long MAX_ITER = 1E+4;
  int rk_order = 3;
  int stages = 3;
  int size = mesh->Nel * (order * order)*4; // overall number of solution points

  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>(CFL,
                                                              MAX_ITER,
                                                              stages,
                                                              rk_order,
                                                              size);
  
  
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
