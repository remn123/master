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
  auto background_msh = std::make_shared<Static_Mesh>(2);
  auto nearbody_msh = std::make_shared<Static_Mesh>(2);

  // background_msh->read_gmsh((cur_path.parent_path() / "resources" / "overset" / "cylinder" / "cylinder_r100_bkgd_64x64.msh").string());
  // nearbody_msh->read_gmsh((cur_path.parent_path() / "resources" / "overset" / "cylinder" / "cylinder_r100_body_64x64.msh").string());
  background_msh->read_gmsh((cur_path.parent_path() / "resources" / "overset" / "cylinder" / "cylinder_r100_bkgd_64x64.msh").string());
  nearbody_msh->read_gmsh((cur_path.parent_path() / "resources" / "overset" / "cylinder" / "cylinder_r100_body_64x64_ho.msh").string());
  
  // Creating kd-tree
  background_msh->create_kdtree();
  nearbody_msh->create_kdtree();

  int order = 3;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  Ghost::Mach = 0.2;
  Ghost::p = 1.0/1.4; 
  Ghost::analytical_solution = FIELDS::CYLINDER_FIELD_MAPPING;

  auto background_msh_ = std::static_pointer_cast<Mesh>(background_msh);
  sd->setup(
    background_msh_, 
    FIELDS::CYLINDER_FIELD_MAPPING
  );

  auto nearbody_msh_ = std::static_pointer_cast<Mesh>(nearbody_msh);
  sd->setup(
    nearbody_msh_, 
    FIELDS::CYLINDER_FIELD_MAPPING
  );

  sd->update_overset(background_msh, nearbody_msh);
  
  sd->communicate_data(background_msh, nearbody_msh);
  sd->communicate_data(nearbody_msh, background_msh);
  
  // auto filename_bkgd = (cur_path.parent_path() / "results" / "overset" / "low"  / "fringes" / "background.vtk" ).string();
  // auto filename_body = (cur_path.parent_path() / "results" / "overset" / "low" /  "fringes" / "near_body.vtk"  ).string();
  
  auto filename_bkgd = (cur_path.parent_path() / "results" / "overset" / "ho_v2"  / "fringes" / "background.vtk" ).string();
  auto filename_body = (cur_path.parent_path() / "results" / "overset" / "ho_v2" /  "fringes" / "near_body.vtk"  ).string();

  background_msh->to_vtk(filename_bkgd);
  nearbody_msh->to_vtk(filename_body);

  // std::cout << "Saving Initial Condition ...\n";
  // filename_bkgd = (cur_path.parent_path() / "results" / "overset" / "low" / "iterations" / "pp_cylinder_bkgd_").string();
  // filename_body = (cur_path.parent_path() / "results" / "overset" / "low" / "iterations" / "pp_cylinder_body_").string();
  filename_bkgd = (cur_path.parent_path() / "results" / "overset" / "ho_v2" / "iterations" / "pp_cylinder_bkgd_").string();
  filename_body = (cur_path.parent_path() / "results" / "overset" / "ho_v2" / "iterations" / "pp_cylinder_body_").string();
  std::string tstamp = std::to_string(0);
  tstamp.insert(tstamp.begin(), 5 - tstamp.length(), '0');
  sd->to_vtk(background_msh_, filename_bkgd + tstamp + std::string{".vtk"});
  sd->to_vtk(nearbody_msh_, filename_body + tstamp + std::string{".vtk"});

  /*
    2) Solver Loop (for each element in Mesh)
       2.1) Apply Boundary Conditions on interfaces (flux points)
       2.2) Interpolate Q and dQ from SPs to FPs
       2.3) Calculate Fluxes at internal FPs
       2.4) Calculate Fluxes at faces FPs (riemann solver and BR2)
       2.5) Interpolate Fluxes derivatives from FPs to SPs
       2.6) Calculate Residue
  */
  std::cout << "Running solver first step for the Background mesh\n";
  sd->solve(background_msh_);
  std::cout << "Running solver first step for the Near-Body mesh\n";
  sd->solve(nearbody_msh_);

  /*
    3) Time Marching Loop
      3.1) Calculate residue norm
      3.2) Check if it's already converged
      3.3) (if not) Apply time iteration then go to (2)
  */
  double CFL = 3.0;
  long MAX_ITER = 3E+4;
  int rk_order = 3;
  int stages = 3;
  int size = (background_msh->Nel + nearbody_msh->Nel)*(order * order)*4; // overall number of solution points

  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>(CFL,
                                                              MAX_ITER,
                                                              stages,
                                                              rk_order,
                                                              size);
  
  std::cout << "Time integration\n";
  time->loop(
    background_msh, 
    nearbody_msh, 
    [&sd](std::shared_ptr<Mesh> & m){sd->solve(m);},
    [&sd](std::shared_ptr<Static_Mesh> & r, const std::shared_ptr<Static_Mesh> & d){sd->communicate_data(r, d);},
    filename_bkgd,
    filename_body,
    [&sd](const std::shared_ptr<Mesh> & m, const std::string & f){sd->to_vtk(m, f);}
  );

  
  time->save(
    background_msh_, 
    filename_bkgd, 
    [&sd](const std::shared_ptr<Mesh> & m, const std::string & f){sd->to_vtk(m, f);}
  );
  time->save(
    nearbody_msh_, 
    filename_body, 
    [&sd](const std::shared_ptr<Mesh> & m, const std::string & f){sd->to_vtk(m, f);}
  );

  return 0;
}
