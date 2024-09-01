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

  background_msh->read_gmsh((cur_path.parent_path() / "resources" / "overset" / "vortex" / "msh" / "vortex_16_bkg.msh").string());
  nearbody_msh->read_gmsh((cur_path.parent_path() / "resources" / "overset" / "vortex" / "msh" / "vortex_5_bdy.msh").string());
  // Creating kd-tree
  background_msh->create_kdtree();
  nearbody_msh->create_kdtree();

  int p = 5;
  int order = p+1;
  fs::path result_folder_path = "results/output/overset/vortex_16_5/p=" + std::to_string(order);
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  
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
  
  auto background_msh_ = std::static_pointer_cast<Mesh>(background_msh);
  sd->setup(
    background_msh_, 
    FIELDS::VORTEX_FIELD_MAPPING
  );

  auto nearbody_msh_ = std::static_pointer_cast<Mesh>(nearbody_msh);
  sd->setup(
    nearbody_msh_, 
    FIELDS::VORTEX_FIELD_MAPPING
  );

  fs::path results_root_path = cur_path.parent_path() / result_folder_path;

  if (!fs::is_directory(results_root_path) || !fs::exists(results_root_path)) 
  {
    fs::create_directories((results_root_path/"fringes"));
    fs::create_directories((results_root_path/"iterations"));
    fs::create_directories((results_root_path/"log"));
  }

  sd->update_overset(background_msh, nearbody_msh);
  
  sd->communicate_data(background_msh, nearbody_msh);
  sd->communicate_data(nearbody_msh, background_msh);
  
  
  //double min_dx = sd->get_min_dx(nearbody_msh_);
  std::cout << "Minimum dx = " << min_dx << "\n";
  
  auto filename_bkgd = (results_root_path / "fringes" / "background.vtk" ).string();
  auto filename_body = (results_root_path / "fringes" / "near_body.vtk"  ).string();

  background_msh->to_vtk(filename_bkgd);
  nearbody_msh->to_vtk(filename_body);

  std::cout << "Saving Initial Condition ...\n";
  filename_bkgd = (results_root_path / "iterations" / "pp_vortex_bkgd_").string();
  filename_body = (results_root_path / "iterations" / "pp_vortex_body_").string();
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
  double CFL = 1.0/1000.0;
  long MAX_ITER = 600;
  int rk_order = 3;
  int stages = 3;
  int size = (background_msh->Nel + nearbody_msh->Nel)*(order * order)*4; // overall number of solution points

  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>(
    CFL, MAX_ITER, stages, rk_order, 
    size, p, Umag, min_dx, results_root_path.string()
  );
  
  std::cout << "Time integration\n";
  // CALLGRIND_START_INSTRUMENTATION;
  // CALLGRIND_TOGGLE_COLLECT;
  try {
    time->loop(
      background_msh, 
      nearbody_msh, 
      [&sd](std::shared_ptr<Mesh> & m){sd->solve(m);},
      [&sd](std::shared_ptr<Mesh> & m, double t)->std::vector<double>{return sd->get_property_error_at_time(m, t);},
      [&sd](std::shared_ptr<Static_Mesh> & r, const std::shared_ptr<Static_Mesh> & d){sd->communicate_data(r, d);},
      filename_bkgd,
      filename_body,
      [&sd](const std::shared_ptr<Mesh> & m, const std::string & f){sd->to_vtk(m, f);}
    );
  } catch (const char* msg) {
    std::cerr << msg;
    throw;
  }
  // CALLGRIND_TOGGLE_COLLECT;
  // CALLGRIND_STOP_INSTRUMENTATION;
  
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
