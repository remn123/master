#include <iostream>
#include <filesystem>
#include <functional>
#include <string>
#include <map>
#include <memory>
//#include <crtdbg.h>

#include <Dummies.h>
#include <Mesh.h>
#include <SD.h>
//#include <Solver.h>
#include <Time.h>

namespace fs = std::filesystem;

// explicit instantiations
// template class SD<Euler>;
// template class SD<NavierStokes>;

// template class Time<Explicit::SSPRungeKutta>;

/* Declaring analytical field solution */

int main()
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "ringleb_v0.msh").string());
  
  int p = 1;
  int order = p+1;
  fs::path result_folder_path = "results/output/single/ringleb_3/p=" + std::to_string(order);
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  //auto field = [](const Node& n){return FIELDS::RINGLEB_FIELD_MAPPING(n);};
  Ghost::analytical_solution = FIELDS::RINGLEB_FIELD_MAPPING;
  
  Ghost::U = 0.5;
  Ghost::V = 0.0;
  Ghost::p = 1.0/1.4;
  Ghost::rho = 1.0;

  sd->setup(
    mesh, 
    FIELDS::RINGLEB_FIELD_MAPPING
  );  
  
  fs::path results_root_path = cur_path.parent_path() / result_folder_path;

  if (!fs::is_directory(results_root_path) || !fs::exists(results_root_path)) 
  {
    fs::create_directories((results_root_path/"boundary"));
    fs::create_directories((results_root_path/"log"));
  }
    

  auto filename1 = (results_root_path / "boundary" / "ringleb_v0.vtk").string();
  std::cout << "Filename1: " << filename1 << "\n";
  mesh->to_vtk(filename1);

  auto filename = (results_root_path / "pp_mesh_ringleb_v0_").string();
  std::cout << "Filename: " << filename << "\n";
  std::string tstamp = std::to_string(0);
  tstamp.insert(tstamp.begin(), 5 - tstamp.length(), '0');
  sd->to_vtk(mesh, filename + tstamp + std::string{".vtk"});
  
  double min_dx = sd->get_min_dx(mesh);
  std::cout << "Minimum dx = " << min_dx << "\n";

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
  double CFL = 0.1;
  long MAX_ITER = 1E+4;
  int rk_order = 3;
  int stages = 3;
  int size = mesh->Nel * (order * order)*4; // overall number of solution points

  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>(
    CFL, MAX_ITER, stages, rk_order, 
    size, p, Ghost::U, min_dx, results_root_path.string()
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
