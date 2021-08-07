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

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "airfoil" / "airfoil.msh").string());

  int order = 3;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  //Ghost::p = (1.0/1.4);
  auto gamma = 1.4;
  auto rho_ref = 1.225;
  auto temp_ref= 288.203086112494;
  auto vel_ref = 340.2939905434710;
  auto gas_ref = temp_ref / std::pow(vel_ref, 2.0);
  auto ene_ref = rho_ref * std::pow(vel_ref, 2.0);
  auto R = 287.0;
  auto Cv = R/(gamma-1.0);

  Ghost::R = R* gas_ref;
  Ghost::Cv = Cv * gas_ref;

  Ghost::Mach = 0.2;
  Ghost::U = 0.2;
  Ghost::V = 0.0;
  //Ghost::T = (1.0 + Ghost::Mach*Ghost::Mach*(gamma-1.0)/2.0));
  Ghost::p = 1.0/gamma;
  Ghost::rho = 1.0;
  Ghost::T = (1.0/(gamma*Ghost::R));
  //Ghost::analytical_solution = FIELDS::CYLINDER_FIELD_MAPPING;

  sd->setup(
    mesh, 
    FIELDS::DEFAULT_FIELD_MAPPING
  );

  std::cout << "Saving Initial Condition ...\n";
  auto filename = (cur_path.parent_path() / "results" / "airfoil" / "pp_airfoil_").string();
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
  double CFL = 0.5;
  long MAX_ITER = 5E+5;
  int rk_order = 3;
  int stages = 3;
  int size = mesh->Nel * (order * order)*4; // overall number of solution points

  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>(CFL,
                                                              MAX_ITER,
                                                              stages,
                                                              rk_order,
                                                              size);
  
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
