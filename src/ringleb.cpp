#include <iostream>
#include <filesystem>
#include <string>
#include <map>
#include <memory>
//#include <crtdbg.h>

#include <Dummies.h>
#include <Mesh.h>
#include <SD.h>
#include <Solver.h>
#include <Time.h>

namespace fs = std::filesystem;

// explicit instantiations
// template class SD<Euler>;
// template class SD<NavierStokes>;

// template class Time<Explicit::SSPRungeKutta>;

/* Declaring analytical field solution */
std::vector<double> ringleb_field (const Node& n)
{
  std::vector<double> vec(4, 0.0);

  // Coordinates
  double x = n.coords[0];
  double y = n.coords[1];
  double q = sqrt(x*x + y*y); // Umag
  double gamma = 1.4;

  double a = sqrt(1.0 - (gamma-1.0)*q*q/2.0);

  double rho = pow(a, (2.0/(gamma-1.0)));
  double u = x;
  double v = y;
  double p = (1.0/gamma)*pow(a, (2.0*gamma/(gamma-1.0)));
  double E = p/(gamma-1.0) + rho*(x*x + y*y)/2.0;

  vec[0] = rho; 
  vec[1] = rho*u; 
  vec[2] = rho*v; 
  vec[3] = E; 
  
  return vec;

}

int main()
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "ringleb.msh").string());

  int order = 2;
  auto sd = std::make_shared<SD<Euler>>(2, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  sd->setup(mesh, ringleb_field);

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
  double CFL = 0.1;
  long MAX_ITER = 1E+3;
  int rk_order = 4;
  int stages = 5;
  int size = mesh->Nel * (order * order); // overall number of solution points

  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>(CFL,
                                                              MAX_ITER,
                                                              stages,
                                                              rk_order,
                                                              size);

  time->loop(mesh, sd->solve);

  auto filename = (cur_path.parent_path() / "results" / "pp_mesh_ringleb.vtk").string();
  time->save(mesh, filename, sd->to_vtk);

  return 0;
}
