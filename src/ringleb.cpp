#include <iostream>
#include <filesystem>
#include <string>
#include <map>
#include <memory>
//#include <crtdbg.h>

#include <Mesh.h>
#include <Solver.h>

namespace fs = std::filesystem;

int main()
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / ".." / "resources" / "ringleb.msh").string());

  int order = 2;
  auto sd = std::make_shared<SD<Euler>>(2, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  sd->setup(mesh);

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
  int stages = 5;
  int rk_order = 4;
  int size = mesh->Nel * (order * order); // overall number of solution points

  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>({CFL,
                                                               MAX_ITER,
                                                               stages,
                                                               rk_order,
                                                               size});

  time->loop(mesh, sd->solve);

  time->save(mesh, sd->to_vtk);

  //std::cin.get();

  return 0;
}
