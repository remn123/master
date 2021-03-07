
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <catch/catch.hpp>
#include <Dummies.h>
#include <Helpers.h>
#include <Poly.h>
#include <Mesh.h>
#include <SD.h>
#include <Time.h>

using namespace Catch::literals;

// explicit instantiations
template class SD<Euler>;
template class SD<NavierStokes>;


namespace fs = std::filesystem;
// Solver
/*
  1) Setup (all element in Mesh)
     1.1) Calculate solution and fluxes points
     1.2) Initialize solution and fluxes

  2) Solver Loop (for each element in Mesh)
     2.1) Apply Boundary Conditions on interfaces (flux points)
     2.2) Interpolate Q and dQ from SPs to FPs
     2.3) Calculate Fluxes at internal FPs
     2.4) Calculate Fluxes at faces FPs (riemann solver and BR2)
     2.5) Interpolate Fluxes derivatives from FPs to SPs
     2.6) Calculate Residue
  
  3) Iterate in time
     3.1) Calculate residue norm
     3.2) Check if it's already converged
     3.3) (if not) Apply time iteration then go to (2)
*/

TEST_CASE("1: Test Solver - Solution Nodes 2nd Order", "[solver]")
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / ".." / "resources" / "mesh1.msh").string());

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
  double CFL = 0.1, Linf = 0.0;
  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>(params...); // TO DO
  long MAX_ITER = 1E+3;
  long iter = 0;
  while (iter <= MAX_ITER)
  {
    time->update(mesh, sd->solve);

    //L1 = mesh->get_residue_norm(0); // L1-norm
    //L2 = mesh->get_residue_norm(1); // L2-norm
    Linf = mesh->get_residue_norm(2); // Linf-norm
    iter++;
    std::cout << "Iter[" << iter << "]: Residue " << log10(Linf) << std::endl;
  }

  time->save(mesh, sd->to_vtk);
}