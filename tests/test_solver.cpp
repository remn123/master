
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

using namespace Catch::literals;

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
       3.1) Apply Boundary Conditions on interfaces (flux points)
      
  */
  double CFL = 0.1;
  auto time = std::make_shared<TD>(2, CFL);

  time->setup();

  time->update(mesh, sd);

  time->save(mesh, sd);
  // REQUIRE( sd.snodes[0] == Approx(-0.577350269189626).margin(1E-15));
  // REQUIRE( sd.snodes[1] == Approx(0.577350269189626).margin(1E-15));
}