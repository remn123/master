#include <iostream>
#include <memory>
#include <sstream>

#include <catch/catch.hpp>
#include <Dummies.h>
#include <Element.h>
#include <Mesh.h>
#include <Node.h>
#include <SD.h>

using namespace Catch::literals;

// explicit instantiations
template class SD<Euler>;
template class SD<NavierStokes>;

// CONSTRUCTORS
TEST_CASE("1: Test SD<Euler> - Constructor", "[sd]")
{
  int order=2;     // second order
  int dimension=2; // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  
  REQUIRE(sd->order == order);
  REQUIRE(sd->dimension == dimension);
  REQUIRE(sd->MU == 0.0);
  REQUIRE(sd->MUv == 0.0);

  REQUIRE(sd->snodes.size() == (order*order));
  REQUIRE(sd->fnodes.size() == dimension);

  for (auto i=0; i<dimension; i++)
    REQUIRE(sd->fnodes[i].size() == (order+1)*order);

  for (auto& sn: sd->snodes)
  {
    REQUIRE(sn.id == -1);
    REQUIRE(sn.coords[0] == 0.0);
    REQUIRE(sn.coords[1] == 0.0);
    REQUIRE(sn.coords[2] == 0.0);
    REQUIRE(sn.left == -1);
    REQUIRE(sn.right == -1);
  }

  for (auto& vec: sd->fnodes)
  { 
    for (auto& fn: vec)
    {
      REQUIRE(fn.id == -1);
      REQUIRE(fn.coords[0] == 0.0);
      REQUIRE(fn.coords[1] == 0.0);
      REQUIRE(fn.coords[2] == 0.0);
      REQUIRE(fn.left == -1);
      REQUIRE(fn.right == -1);
    }
  }
}


TEST_CASE("2: Test SD<Euler> - create_nodes", "[sd]")
{
  int order=2;     // second order
  int dimension=2; // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  
  sd->create_nodes();

  // SNODES
  // 0
  REQUIRE(sd->snodes[0].id == -1);
  REQUIRE(sd->snodes[0].coords[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[0].coords[1] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[0].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->snodes[0].left == -1);
  REQUIRE(sd->snodes[0].right == -1);

  // 1
  REQUIRE(sd->snodes[1].id == -1);
  REQUIRE(sd->snodes[1].coords[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[1].coords[1] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[1].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->snodes[1].left == -1);
  REQUIRE(sd->snodes[1].right == -1);

  // 2
  REQUIRE(sd->snodes[2].id == -1);
  REQUIRE(sd->snodes[2].coords[0] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[2].coords[1] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[2].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->snodes[2].left == -1);
  REQUIRE(sd->snodes[2].right == -1);

  // 3
  REQUIRE(sd->snodes[3].id == -1);
  REQUIRE(sd->snodes[3].coords[0] == Approx(0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[3].coords[1] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[3].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->snodes[3].left == -1);
  REQUIRE(sd->snodes[3].right == -1);

  // FNODES - x
  // 0
  REQUIRE(sd->fnodes[0][0].id == -1);
  REQUIRE(sd->fnodes[0][0].coords[0] == -1.0);
  REQUIRE(sd->fnodes[0][0].coords[1] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][0].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[0][0].left == -1);
  REQUIRE(sd->fnodes[0][0].right == -1);

  // 1
  REQUIRE(sd->fnodes[0][1].id == -1);
  REQUIRE(sd->fnodes[0][1].coords[0] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[0][1].coords[1] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][1].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[0][1].left == -1);
  REQUIRE(sd->fnodes[0][1].right == -1);

  // 2
  REQUIRE(sd->fnodes[0][2].id == -1);
  REQUIRE(sd->fnodes[0][2].coords[0] == 1.0);
  REQUIRE(sd->fnodes[0][2].coords[1] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][2].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[0][2].left == -1);
  REQUIRE(sd->fnodes[0][2].right == -1);

  // 3
  REQUIRE(sd->fnodes[0][3].id == -1);
  REQUIRE(sd->fnodes[0][3].coords[0] == -1.0);
  REQUIRE(sd->fnodes[0][3].coords[1] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][3].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[0][3].left == -1);
  REQUIRE(sd->fnodes[0][3].right == -1);

  // 4
  REQUIRE(sd->fnodes[0][4].id == -1);
  REQUIRE(sd->fnodes[0][4].coords[0] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[0][4].coords[1] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][4].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[0][4].left == -1);
  REQUIRE(sd->fnodes[0][4].right == -1);

  // 5
  REQUIRE(sd->fnodes[0][5].id == -1);
  REQUIRE(sd->fnodes[0][5].coords[0] == 1.0);
  REQUIRE(sd->fnodes[0][5].coords[1] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][5].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[0][5].left == -1);
  REQUIRE(sd->fnodes[0][5].right == -1);

  // FNODES - y
  // 0
  REQUIRE(sd->fnodes[1][0].id == -1);
  REQUIRE(sd->fnodes[1][0].coords[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][0].coords[1] == -1.0);
  REQUIRE(sd->fnodes[1][0].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[1][0].left == -1);
  REQUIRE(sd->fnodes[1][0].right == -1);

  // 1
  REQUIRE(sd->fnodes[1][1].id == -1);
  REQUIRE(sd->fnodes[1][1].coords[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][1].coords[1] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[1][1].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[1][1].left == -1);
  REQUIRE(sd->fnodes[1][1].right == -1);

  // 2
  REQUIRE(sd->fnodes[1][2].id == -1);
  REQUIRE(sd->fnodes[1][2].coords[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][2].coords[1] == 1.0);
  REQUIRE(sd->fnodes[1][2].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[1][2].left == -1);
  REQUIRE(sd->fnodes[1][2].right == -1);

  // 3
  REQUIRE(sd->fnodes[1][3].id == -1);
  REQUIRE(sd->fnodes[1][3].coords[0] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][3].coords[1] == -1.0);
  REQUIRE(sd->fnodes[1][3].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[1][3].left == -1);
  REQUIRE(sd->fnodes[1][3].right == -1);

  // 4
  REQUIRE(sd->fnodes[1][4].id == -1);
  REQUIRE(sd->fnodes[1][4].coords[0] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][4].coords[1] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[1][4].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[1][4].left == -1);
  REQUIRE(sd->fnodes[1][4].right == -1);

  // 5
  REQUIRE(sd->fnodes[1][5].id == -1);
  REQUIRE(sd->fnodes[1][5].coords[0] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][5].coords[1] == 1.0);
  REQUIRE(sd->fnodes[1][5].coords[2] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[1][5].left == -1);
  REQUIRE(sd->fnodes[1][5].right == -1);
}