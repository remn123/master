
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

// explicit instances
// template class SD<Euler>;
// template class SD<NavierStokes>;

// CONSTRUCTORS
TEST_CASE("1: Test Quadrangle - constructor", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});

  REQUIRE(q->id == 0);
  REQUIRE(q->boundary == 0);
  REQUIRE(q->fringe == 0);
  //REQUIRE(q->J[0][0] == 0.0);
  REQUIRE(q->nodes.size() == 4);
  REQUIRE(q->nodes[0] == 0);
  REQUIRE(q->nodes[1] == 1);
  REQUIRE(q->nodes[2] == 2);
  REQUIRE(q->nodes[3] == 3);

}

TEST_CASE("3: Test Quadrangle - constructor", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});

  REQUIRE(q->id == 0);
  REQUIRE(q->boundary == 0);
  REQUIRE(q->fringe == 0);
  //REQUIRE(q->J[0][0] == 0.0);
  REQUIRE(q->nodes.size() == 4);
  REQUIRE(q->nodes[0] == 0);
  REQUIRE(q->nodes[1] == 1);
  REQUIRE(q->nodes[2] == 2);
  REQUIRE(q->nodes[3] == 3);
  
}

TEST_CASE("2: Test Quadrangle - calculate_jacobian", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});
  std::vector<Vertice> enodes = {{0.724138, 0.724138, 0.0},
                                 {0.758621, 0.724138, 0.0},
                                 {0.758621, 0.758621, 0.0},
                                 {0.724138, 0.758621, 0.0}};

  REQUIRE(enodes[0].coords[0] == Approx(0.724138).margin(1E-15));
  REQUIRE(enodes[0].coords[1] == Approx(0.724138).margin(1E-15));

  REQUIRE(enodes[1].coords[0] == Approx(0.758621).margin(1E-15));
  REQUIRE(enodes[1].coords[1] == Approx(0.724138).margin(1E-15));

  REQUIRE(enodes[2].coords[0] == Approx(0.758621).margin(1E-15));
  REQUIRE(enodes[2].coords[1] == Approx(0.758621).margin(1E-15));

  REQUIRE(enodes[3].coords[0] == Approx(0.724138).margin(1E-15));
  REQUIRE(enodes[3].coords[1] == Approx(0.758621).margin(1E-15));

  int order = 2;                                           // second order
  int dimension = 2;                                       // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  sd->create_nodes();
  q->allocate_jacobian(order);
  q->calculate_jacobian(sd->snodes, sd->fnodes, enodes);

  //REQUIRE(q->J[0][0] == Approx(0.0002972693).margin(1E-15)); // Jacobian
  REQUIRE(q->J[0][0] == Approx(0.0047563092).margin(1E-15)); // Jacobian

  /* *************** SP ****************** */
  // SP 1
  // Jm
  REQUIRE(q->Jm[0][0][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[0][0][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[0][0][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[0][0][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][0][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[0][0][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[0][0][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[0][0][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // SP 2
  // Jm
  REQUIRE(q->Jm[0][1][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[0][1][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[0][1][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[0][1][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][1][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[0][1][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[0][1][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[0][1][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // SP 3
  // Jm
  REQUIRE(q->Jm[0][2][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[0][2][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[0][2][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[0][2][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][2][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[0][2][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[0][2][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[0][2][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // SP 4
  // Jm
  REQUIRE(q->Jm[0][3][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[0][3][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[0][3][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[0][3][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][3][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[0][3][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[0][3][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[0][3][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPs
  /* *************** FP X ****************** */
  // FPx 1
  // Jm
  REQUIRE(q->Jm[1][0][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][0][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][0][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][0][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][0][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][0][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][0][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][0][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPx 2
  // Jm
  REQUIRE(q->Jm[1][1][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][1][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][1][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][1][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][1][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][1][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][1][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][1][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPx 3
  // Jm
  REQUIRE(q->Jm[1][2][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][2][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][2][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][2][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][2][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][2][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][2][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][2][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPx 4
  // Jm
  REQUIRE(q->Jm[1][3][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][3][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][3][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][3][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][3][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][3][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][3][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][3][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPx 5
  // Jm
  REQUIRE(q->Jm[1][4][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][4][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][4][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][4][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][4][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][4][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][4][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][4][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPx 6
  // Jm
  REQUIRE(q->Jm[1][5][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][5][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][5][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][5][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][5][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][5][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][5][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][5][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  /* *************** FP Y ****************** */
  // FPy 1
  // Jm
  REQUIRE(q->Jm[2][0][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][0][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][0][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][0][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][0][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][0][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][0][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][0][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPy 2
  // Jm
  REQUIRE(q->Jm[2][1][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][1][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][1][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][1][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][1][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][1][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][1][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][1][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPy 3
  // Jm
  REQUIRE(q->Jm[2][2][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][2][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][2][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][2][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][2][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][2][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][2][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][2][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPy 4
  // Jm
  REQUIRE(q->Jm[2][3][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][3][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][3][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][3][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][3][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][3][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][3][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][3][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPy 5
  // Jm
  REQUIRE(q->Jm[2][4][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][4][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][4][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][4][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][4][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][4][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][4][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][4][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy

  // FPy 6
  // Jm
  REQUIRE(q->Jm[2][5][0] == Approx(0.068966).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][5][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][5][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][5][3] == Approx(0.068966).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][5][0] == Approx(14.4998985007).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][5][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][5][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][5][3] == Approx(14.4998985007).margin(1E-15)); // deta_dy
}

// CONSTRUCTORS
// TEST_CASE("1: Test Quadrangle - constructor", "[elems]")
// {
//   auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});

//   REQUIRE(q->id == 0);
//   REQUIRE(q->boundary == 0);
//   REQUIRE(q->fringe == 0);
//   REQUIRE(q->J[0][0] == 0.0);
//   REQUIRE(q->nodes.size() == 4);
//   REQUIRE(q->nodes[0] == 0);
//   REQUIRE(q->nodes[1] == 1);
//   REQUIRE(q->nodes[2] == 2);
//   REQUIRE(q->nodes[3] == 3);
// }
