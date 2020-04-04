#include <iostream>
#include <memory>
#include <sstream>

#include <catch/catch.hpp>
#include <Dummies.h>
#include <Mesh.h>
#include <SD.h>

using namespace Catch::literals;

// CONSTRUCTORS
TEST_CASE("1: Test Quadrangle - constructor", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1 2 3 4"});

  REQUIRE(q->id == 1);
  REQUIRE(q->boundary == 0);
  REQUIRE(q->fringe == 0);
  REQUIRE(q->J == 0.0);
  REQUIRE(q->nodes.size() == 4);
  REQUIRE(q->nodes[0] == 0);
  REQUIRE(q->nodes[1] == 1);
  REQUIRE(q->nodes[2] == 2);
  REQUIRE(q->nodes[3] == 3);
}

TEST_CASE("2: Test Quadrangle - calculate_jacobian", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1 2 3 4"});

  int order=2; // second order
  auto  sd = std::make_shared<SD<Euler>>(2, 2); // (int order, int dimension)

  std::vector<Node> enodes = {};

  q->calculate_jacobian(sd->snodes, sd->fnodes, enodes);

  REQUIRE(q->J == 0.0); // Jacobian
  
  /* *************** SP ****************** */
  // SP 1
  // Jm
  REQUIRE(q->Jm[0][0][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[0][0][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[0][0][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[0][0][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][0][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[0][0][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[0][0][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[0][0][3] == 0.0); // deta_dy

  // SP 2
  // Jm
  REQUIRE(q->Jm[0][1][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[0][1][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[0][1][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[0][1][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][1][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[0][1][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[0][1][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[0][1][3] == 0.0); // deta_dy

  // SP 3
  // Jm
  REQUIRE(q->Jm[0][2][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[0][2][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[0][2][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[0][2][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][2][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[0][2][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[0][2][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[0][2][3] == 0.0); // deta_dy

  // SP 4
  // Jm
  REQUIRE(q->Jm[0][3][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[0][3][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[0][3][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[0][3][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][3][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[0][3][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[0][3][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[0][3][3] == 0.0); // deta_dy

  // FPs 
  /* *************** FP X ****************** */
  // FPx 1
  // Jm
  REQUIRE(q->Jm[1][0][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[1][0][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[1][0][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[1][0][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][0][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[1][0][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[1][0][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[1][0][3] == 0.0); // deta_dy

  // FPx 2
  // Jm
  REQUIRE(q->Jm[1][1][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[1][1][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[1][1][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[1][1][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][1][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[1][1][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[1][1][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[1][1][3] == 0.0); // deta_dy

  // FPx 3
  // Jm
  REQUIRE(q->Jm[1][2][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[1][2][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[1][2][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[1][2][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][2][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[1][2][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[1][2][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[1][2][3] == 0.0); // deta_dy

  // FPx 4
  // Jm
  REQUIRE(q->Jm[1][3][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[1][3][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[1][3][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[1][3][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][3][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[1][3][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[1][3][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[1][3][3] == 0.0); // deta_dy

  // FPx 5
  // Jm
  REQUIRE(q->Jm[1][4][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[1][4][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[1][4][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[1][4][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][4][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[1][4][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[1][4][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[1][4][3] == 0.0); // deta_dy

  // FPx 6
  // Jm
  REQUIRE(q->Jm[1][5][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[1][5][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[1][5][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[1][5][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][5][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[1][5][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[1][5][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[1][5][3] == 0.0); // deta_dy

  /* *************** FP Y ****************** */
  // FPy 1
  // Jm
  REQUIRE(q->Jm[2][0][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[2][0][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[2][0][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[2][0][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][0][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[2][0][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[2][0][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[2][0][3] == 0.0); // deta_dy

  // FPy 2
  // Jm
  REQUIRE(q->Jm[2][1][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[2][1][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[2][1][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[2][1][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][1][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[2][1][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[2][1][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[2][1][3] == 0.0); // deta_dy

  // FPy 3
  // Jm
  REQUIRE(q->Jm[2][2][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[2][2][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[2][2][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[2][2][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][2][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[2][2][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[2][2][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[2][2][3] == 0.0); // deta_dy

  // FPy 4
  // Jm
  REQUIRE(q->Jm[2][3][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[2][3][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[2][3][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[2][3][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][3][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[2][3][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[2][3][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[2][3][3] == 0.0); // deta_dy

  // FPy 5
  // Jm
  REQUIRE(q->Jm[2][4][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[2][4][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[2][4][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[2][4][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][4][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[2][4][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[2][4][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[2][4][3] == 0.0); // deta_dy

  // FPy 6
  // Jm
  REQUIRE(q->Jm[2][5][0] == 0.0); // dx_dcsi
  REQUIRE(q->Jm[2][5][1] == 0.0); // dx_deta
  REQUIRE(q->Jm[2][5][2] == 0.0); // dy_dcsi
  REQUIRE(q->Jm[2][5][3] == 0.0); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][5][0] == 0.0); // dcsi_dx
  REQUIRE(q->Ji[2][5][1] == 0.0); // dcsi_dy
  REQUIRE(q->Ji[2][5][2] == 0.0); // deta_dx
  REQUIRE(q->Ji[2][5][3] == 0.0); // deta_dy
}


