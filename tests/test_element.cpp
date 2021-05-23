
#include <iostream>
#include <memory>
#include <sstream>

#include <catch/catch.hpp>
#include <Dummies.h>
#include <Element.h>
#include <Field.h>
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

  REQUIRE(q->J[0][0] == Approx(0.0002972693).margin(1E-15)); // Jacobian
  //REQUIRE(q->J[0][0] == Approx(0.0047563092).margin(1E-15)); // Jacobian

  /* *************** SP ****************** */
  // SP 1
  // Jm
  REQUIRE(q->Jm[0][0][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[0][0][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[0][0][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[0][0][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][0][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[0][0][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[0][0][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[0][0][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // SP 2
  // Jm
  REQUIRE(q->Jm[0][1][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[0][1][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[0][1][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[0][1][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][1][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[0][1][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[0][1][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[0][1][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // SP 3
  // Jm
  REQUIRE(q->Jm[0][2][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[0][2][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[0][2][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[0][2][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][2][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[0][2][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[0][2][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[0][2][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // SP 4
  // Jm
  REQUIRE(q->Jm[0][3][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[0][3][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[0][3][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[0][3][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[0][3][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[0][3][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[0][3][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[0][3][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPs
  /* *************** FP X ****************** */
  // FPx 1
  // Jm
  REQUIRE(q->Jm[1][0][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][0][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][0][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][0][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][0][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][0][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][0][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][0][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPx 2
  // Jm
  REQUIRE(q->Jm[1][1][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][1][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][1][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][1][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][1][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][1][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][1][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][1][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPx 3
  // Jm
  REQUIRE(q->Jm[1][2][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][2][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][2][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][2][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][2][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][2][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][2][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][2][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPx 4
  // Jm
  REQUIRE(q->Jm[1][3][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][3][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][3][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][3][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][3][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][3][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][3][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][3][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPx 5
  // Jm
  REQUIRE(q->Jm[1][4][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][4][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][4][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][4][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][4][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][4][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][4][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][4][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPx 6
  // Jm
  REQUIRE(q->Jm[1][5][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[1][5][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[1][5][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[1][5][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[1][5][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[1][5][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[1][5][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[1][5][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  /* *************** FP Y ****************** */
  // FPy 1
  // Jm
  REQUIRE(q->Jm[2][0][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][0][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][0][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][0][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][0][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][0][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][0][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][0][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPy 2
  // Jm
  REQUIRE(q->Jm[2][1][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][1][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][1][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][1][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][1][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][1][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][1][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][1][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPy 3
  // Jm
  REQUIRE(q->Jm[2][2][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][2][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][2][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][2][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][2][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][2][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][2][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][2][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPy 4
  // Jm
  REQUIRE(q->Jm[2][3][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][3][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][3][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][3][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][3][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][3][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][3][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][3][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPy 5
  // Jm
  REQUIRE(q->Jm[2][4][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][4][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][4][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][4][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][4][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][4][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][4][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][4][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy

  // FPy 6
  // Jm
  REQUIRE(q->Jm[2][5][0] == Approx(0.0172415).margin(1E-15)); // dx_dcsi
  REQUIRE(q->Jm[2][5][1] == 0.0);                             // dx_deta
  REQUIRE(q->Jm[2][5][2] == 0.0);                             // dy_dcsi
  REQUIRE(q->Jm[2][5][3] == Approx(0.0172415).margin(1E-15)); // dy_deta
  // Ji
  REQUIRE(q->Ji[2][5][0] == Approx(57.9995940028).margin(1E-15)); // dcsi_dx
  REQUIRE(q->Ji[2][5][1] == 0.0);                                 // dcsi_dy
  REQUIRE(q->Ji[2][5][2] == 0.0);                                 // deta_dx
  REQUIRE(q->Ji[2][5][3] == Approx(57.9995940028).margin(1E-15)); // deta_dy
}


TEST_CASE("3: Test Quadrangle - get_normals", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});
  std::vector<Vertice> enodes = {{0.0, -1.0, 0.0},
                                 {1.0, 0.0, 0.0},
                                 {0.0, 1.0, 0.0},
                                 {-1.0, 0.0, 0.0}};
  double area = 2.0;

  int order = 2;                                           // second order
  int dimension = 2;                                       // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  sd->create_nodes();
  q->allocate_jacobian(order);
  q->calculate_jacobian(sd->snodes, sd->fnodes, enodes);
  
  REQUIRE(q->Jm[0][0][0] == Approx(0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][1] == Approx(-0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][2] == Approx(0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][3] == Approx(0.5).margin(1E-15));
  REQUIRE(q->J[0][0] == Approx(area/4.0).margin(1E-15)); // A = 4*det(J) [standard square (-1, 1)]

  std::shared_ptr<Element> e = q;
  sd->initialize_element_properties(e, enodes, FIELDS::DEFAULT_FIELD_MAPPING);

  int f_index=0, dir=0, local_ed=0;
  std::vector<std::vector<double>> expected_norms = {{std::sqrt(2.0)/2.0, -std::sqrt(2.0)/2.0, 0.0},
                                                     {std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0},
                                                     {-std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0},
                                                     {-std::sqrt(2.0)/2.0, -std::sqrt(2.0)/2.0, 0.0}};
  
  std::vector<std::vector<std::vector<double>>> observed_norms = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  
  for (auto& ed : q->edges)
  {
    f_index=0;
    for (auto& node : ed.fnodes)
    {
      dir = (local_ed == 0 || local_ed == 2) ? 1 : 0; // 0: x, 1: y
      auto n = q->transform(node, enodes);
      auto norm = q->get_normal_vector(dir+1, node.local, local_ed);

      observed_norms[local_ed][f_index][0] = norm[0];
      observed_norms[local_ed][f_index][1] = norm[1];
      observed_norms[local_ed][f_index][2] = 0.0;
      f_index++;
    }
    local_ed++;
  }

  REQUIRE(observed_norms[0][0][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][0][1] == Approx(expected_norms[0][1]).margin(1E-15));
  REQUIRE(observed_norms[0][1][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][1][1] == Approx(expected_norms[0][1]).margin(1E-15));

  REQUIRE(observed_norms[1][0][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][0][1] == Approx(expected_norms[1][1]).margin(1E-15));
  REQUIRE(observed_norms[1][1][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][1][1] == Approx(expected_norms[1][1]).margin(1E-15));

  REQUIRE(observed_norms[2][0][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][0][1] == Approx(expected_norms[2][1]).margin(1E-15));
  REQUIRE(observed_norms[2][1][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][1][1] == Approx(expected_norms[2][1]).margin(1E-15));

  REQUIRE(observed_norms[3][0][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][0][1] == Approx(expected_norms[3][1]).margin(1E-15));
  REQUIRE(observed_norms[3][1][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][1][1] == Approx(expected_norms[3][1]).margin(1E-15));

}

TEST_CASE("4: Test QuadrangleHO - get_normals", "[elems]")
{
  auto q = std::make_shared<QuadrangleHO>(std::vector<std::string>{"1", "2", "3", "4"});
  std::vector<Vertice> enodes = {{0.0, -1.0, 0.0},
                                 {1.0, 0.0, 0.0},
                                 {0.0, 1.0, 0.0},
                                 {-1.0, 0.0, 0.0}};
  double area = 2.0;

  int order = 2;                                           // second order
  int dimension = 2;                                       // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  sd->create_nodes();
  q->allocate_jacobian(order);
  q->calculate_jacobian(sd->snodes, sd->fnodes, enodes);
  
  REQUIRE(q->NUM_NODES == 4);
  REQUIRE(q->NUM_VERTICES == 4);
  REQUIRE(q->ORDER == 1);
  REQUIRE(q->ce_space.size() == 2);
  REQUIRE(q->ce_space[0] == Approx(-1.0).margin(1E-15));
  REQUIRE(q->ce_space[1] == Approx(1.0).margin(1E-15));

  
  auto coords = q->computational_map.find(0)->second;
  auto i = (unsigned int) coords[0];
  auto j = (unsigned int) coords[1];
  REQUIRE(i == 0);
  REQUIRE(j == 0);

  coords = q->computational_map.find(1)->second;
  i = (unsigned int) coords[0];
  j = (unsigned int) coords[1];
  REQUIRE(i == 1);
  REQUIRE(j == 0);

  coords = q->computational_map.find(2)->second;
  i = (unsigned int) coords[0];
  j = (unsigned int) coords[1];
  REQUIRE(i == 1);
  REQUIRE(j == 1);

  coords = q->computational_map.find(3)->second;
  i = (unsigned int) coords[0];
  j = (unsigned int) coords[1];
  REQUIRE(i == 0);
  REQUIRE(j == 1);

  auto n = q->transform(Node{-1.0, -1.0, 0.0}, enodes);
  REQUIRE(n.coords[0] == Approx(0.0).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(-1.0).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(0.0).margin(1E-15));
  n = q->transform(Node{1.0, -1.0, 0.0}, enodes);
  REQUIRE(n.coords[0] == Approx(1.0).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(0.0).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(0.0).margin(1E-15));
  n = q->transform(Node{1.0, 1.0, 0.0}, enodes);
  REQUIRE(n.coords[0] == Approx(0.0).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(1.0).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(0.0).margin(1E-15));
  n = q->transform(Node{-1.0, 1.0, 0.0}, enodes);
  REQUIRE(n.coords[0] == Approx(-1.0).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(0.0).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(0.0).margin(1E-15));
  n = q->transform(Node{0.0, 0.0, 0.0}, enodes);
  REQUIRE(n.coords[0] == Approx(0.0).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(0.0).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(0.0).margin(1E-15));
  n = q->transform(Node{1.0, 0.0, 0.0}, enodes);
  REQUIRE(n.coords[0] == Approx(0.5).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(0.5).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(0.0).margin(1E-15));

  REQUIRE(q->J[0][0] == Approx(area/4.0).margin(1E-15)); // A = 4*det(J) [standard square (-1, 1)]
  REQUIRE(q->Jm[0][0][0] == Approx(0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][1] == Approx(-0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][2] == Approx(0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][3] == Approx(0.5).margin(1E-15));

  std::shared_ptr<Element> e = q;
  sd->initialize_element_properties(e, enodes, FIELDS::DEFAULT_FIELD_MAPPING);

  int f_index=0, dir=0, local_ed=0;
  std::vector<std::vector<double>> expected_norms = {{std::sqrt(2.0)/2.0, -std::sqrt(2.0)/2.0, 0.0},
                                                     {std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0},
                                                     {-std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0},
                                                     {-std::sqrt(2.0)/2.0, -std::sqrt(2.0)/2.0, 0.0}};
   std::vector<std::vector<std::vector<double>>> observed_norms = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  
  for (auto& ed : q->edges)
  {
    f_index=0;
    for (auto& node : ed.fnodes)
    {
      dir = (local_ed == 0 || local_ed == 2) ? 1 : 0; // 0: x, 1: y
      auto n = q->transform(node, enodes);
      auto norm = q->get_normal_vector(dir+1, node.local, local_ed);
      
      observed_norms[local_ed][f_index][0] = norm[0];
      observed_norms[local_ed][f_index][1] = norm[1];
      observed_norms[local_ed][f_index][2] = 0.0;
      f_index++;
    }
    local_ed++;
  }

  REQUIRE(observed_norms[0][0][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][0][1] == Approx(expected_norms[0][1]).margin(1E-15));
  REQUIRE(observed_norms[0][1][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][1][1] == Approx(expected_norms[0][1]).margin(1E-15));

  REQUIRE(observed_norms[1][0][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][0][1] == Approx(expected_norms[1][1]).margin(1E-15));
  REQUIRE(observed_norms[1][1][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][1][1] == Approx(expected_norms[1][1]).margin(1E-15));

  REQUIRE(observed_norms[2][0][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][0][1] == Approx(expected_norms[2][1]).margin(1E-15));
  REQUIRE(observed_norms[2][1][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][1][1] == Approx(expected_norms[2][1]).margin(1E-15));

  REQUIRE(observed_norms[3][0][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][0][1] == Approx(expected_norms[3][1]).margin(1E-15));
  REQUIRE(observed_norms[3][1][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][1][1] == Approx(expected_norms[3][1]).margin(1E-15));
}


TEST_CASE("5: Test Quadrangle - get_normals", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});
  std::vector<Vertice> enodes = {{0.0, 1.0, 0.0},
                                 {-1.0, 0.0, 0.0},
                                 {0.0, -1.0, 0.0},
                                 {1.0, 0.0, 0.0}};
  double area = 2.0;

  int order = 2;                                           // second order
  int dimension = 2;                                       // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  sd->create_nodes();
  q->allocate_jacobian(order);
  q->calculate_jacobian(sd->snodes, sd->fnodes, enodes);
  
  REQUIRE(q->Jm[0][0][0] == Approx(-0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][1] == Approx(0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][2] == Approx(-0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][3] == Approx(-0.5).margin(1E-15));
  REQUIRE(q->J[0][0] == Approx(area/4.0).margin(1E-15)); // A = 4*det(J) [standard square (-1, 1)]

  std::shared_ptr<Element> e = q;
  sd->initialize_element_properties(e, enodes, FIELDS::DEFAULT_FIELD_MAPPING);

  int f_index=0, dir=0, local_ed=0;
  std::vector<std::vector<double>> expected_norms = {{-std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0},
                                                     {-std::sqrt(2.0)/2.0, -std::sqrt(2.0)/2.0, 0.0},
                                                     {std::sqrt(2.0)/2.0, -std::sqrt(2.0)/2.0, 0.0},
                                                     {std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0}};
  
  std::vector<std::vector<std::vector<double>>> observed_norms = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  
  for (auto& ed : q->edges)
  {
    f_index=0;
    for (auto& node : ed.fnodes)
    {
      dir = (local_ed == 0 || local_ed == 2) ? 1 : 0; // 0: x, 1: y
      auto n = q->transform(node, enodes);
      auto norm = q->get_normal_vector(dir+1, node.local, local_ed);

      observed_norms[local_ed][f_index][0] = norm[0];
      observed_norms[local_ed][f_index][1] = norm[1];
      observed_norms[local_ed][f_index][2] = 0.0;
      f_index++;
    }
    local_ed++;
  }

  REQUIRE(observed_norms[0][0][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][0][1] == Approx(expected_norms[0][1]).margin(1E-15));
  REQUIRE(observed_norms[0][1][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][1][1] == Approx(expected_norms[0][1]).margin(1E-15));

  REQUIRE(observed_norms[1][0][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][0][1] == Approx(expected_norms[1][1]).margin(1E-15));
  REQUIRE(observed_norms[1][1][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][1][1] == Approx(expected_norms[1][1]).margin(1E-15));

  REQUIRE(observed_norms[2][0][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][0][1] == Approx(expected_norms[2][1]).margin(1E-15));
  REQUIRE(observed_norms[2][1][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][1][1] == Approx(expected_norms[2][1]).margin(1E-15));

  REQUIRE(observed_norms[3][0][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][0][1] == Approx(expected_norms[3][1]).margin(1E-15));
  REQUIRE(observed_norms[3][1][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][1][1] == Approx(expected_norms[3][1]).margin(1E-15));

}

TEST_CASE("6: Test Quadrangle - get_normals", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});
  std::vector<Vertice> enodes = {{-1.0, 0.0, 0.0},
                                 {0.0, -1.0, 0.0},
                                 {1.0, 0.0, 0.0},
                                 {0.0, 1.0, 0.0}};
  double area = 2.0;

  int order = 2;                                           // second order
  int dimension = 2;                                       // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  sd->create_nodes();
  q->allocate_jacobian(order);
  q->calculate_jacobian(sd->snodes, sd->fnodes, enodes);
  
  REQUIRE(q->Jm[0][0][0] == Approx(0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][1] == Approx(0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][2] == Approx(-0.5).margin(1E-15));
  REQUIRE(q->Jm[0][0][3] == Approx(0.5).margin(1E-15));
  REQUIRE(q->J[0][0] == Approx(area/4.0).margin(1E-15)); // A = 4*det(J) [standard square (-1, 1)]

  std::shared_ptr<Element> e = q;
  sd->initialize_element_properties(e, enodes, FIELDS::DEFAULT_FIELD_MAPPING);

  int f_index=0, dir=0, local_ed=0;
  std::vector<std::vector<double>> expected_norms = {{-std::sqrt(2.0)/2.0, -std::sqrt(2.0)/2.0, 0.0},
                                                     {std::sqrt(2.0)/2.0, -std::sqrt(2.0)/2.0, 0.0},
                                                     {std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0},
                                                     {-std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0}};
  
  std::vector<std::vector<std::vector<double>>> observed_norms = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  
  for (auto& ed : q->edges)
  {
    f_index=0;
    for (auto& node : ed.fnodes)
    {
      dir = (local_ed == 0 || local_ed == 2) ? 1 : 0; // 0: x, 1: y
      auto n = q->transform(node, enodes);
      auto norm = q->get_normal_vector(dir+1, node.local, local_ed);

      observed_norms[local_ed][f_index][0] = norm[0];
      observed_norms[local_ed][f_index][1] = norm[1];
      observed_norms[local_ed][f_index][2] = 0.0;
      f_index++;
    }
    local_ed++;
  }

  REQUIRE(observed_norms[0][0][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][0][1] == Approx(expected_norms[0][1]).margin(1E-15));
  REQUIRE(observed_norms[0][1][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][1][1] == Approx(expected_norms[0][1]).margin(1E-15));

  REQUIRE(observed_norms[1][0][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][0][1] == Approx(expected_norms[1][1]).margin(1E-15));
  REQUIRE(observed_norms[1][1][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][1][1] == Approx(expected_norms[1][1]).margin(1E-15));

  REQUIRE(observed_norms[2][0][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][0][1] == Approx(expected_norms[2][1]).margin(1E-15));
  REQUIRE(observed_norms[2][1][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][1][1] == Approx(expected_norms[2][1]).margin(1E-15));

  REQUIRE(observed_norms[3][0][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][0][1] == Approx(expected_norms[3][1]).margin(1E-15));
  REQUIRE(observed_norms[3][1][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][1][1] == Approx(expected_norms[3][1]).margin(1E-15));

}


TEST_CASE("7: Test Quadrangle - get_normals", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});
  std::vector<Vertice> enodes = {{0.5, 0.5, 0.0},
                                 {0.75, 0.25, 0.0},
                                 {1.0, 0.6, 0.0},
                                 {0.8, 1.0, 0.0}};
  double area = 2.0;

  int order = 2;                                           // second order
  int dimension = 2;                                       // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  sd->create_nodes();
  q->allocate_jacobian(order);
  q->calculate_jacobian(sd->snodes, sd->fnodes, enodes);
  
  std::shared_ptr<Element> e = q;
  sd->initialize_element_properties(e, enodes, FIELDS::DEFAULT_FIELD_MAPPING);

  int f_index=0, dir=0, local_ed=0;
  
  std::vector<std::vector<double>> expected_norms = {{-0.7071067811865474, -0.7071067811865475, 0.0},
                                                     {0.8137334712067351, -0.5812381937190967, 0.0},
                                                     {0.8944271909999157, 0.4472135954999581, 0.0},
                                                     {-0.8574929257125441, 0.5144957554275265, 0.0}};
  
  std::vector<std::vector<std::vector<double>>> observed_norms = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  auto i_exp = 0;
  for (auto en : expected_norms)
  {
    auto v = expected_norms[i_exp];
    auto norm = std::sqrt(v[0]*v[0] + v[1]*v[1]);

    expected_norms[i_exp][0] = v[0]/norm;
    expected_norms[i_exp][1] = v[1]/norm;
    expected_norms[i_exp][2] = 0.0;

    i_exp++;
  }

  for (auto& ed : q->edges)
  {
    f_index=0;
    for (auto& node : ed.fnodes)
    {
      dir = (local_ed == 0 || local_ed == 2) ? 1 : 0; // 0: x, 1: y
      auto norm = q->get_normal_vector(dir+1, node.local, local_ed);

      observed_norms[local_ed][f_index][0] = norm[0];
      observed_norms[local_ed][f_index][1] = norm[1];
      observed_norms[local_ed][f_index][2] = 0.0;

      std::cout << "local_ed = " << local_ed << "\n";
      std::cout << "fnode = " << node.coords[0] << "; " << node.coords[1] << "\n";
      std::cout << observed_norms[local_ed][f_index][0] << "; " << observed_norms[local_ed][f_index][1] << "\n";
      
      f_index++;
    }
    local_ed++;
  }

  REQUIRE(observed_norms[0][0][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][0][1] == Approx(expected_norms[0][1]).margin(1E-15));
  REQUIRE(observed_norms[0][1][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][1][1] == Approx(expected_norms[0][1]).margin(1E-15));

  REQUIRE(observed_norms[1][0][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][0][1] == Approx(expected_norms[1][1]).margin(1E-15));
  REQUIRE(observed_norms[1][1][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][1][1] == Approx(expected_norms[1][1]).margin(1E-15));

  REQUIRE(observed_norms[2][0][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][0][1] == Approx(expected_norms[2][1]).margin(1E-15));
  REQUIRE(observed_norms[2][1][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][1][1] == Approx(expected_norms[2][1]).margin(1E-15));

  REQUIRE(observed_norms[3][0][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][0][1] == Approx(expected_norms[3][1]).margin(1E-15));
  REQUIRE(observed_norms[3][1][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][1][1] == Approx(expected_norms[3][1]).margin(1E-15));

}

TEST_CASE("8: Test Quadrangle - get_normals", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});
  std::vector<Vertice> enodes = {{-0.5, 0.5, 0.0},
                                 {-0.25, 0.75, 0.0},
                                 {-0.6, 1.0, 0.0},
                                 {-1.0, 0.8, 0.0}};
  double area = 2.0;

  int order = 2;                                           // second order
  int dimension = 2;                                       // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  sd->create_nodes();
  q->allocate_jacobian(order);
  q->calculate_jacobian(sd->snodes, sd->fnodes, enodes);
  
  std::shared_ptr<Element> e = q;
  sd->initialize_element_properties(e, enodes, FIELDS::DEFAULT_FIELD_MAPPING);

  int f_index=0, dir=0, local_ed=0;
  
  std::vector<std::vector<double>> expected_norms = {{0.7071067811865475, -0.7071067811865472,0.0},
                                                     {0.5812381937190965, 0.813733471206735,0.0},
                                                     {-0.44721359549995787, 0.894427190999916,0.0},
                                                     {-0.5144957554275267, -0.8574929257125442,0.0}};
  
  std::vector<std::vector<std::vector<double>>> observed_norms = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  auto i_exp = 0;
  for (auto en : expected_norms)
  {
    auto v = expected_norms[i_exp];
    auto norm = std::sqrt(v[0]*v[0] + v[1]*v[1]);

    expected_norms[i_exp][0] = v[0]/norm;
    expected_norms[i_exp][1] = v[1]/norm;
    expected_norms[i_exp][2] = 0.0;

    i_exp++;
  }

  for (auto& ed : q->edges)
  {
    f_index=0;
    for (auto& node : ed.fnodes)
    {
      dir = (local_ed == 0 || local_ed == 2) ? 1 : 0; // 0: x, 1: y
      auto n = q->transform(node, enodes);
      auto norm = q->get_normal_vector(dir+1, node.local, local_ed);

      observed_norms[local_ed][f_index][0] = norm[0];
      observed_norms[local_ed][f_index][1] = norm[1];
      observed_norms[local_ed][f_index][2] = 0.0;
      f_index++;
    }
    local_ed++;
  }

  REQUIRE(observed_norms[0][0][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][0][1] == Approx(expected_norms[0][1]).margin(1E-15));
  REQUIRE(observed_norms[0][1][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][1][1] == Approx(expected_norms[0][1]).margin(1E-15));

  REQUIRE(observed_norms[1][0][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][0][1] == Approx(expected_norms[1][1]).margin(1E-15));
  REQUIRE(observed_norms[1][1][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][1][1] == Approx(expected_norms[1][1]).margin(1E-15));

  REQUIRE(observed_norms[2][0][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][0][1] == Approx(expected_norms[2][1]).margin(1E-15));
  REQUIRE(observed_norms[2][1][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][1][1] == Approx(expected_norms[2][1]).margin(1E-15));

  REQUIRE(observed_norms[3][0][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][0][1] == Approx(expected_norms[3][1]).margin(1E-15));
  REQUIRE(observed_norms[3][1][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][1][1] == Approx(expected_norms[3][1]).margin(1E-15));

}


TEST_CASE("9: Test Quadrangle - get_normals", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});
  std::vector<Vertice> enodes = {{-0.5, -0.5, 0.0},
                                 {-0.75, -0.25, 0.0},
                                 {-1.0, -0.6, 0.0},
                                 {-0.8, -1.0, 0.0}};
  double area = 2.0;

  int order = 2;                                           // second order
  int dimension = 2;                                       // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  sd->create_nodes();
  q->allocate_jacobian(order);
  q->calculate_jacobian(sd->snodes, sd->fnodes, enodes);
  
  std::shared_ptr<Element> e = q;
  sd->initialize_element_properties(e, enodes, FIELDS::DEFAULT_FIELD_MAPPING);

  int f_index=0, dir=0, local_ed=0;

  std::vector<std::vector<double>> expected_norms = {{0.7071067811865475, 0.7071067811865475, 0.0},
                                                     {-0.813733471206735, 0.5812381937190965, 0.0},
                                                     {-0.8944271909999162, -0.44721359549995765, 0.0},
                                                     {0.8574929257125439, -0.5144957554275268, 0.0}};
  
  std::vector<std::vector<std::vector<double>>> observed_norms = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  auto i_exp = 0;
  for (auto en : expected_norms)
  {
    auto v = expected_norms[i_exp];
    auto norm = std::sqrt(v[0]*v[0] + v[1]*v[1]);

    expected_norms[i_exp][0] = v[0]/norm;
    expected_norms[i_exp][1] = v[1]/norm;
    expected_norms[i_exp][2] = 0.0;

    i_exp++;
  }

  for (auto& ed : q->edges)
  {
    f_index=0;
    for (auto& node : ed.fnodes)
    {
      dir = (local_ed == 0 || local_ed == 2) ? 1 : 0; // 0: x, 1: y
      auto n = q->transform(node, enodes);
      auto norm = q->get_normal_vector(dir+1, node.local, local_ed);

      observed_norms[local_ed][f_index][0] = norm[0];
      observed_norms[local_ed][f_index][1] = norm[1];
      observed_norms[local_ed][f_index][2] = 0.0;
      f_index++;
    }
    local_ed++;
  }

  REQUIRE(observed_norms[0][0][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][0][1] == Approx(expected_norms[0][1]).margin(1E-15));
  REQUIRE(observed_norms[0][1][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][1][1] == Approx(expected_norms[0][1]).margin(1E-15));

  REQUIRE(observed_norms[1][0][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][0][1] == Approx(expected_norms[1][1]).margin(1E-15));
  REQUIRE(observed_norms[1][1][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][1][1] == Approx(expected_norms[1][1]).margin(1E-15));

  REQUIRE(observed_norms[2][0][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][0][1] == Approx(expected_norms[2][1]).margin(1E-15));
  REQUIRE(observed_norms[2][1][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][1][1] == Approx(expected_norms[2][1]).margin(1E-15));

  REQUIRE(observed_norms[3][0][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][0][1] == Approx(expected_norms[3][1]).margin(1E-15));
  REQUIRE(observed_norms[3][1][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][1][1] == Approx(expected_norms[3][1]).margin(1E-15));

}

TEST_CASE("10: Test Quadrangle - get_normals", "[elems]")
{
  auto q = std::make_shared<Quadrangle>(std::vector<std::string>{"1", "2", "3", "4"});
  std::vector<Vertice> enodes = {{0.5, -0.5, 0.0},
                                 {0.25, -0.75, 0.0},
                                 {0.6, -1.0, 0.0},
                                 {1.0, -0.8, 0.0}};
  double area = 2.0;

  int order = 2;                                           // second order
  int dimension = 2;                                       // 2D
  auto sd = std::make_shared<SD<Euler>>(order, dimension); // (int order, int dimension)
  sd->create_nodes();
  q->allocate_jacobian(order);
  q->calculate_jacobian(sd->snodes, sd->fnodes, enodes);
  
  std::shared_ptr<Element> e = q;
  sd->initialize_element_properties(e, enodes, FIELDS::DEFAULT_FIELD_MAPPING);

  int f_index=0, dir=0, local_ed=0;

  std::vector<std::vector<double>> expected_norms = {{-0.7071067811865475, 0.7071067811865474, 0.0},
                                                     {-0.5812381937190965, -0.8137334712067349, 0.0},
                                                     {0.4472135954999572, -0.8944271909999162, 0.0},
                                                     {0.5144957554275269, 0.857492925712544, 0.0}};
  
  std::vector<std::vector<std::vector<double>>> observed_norms = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  auto i_exp = 0;
  for (auto en : expected_norms)
  {
    auto v = expected_norms[i_exp];
    auto norm = std::sqrt(v[0]*v[0] + v[1]*v[1]);

    expected_norms[i_exp][0] = v[0]/norm;
    expected_norms[i_exp][1] = v[1]/norm;
    expected_norms[i_exp][2] = 0.0;

    i_exp++;
  }

  for (auto& ed : q->edges)
  {
    f_index=0;
    for (auto& node : ed.fnodes)
    {
      dir = (local_ed == 0 || local_ed == 2) ? 1 : 0; // 0: x, 1: y
      auto n = q->transform(node, enodes);
      auto norm = q->get_normal_vector(dir+1, node.local, local_ed);

      observed_norms[local_ed][f_index][0] = norm[0];
      observed_norms[local_ed][f_index][1] = norm[1];
      observed_norms[local_ed][f_index][2] = 0.0;
      f_index++;
    }
    local_ed++;
  }

  REQUIRE(observed_norms[0][0][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][0][1] == Approx(expected_norms[0][1]).margin(1E-15));
  REQUIRE(observed_norms[0][1][0] == Approx(expected_norms[0][0]).margin(1E-15));
  REQUIRE(observed_norms[0][1][1] == Approx(expected_norms[0][1]).margin(1E-15));

  REQUIRE(observed_norms[1][0][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][0][1] == Approx(expected_norms[1][1]).margin(1E-15));
  REQUIRE(observed_norms[1][1][0] == Approx(expected_norms[1][0]).margin(1E-15));
  REQUIRE(observed_norms[1][1][1] == Approx(expected_norms[1][1]).margin(1E-15));

  REQUIRE(observed_norms[2][0][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][0][1] == Approx(expected_norms[2][1]).margin(1E-15));
  REQUIRE(observed_norms[2][1][0] == Approx(expected_norms[2][0]).margin(1E-15));
  REQUIRE(observed_norms[2][1][1] == Approx(expected_norms[2][1]).margin(1E-15));

  REQUIRE(observed_norms[3][0][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][0][1] == Approx(expected_norms[3][1]).margin(1E-15));
  REQUIRE(observed_norms[3][1][0] == Approx(expected_norms[3][0]).margin(1E-15));
  REQUIRE(observed_norms[3][1][1] == Approx(expected_norms[3][1]).margin(1E-15));

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
