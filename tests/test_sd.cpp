#include <filesystem>
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
      namespace fs = std::filesystem;

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
  }

  for (auto& vec: sd->fnodes)
  { 
    for (auto& fn: vec)
    {
      REQUIRE(fn.id == -1);
      REQUIRE(fn.coords[0] == 0.0);
      REQUIRE(fn.coords[1] == 0.0);
      REQUIRE(fn.coords[2] == 0.0);
    }
  }
}


TEST_CASE("2: Test SD<Euler> - create_nodes", "[sd2]")
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
  
  // 1
  REQUIRE(sd->snodes[1].id == -1);
  REQUIRE(sd->snodes[1].coords[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[1].coords[1] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[1].coords[2] == Approx(0.0).margin(1E-15));
  
  // 2
  REQUIRE(sd->snodes[2].id == -1);
  REQUIRE(sd->snodes[2].coords[0] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[2].coords[1] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[2].coords[2] == Approx(0.0).margin(1E-15));
  
  // 3
  REQUIRE(sd->snodes[3].id == -1);
  REQUIRE(sd->snodes[3].coords[0] == Approx(0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[3].coords[1] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->snodes[3].coords[2] == Approx(0.0).margin(1E-15));
  
  // FNODES - x
  // 0
  REQUIRE(sd->fnodes[0][0].id == -1);
  REQUIRE(sd->fnodes[0][0].coords[0] == -1.0);
  REQUIRE(sd->fnodes[0][0].coords[1] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][0].coords[2] == Approx(0.0).margin(1E-15));
  
  // 1
  REQUIRE(sd->fnodes[0][1].id == -1);
  REQUIRE(sd->fnodes[0][1].coords[0] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[0][1].coords[1] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][1].coords[2] == Approx(0.0).margin(1E-15));
  
  // 2
  REQUIRE(sd->fnodes[0][2].id == -1);
  REQUIRE(sd->fnodes[0][2].coords[0] == 1.0);
  REQUIRE(sd->fnodes[0][2].coords[1] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][2].coords[2] == Approx(0.0).margin(1E-15));
  
  // 3
  REQUIRE(sd->fnodes[0][3].id == -1);
  REQUIRE(sd->fnodes[0][3].coords[0] == -1.0);
  REQUIRE(sd->fnodes[0][3].coords[1] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][3].coords[2] == Approx(0.0).margin(1E-15));
  
  // 4
  REQUIRE(sd->fnodes[0][4].id == -1);
  REQUIRE(sd->fnodes[0][4].coords[0] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[0][4].coords[1] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][4].coords[2] == Approx(0.0).margin(1E-15));
  
  // 5
  REQUIRE(sd->fnodes[0][5].id == -1);
  REQUIRE(sd->fnodes[0][5].coords[0] == 1.0);
  REQUIRE(sd->fnodes[0][5].coords[1] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[0][5].coords[2] == Approx(0.0).margin(1E-15));
  
  // FNODES - y
  // 0
  REQUIRE(sd->fnodes[1][0].id == -1);
  REQUIRE(sd->fnodes[1][0].coords[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][0].coords[1] == -1.0);
  REQUIRE(sd->fnodes[1][0].coords[2] == Approx(0.0).margin(1E-15));
  
  // 1
  REQUIRE(sd->fnodes[1][1].id == -1);
  REQUIRE(sd->fnodes[1][1].coords[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][1].coords[1] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[1][1].coords[2] == Approx(0.0).margin(1E-15));
  
  // 2
  REQUIRE(sd->fnodes[1][2].id == -1);
  REQUIRE(sd->fnodes[1][2].coords[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][2].coords[1] == 1.0);
  REQUIRE(sd->fnodes[1][2].coords[2] == Approx(0.0).margin(1E-15));
  
  // 3
  REQUIRE(sd->fnodes[1][3].id == -1);
  REQUIRE(sd->fnodes[1][3].coords[0] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][3].coords[1] == -1.0);
  REQUIRE(sd->fnodes[1][3].coords[2] == Approx(0.0).margin(1E-15));
  
  // 4
  REQUIRE(sd->fnodes[1][4].id == -1);
  REQUIRE(sd->fnodes[1][4].coords[0] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][4].coords[1] == Approx(0.0).margin(1E-15));
  REQUIRE(sd->fnodes[1][4].coords[2] == Approx(0.0).margin(1E-15));
  
  // 5
  REQUIRE(sd->fnodes[1][5].id == -1);
  REQUIRE(sd->fnodes[1][5].coords[0] == Approx(+0.577350269189626).margin(1E-15));
  REQUIRE(sd->fnodes[1][5].coords[1] == 1.0);
  REQUIRE(sd->fnodes[1][5].coords[2] == Approx(0.0).margin(1E-15));
}

TEST_CASE("3: Test SD<Euler> - setup", "[sd3gaussian]")
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "mesh_test_gaussian.msh").string());

  int order = 2;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  
  sd->setup(mesh, FIELDS::GAUSSIAN_FIELD_MAPPING);
  
  // for (auto& e : mesh->elems) 
  // {
  //   std::cout << "Element: " << e->id << std::endl;
  //   long index = 0;
  //   for (auto& sol : e->physical->Qsp) 
  //   {
  //     std::cout << "  Node   : " << index << std::endl;
  //     std::cout << "    pho  : " << sol[0] << std::endl;
  //     std::cout << "    pho*u: " << sol[1] << std::endl;
  //     std::cout << "    pho*v: " << sol[2] << std::endl;
  //     std::cout << "    E    : " << sol[3] << std::endl;
  //     index++;
  //   }
  //   std::cin.get();
  // }

  sd->to_vtk(mesh, 
            (cur_path.parent_path() / "results" / "pp_mesh_test_gaussian.vtk").string());
}


TEST_CASE("4: Test SD<Euler> - setup", "[sd3uniform]")
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "mesh_test_0_order.msh").string());

  int order = 2;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  
  sd->setup(mesh, FIELDS::DEFAULT_FIELD_MAPPING);

  // Element 0:
  auto e = mesh->elems[0];
  // snodes
  // sp 0
  auto n = e->transform(sd->snodes[0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.3943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 1
  n = e->transform(sd->snodes[1], mesh->nodes);
  

  REQUIRE(n.coords[0] == Approx(+0.1056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 2
  n = e->transform(sd->snodes[2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.3943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 3
  n = e->transform(sd->snodes[3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.1056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  
  // x-fnodes
  // fp 0
  n = e->transform(sd->fnodes[0][0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.3943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));

  // fp 1
  n = e->transform(sd->fnodes[0][1], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.3943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 2
  n = e->transform(sd->fnodes[0][2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.3943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 3
  n = e->transform(sd->fnodes[0][3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.1056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 4
  n = e->transform(sd->fnodes[0][4], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.1056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 5
  n = e->transform(sd->fnodes[0][5], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.1056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));

  // y-fnodes
  // fp 0
  n = e->transform(sd->fnodes[1][0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 1
  n = e->transform(sd->fnodes[1][1], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.250).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 2
  n = e->transform(sd->fnodes[1][2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.000).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 3
  n = e->transform(sd->fnodes[1][3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 4
  n = e->transform(sd->fnodes[1][4], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.250).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 5
  n = e->transform(sd->fnodes[1][5], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.000).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));

  // Neighboors
  // edge 0
  auto& ed = e->edges[0];
  REQUIRE(ed.nodes[0] == 8);
  REQUIRE(ed.nodes[1] == 6);
  REQUIRE(ed.ghost == -1);
  REQUIRE(ed.left == 0);
  REQUIRE(ed.right == 1);
  REQUIRE(ed.lr_edge == 3);
  REQUIRE(ed.boundary == 0);
  auto& fn = ed.fnodes[0];
  REQUIRE(fn.local == 0);
  REQUIRE(fn.right == 0);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 3);
  REQUIRE(fn.right == 3);
  // Normal
  auto norm = e->get_normal_vector(2, 0, 0); // y-fnode[0] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));
  
  norm = e->get_normal_vector(2, 3, 0);      // y-fnode[3] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));

  // edge 1
  ed = e->edges[1];
  REQUIRE(ed.nodes[0] == 6);
  REQUIRE(ed.nodes[1] == 3);
  REQUIRE(ed.ghost == 0);
  REQUIRE(ed.left == 0);
  REQUIRE(ed.right == -1);
  REQUIRE(ed.lr_edge == -1);
  REQUIRE(ed.boundary == 1);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 2);
  REQUIRE(fn.right == 2);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 5);
  REQUIRE(fn.right == 5);
  // Normal
  norm = e->get_normal_vector(1, 2, 1); // x-fnode[2] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+1.0).margin(1E-15));
  
  norm = e->get_normal_vector(1, 5, 1);      // x-fnode[5] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+1.0).margin(1E-15));
  
  // edge 2
  ed = e->edges[2];
  REQUIRE(ed.nodes[0] == 3);
  REQUIRE(ed.nodes[1] == 7);
  REQUIRE(ed.ghost == 1);
  REQUIRE(ed.left == 0);
  REQUIRE(ed.right == -1);
  REQUIRE(ed.lr_edge == -1);
  REQUIRE(ed.boundary == 1);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 5);
  REQUIRE(fn.right == 5);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 2);
  REQUIRE(fn.right == 2);
  // Normal
  norm = e->get_normal_vector(2, 5, 2); // y-fnode[5] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(-1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));
  
  norm = e->get_normal_vector(2, 2, 2);      // y-fnode[2] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(-1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));

  // edge 3
  ed = e->edges[3];
  REQUIRE(ed.nodes[0] == 7);
  REQUIRE(ed.nodes[1] == 8);
  REQUIRE(ed.ghost == -1);
  REQUIRE(ed.left == 0);
  REQUIRE(ed.right == 3);
  REQUIRE(ed.lr_edge == 3);
  REQUIRE(ed.boundary == 0);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 3);
  REQUIRE(fn.right == 0);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 0);
  REQUIRE(fn.right == 3);
  // Normal
  norm = e->get_normal_vector(1, 3, 3); // x-fnode[3] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(-1.0).margin(1E-15));
  
  norm = e->get_normal_vector(1, 0, 3);      // x-fnode[0] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(-1.0).margin(1E-15));

  std::cout << "Element 1)" << "\n";

  // Element 1:
  e = mesh->elems[1];
  // snodes
  // sp 0
  n = e->transform(sd->snodes[0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 1
  n = e->transform(sd->snodes[1], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 2
  n = e->transform(sd->snodes[2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 3
  n = e->transform(sd->snodes[3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  
  // x-fnodes
  // fp 0
  n = e->transform(sd->fnodes[0][0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 1
  n = e->transform(sd->fnodes[0][1], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 2
  n = e->transform(sd->fnodes[0][2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 3
  n = e->transform(sd->fnodes[0][3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 4
  n = e->transform(sd->fnodes[0][4], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 5
  n = e->transform(sd->fnodes[0][5], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));

  // y-fnodes
  // fp 0
  n = e->transform(sd->fnodes[1][0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 1
  n = e->transform(sd->fnodes[1][1], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 2
  n = e->transform(sd->fnodes[1][2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 3
  n = e->transform(sd->fnodes[1][3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 4
  n = e->transform(sd->fnodes[1][4], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 5
  n = e->transform(sd->fnodes[1][5], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));

  // Neighboors
  // edge 0
  ed = e->edges[0];
  REQUIRE(ed.nodes[0] == 8);
  REQUIRE(ed.nodes[1] == 5);
  REQUIRE(ed.ghost == -1);
  REQUIRE(ed.left == 1);
  REQUIRE(ed.right == 2);
  REQUIRE(ed.lr_edge == 3);
  REQUIRE(ed.boundary == 0);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 0);
  REQUIRE(fn.right == 0);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 3);
  REQUIRE(fn.right == 3);
  // Normal
  norm = e->get_normal_vector(2, 0, 0); // y-fnode[0] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(-1.0).margin(1E-15));
  
  norm = e->get_normal_vector(2, 3, 0);      // y-fnode[3] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(-1.0).margin(1E-15));

  // edge 1
  ed = e->edges[1];
  REQUIRE(ed.nodes[0] == 5);
  REQUIRE(ed.nodes[1] == 2);
  REQUIRE(ed.ghost == 2);
  REQUIRE(ed.left == 1);
  REQUIRE(ed.right == -1);
  REQUIRE(ed.lr_edge == -1);
  REQUIRE(ed.boundary == 1);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 2);
  REQUIRE(fn.right == 2);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 5);
  REQUIRE(fn.right == 5);
  // Normal
  norm = e->get_normal_vector(1, 2, 1); // x-fnode[2] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));
  
  norm = e->get_normal_vector(1, 5, 1);      // x-fnode[5] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));

  // edge 2
  ed = e->edges[2];
  REQUIRE(ed.nodes[0] == 2);
  REQUIRE(ed.nodes[1] == 6);
  REQUIRE(ed.ghost == 3);
  REQUIRE(ed.left == 1);
  REQUIRE(ed.right == -1);
  REQUIRE(ed.lr_edge == -1);
  REQUIRE(ed.boundary == 1);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 5);
  REQUIRE(fn.right == 5);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 2);
  REQUIRE(fn.right == 2);
  // Normal
  norm = e->get_normal_vector(2, 5, 2); // y-fnode[5] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+1.0).margin(1E-15));
  
  norm = e->get_normal_vector(2, 2, 2);      // y-fnode[2] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+1.0).margin(1E-15));

  // edge 3
  ed = e->edges[3];
  REQUIRE(ed.nodes[0] == 6);
  REQUIRE(ed.nodes[1] == 8);
  REQUIRE(ed.ghost == -1);
  REQUIRE(ed.left == 1);
  REQUIRE(ed.right == 0);
  REQUIRE(ed.lr_edge == 0);
  REQUIRE(ed.boundary == 0);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 3);
  REQUIRE(fn.right == 3);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 0);
  REQUIRE(fn.right == 0);
  // Normal
  norm = e->get_normal_vector(1, 3, 3); // x-fnode[3] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(-1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));
  
  norm = e->get_normal_vector(1, 0, 3);      // x-fnode[0] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(-1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));



  sd->to_vtk(mesh, 
            (cur_path.parent_path() / "results" / "pp_mesh_test_0_order.vtk").string());
}

TEST_CASE("5: Test SD<Euler> - setup", "[sd3linear]")
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "mesh_test_linear.msh").string());

  int order = 2;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  
  sd->setup(
    mesh, 
    FIELDS::LINEAR_FIELD_MAPPING
  );

  sd->to_vtk(mesh, 
            (cur_path.parent_path() / "results" / "pp_mesh_test_linear.vtk").string());
}

TEST_CASE("6: Test SD<Euler> - setup", "[sd3uniform2]")
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "mesh_test_incl.msh").string());

  int order = 2;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  
  sd->setup(mesh, FIELDS::DEFAULT_FIELD_MAPPING);

  // Element 0:
  auto e = mesh->elems[0];
  // Neighboors
  // edge 0
  auto& ed = e->edges[0];
  // REQUIRE(ed.nodes[0] == 0);
  // REQUIRE(ed.nodes[1] == 4);
  // REQUIRE(ed.ghost == -1);
  // REQUIRE(ed.left == 0);
  // REQUIRE(ed.right == 1);
  // REQUIRE(ed.lr_edge == 3);
  // REQUIRE(ed.boundary == 0);
  // auto& fn = ed.fnodes[0];
  // REQUIRE(fn.local == 0);
  // REQUIRE(fn.right == 0);
  // fn = ed.fnodes[1];
  // REQUIRE(fn.local == 3);
  // REQUIRE(fn.right == 3);
  // Normal
  auto norm = e->get_normal_vector(2, 0, 0); // y-fnode[0] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  
  // Before
  double old_nx = 1.0;
  double old_ny = 0.0;
  // After rotation
  double theta = -M_PI/4;
  double nx = cos(theta)*old_nx - sin(theta)*old_ny;
  double ny = sin(theta)*old_nx + cos(theta)*old_ny;
  
  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));
  
  norm = e->get_normal_vector(2, 3, 0);      // y-fnode[3] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));

  // edge 1
  ed = e->edges[1];
  // REQUIRE(ed.nodes[0] == 6);
  // REQUIRE(ed.nodes[1] == 3);
  // REQUIRE(ed.ghost == 0);
  // REQUIRE(ed.left == 0);
  // REQUIRE(ed.right == -1);
  // REQUIRE(ed.lr_edge == -1);
  // REQUIRE(ed.boundary == 1);
  // fn = ed.fnodes[0];
  // REQUIRE(fn.local == 2);
  // REQUIRE(fn.right == 2);
  // fn = ed.fnodes[1];
  // REQUIRE(fn.local == 5);
  // REQUIRE(fn.right == 5);
  // Normal
  norm = e->get_normal_vector(1, 2, 1); // x-fnode[2] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  // Before
  old_nx = 0.0;
  old_ny = 1.0;
  // After rotation
  theta = -M_PI/4;
  nx = cos(theta)*old_nx - sin(theta)*old_ny;
  ny = sin(theta)*old_nx + cos(theta)*old_ny;

  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));
  
  norm = e->get_normal_vector(1, 5, 1);      // x-fnode[5] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));
  
  // edge 2
  ed = e->edges[2];
  // REQUIRE(ed.nodes[0] == 3);
  // REQUIRE(ed.nodes[1] == 7);
  // REQUIRE(ed.ghost == 1);
  // REQUIRE(ed.left == 0);
  // REQUIRE(ed.right == -1);
  // REQUIRE(ed.lr_edge == -1);
  // REQUIRE(ed.boundary == 1);
  // fn = ed.fnodes[0];
  // REQUIRE(fn.local == 5);
  // REQUIRE(fn.right == 5);
  // fn = ed.fnodes[1];
  // REQUIRE(fn.local == 2);
  // REQUIRE(fn.right == 2);
  // Normal
  norm = e->get_normal_vector(2, 5, 2); // y-fnode[5] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  // Before
  old_nx = -1.0;
  old_ny = 0.0;
  // After rotation
  theta = -M_PI/4;
  nx = cos(theta)*old_nx - sin(theta)*old_ny;
  ny = sin(theta)*old_nx + cos(theta)*old_ny;

  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));
  
  norm = e->get_normal_vector(2, 2, 2);      // y-fnode[2] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));

  // edge 3
  ed = e->edges[3];
  // REQUIRE(ed.nodes[0] == 7);
  // REQUIRE(ed.nodes[1] == 8);
  // REQUIRE(ed.ghost == -1);
  // REQUIRE(ed.left == 0);
  // REQUIRE(ed.right == 3);
  // REQUIRE(ed.lr_edge == 3);
  // REQUIRE(ed.boundary == 0);
  // fn = ed.fnodes[0];
  // REQUIRE(fn.local == 3);
  // REQUIRE(fn.right == 0);
  // fn = ed.fnodes[1];
  // REQUIRE(fn.local == 0);
  // REQUIRE(fn.right == 3);
  // Normal
  norm = e->get_normal_vector(1, 3, 3); // x-fnode[3] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  // Before
  old_nx = 0.0;
  old_ny = -1.0;
  // After rotation
  theta = -M_PI/4;
  nx = cos(theta)*old_nx - sin(theta)*old_ny;
  ny = sin(theta)*old_nx + cos(theta)*old_ny;

  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));
  
  norm = e->get_normal_vector(1, 0, 3);      // x-fnode[0] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));

  std::cout << "Element 1)" << "\n";

  // Element 1:
  e = mesh->elems[1];
  // Neighboors
  // edge 0
  ed = e->edges[0];
  // REQUIRE(ed.nodes[0] == 8);
  // REQUIRE(ed.nodes[1] == 5);
  // REQUIRE(ed.ghost == -1);
  // REQUIRE(ed.left == 1);
  // REQUIRE(ed.right == 2);
  // REQUIRE(ed.lr_edge == 3);
  // REQUIRE(ed.boundary == 0);
  // fn = ed.fnodes[0];
  // REQUIRE(fn.local == 0);
  // REQUIRE(fn.right == 0);
  // fn = ed.fnodes[1];
  // REQUIRE(fn.local == 3);
  // REQUIRE(fn.right == 3);
  // Normal
  norm = e->get_normal_vector(2, 0, 0); // y-fnode[0] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  // Before
  old_nx = 0.0;
  old_ny = 1.0;
  // After rotation
  theta = -M_PI/4;
  nx = cos(theta)*old_nx - sin(theta)*old_ny;
  ny = sin(theta)*old_nx + cos(theta)*old_ny;

  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));
  
  norm = e->get_normal_vector(2, 3, 0);      // y-fnode[3] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));

  // edge 1
  ed = e->edges[1];
  // REQUIRE(ed.nodes[0] == 5);
  // REQUIRE(ed.nodes[1] == 2);
  // REQUIRE(ed.ghost == 2);
  // REQUIRE(ed.left == 1);
  // REQUIRE(ed.right == -1);
  // REQUIRE(ed.lr_edge == -1);
  // REQUIRE(ed.boundary == 1);
  // fn = ed.fnodes[0];
  // REQUIRE(fn.local == 2);
  // REQUIRE(fn.right == 2);
  // fn = ed.fnodes[1];
  // REQUIRE(fn.local == 5);
  // REQUIRE(fn.right == 5);
  // Normal
  norm = e->get_normal_vector(1, 2, 1); // x-fnode[2] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  // Before
  old_nx = -1.0;
  old_ny = 0.0;
  // After rotation
  theta = -M_PI/4;
  nx = cos(theta)*old_nx - sin(theta)*old_ny;
  ny = sin(theta)*old_nx + cos(theta)*old_ny;

  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));
  
  norm = e->get_normal_vector(1, 5, 1);      // x-fnode[5] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));

  // edge 2
  ed = e->edges[2];
  // REQUIRE(ed.nodes[0] == 2);
  // REQUIRE(ed.nodes[1] == 6);
  // REQUIRE(ed.ghost == 3);
  // REQUIRE(ed.left == 1);
  // REQUIRE(ed.right == -1);
  // REQUIRE(ed.lr_edge == -1);
  // REQUIRE(ed.boundary == 1);
  // fn = ed.fnodes[0];
  // REQUIRE(fn.local == 5);
  // REQUIRE(fn.right == 5);
  // fn = ed.fnodes[1];
  // REQUIRE(fn.local == 2);
  // REQUIRE(fn.right == 2);
  // Normal
  norm = e->get_normal_vector(2, 5, 2); // y-fnode[5] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  // Before
  old_nx = 0.0;
  old_ny = -1.0;
  // After rotation
  theta = -M_PI/4;
  nx = cos(theta)*old_nx - sin(theta)*old_ny;
  ny = sin(theta)*old_nx + cos(theta)*old_ny;

  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));
  
  norm = e->get_normal_vector(2, 2, 2);      // y-fnode[2] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));

  // edge 3
  ed = e->edges[3];
  // REQUIRE(ed.nodes[0] == 6);
  // REQUIRE(ed.nodes[1] == 8);
  // REQUIRE(ed.ghost == -1);
  // REQUIRE(ed.left == 1);
  // REQUIRE(ed.right == 0);
  // REQUIRE(ed.lr_edge == 0);
  // REQUIRE(ed.boundary == 0);
  // fn = ed.fnodes[0];
  // REQUIRE(fn.local == 3);
  // REQUIRE(fn.right == 3);
  // fn = ed.fnodes[1];
  // REQUIRE(fn.local == 0);
  // REQUIRE(fn.right == 0);
  // Normal
  norm = e->get_normal_vector(1, 3, 3); // x-fnode[3] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  // Before
  old_nx = 1.0;
  old_ny = 0.0;
  // After rotation
  theta = -M_PI/4;
  nx = cos(theta)*old_nx - sin(theta)*old_ny;
  ny = sin(theta)*old_nx + cos(theta)*old_ny;
  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));
  
  norm = e->get_normal_vector(1, 0, 3);      // x-fnode[0] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(nx).margin(1E-15));
  REQUIRE(norm[1] == Approx(ny).margin(1E-15));



  sd->to_vtk(mesh, 
            (cur_path.parent_path() / "results" / "pp_mesh_test_incl.vtk").string());
}



// std::vector<double> ringleb_field (const Node& n)
// {
//   std::vector<double> vec(4, 0.0);

//   // Coordinates
//   double x = n.coords[0];
//   double y = n.coords[1];
//   double q = sqrt(x*x + y*y); // Umag
//   double gamma = 1.4;

//   double a = sqrt(1.0 - (gamma-1.0)*q*q/2.0);

//   double rho = pow(a, (2.0/(gamma-1.0)));
//   double u = x;
//   double v = y;
//   double p = (1.0/gamma)*pow(a, (2.0*gamma/(gamma-1.0)));
//   double E = p/(gamma-1.0) + rho*(x*x + y*y)/2.0;
  
//   vec[0] = rho; 
//   vec[1] = rho*u; 
//   vec[2] = rho*v; 
//   vec[3] = E; 
  
//   return vec;

// }



TEST_CASE("7: Test SD<Euler> - setup", "[sd_ringleb]")
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "ringleb_v0.msh").string());

  int order = 10;
  auto sd = std::make_shared<SD<Euler>>(order, 2);

  
  //sd->setup(mesh, ringleb_field);
  //sd->setup(mesh, FIELDS::LINEAR_FIELD_MAPPING);
  sd->setup(mesh, FIELDS::RINGLEB_FIELD_MAPPING);

  sd->to_vtk(mesh, 
            (cur_path.parent_path() / "results" / "pp_mesh_test_ringleb.vtk").string());
}

TEST_CASE("8: Test SD<Euler> - setup", "[sduniform_2order]")
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "mesh_test_2order.msh").string());

  int order = 2;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  
  sd->setup(mesh, FIELDS::DEFAULT_FIELD_MAPPING);

  // Element 0:
  auto e = mesh->elems[0];
  // snodes
  // sp 0
  auto n = e->transform(sd->snodes[0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.3943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 1
  n = e->transform(sd->snodes[1], mesh->nodes);
  

  REQUIRE(n.coords[0] == Approx(+0.1056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 2
  n = e->transform(sd->snodes[2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.3943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 3
  n = e->transform(sd->snodes[3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.1056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  
  // x-fnodes
  // fp 0
  n = e->transform(sd->fnodes[0][0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.3943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));

  // fp 1
  n = e->transform(sd->fnodes[0][1], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.3943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 2
  n = e->transform(sd->fnodes[0][2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.3943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 3
  n = e->transform(sd->fnodes[0][3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.1056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 4
  n = e->transform(sd->fnodes[0][4], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.1056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 5
  n = e->transform(sd->fnodes[0][5], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.1056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));

  // y-fnodes
  // fp 0
  n = e->transform(sd->fnodes[1][0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 1
  n = e->transform(sd->fnodes[1][1], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.250).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 2
  n = e->transform(sd->fnodes[1][2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.000).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 3
  n = e->transform(sd->fnodes[1][3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 4
  n = e->transform(sd->fnodes[1][4], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.250).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 5
  n = e->transform(sd->fnodes[1][5], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.000).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));

  // Neighboors
  // edge 0
  auto ed = e->edges[0];
  REQUIRE(ed.nodes[0] == 16);
  REQUIRE(ed.nodes[1] == 10);
  REQUIRE(ed.ghost == -1);
  REQUIRE(ed.left == 0);
  REQUIRE(ed.right == 1);
  REQUIRE(ed.lr_edge == 3);
  REQUIRE(ed.boundary == 0);
  auto& fn = ed.fnodes[0];
  REQUIRE(fn.local == 0);
  REQUIRE(fn.right == 0);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 3);
  REQUIRE(fn.right == 3);
  // Normal
  auto norm = e->get_normal_vector(2, 0, 0); // y-fnode[0] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));
  
  norm = e->get_normal_vector(2, 3, 0);      // y-fnode[3] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));

  // edge 1
  ed = e->edges[1];
  REQUIRE(ed.nodes[0] == 10);
  REQUIRE(ed.nodes[1] == 3);
  REQUIRE(ed.ghost == 0);
  REQUIRE(ed.left == 0);
  REQUIRE(ed.right == -1);
  REQUIRE(ed.lr_edge == -1);
  REQUIRE(ed.boundary == 1);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 2);
  REQUIRE(fn.right == 2);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 5);
  REQUIRE(fn.right == 5);
  // Normal
  norm = e->get_normal_vector(1, 2, 1); // x-fnode[2] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+1.0).margin(1E-15));
  
  norm = e->get_normal_vector(1, 5, 1);      // x-fnode[5] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+1.0).margin(1E-15));
  
  // edge 2
  ed = e->edges[2];
  REQUIRE(ed.nodes[0] == 3);
  REQUIRE(ed.nodes[1] == 13);
  REQUIRE(ed.ghost == 1);
  REQUIRE(ed.left == 0);
  REQUIRE(ed.right == -1);
  REQUIRE(ed.lr_edge == -1);
  REQUIRE(ed.boundary == 1);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 5);
  REQUIRE(fn.right == 5);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 2);
  REQUIRE(fn.right == 2);
  // Normal
  norm = e->get_normal_vector(2, 5, 2); // y-fnode[5] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(-1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));
  
  norm = e->get_normal_vector(2, 2, 2);      // y-fnode[2] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(-1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));

  // edge 3
  ed = e->edges[3];
  REQUIRE(ed.nodes[0] == 13);
  REQUIRE(ed.nodes[1] == 16);
  REQUIRE(ed.ghost == -1);
  REQUIRE(ed.left == 0);
  REQUIRE(ed.right == 3);
  REQUIRE(ed.lr_edge == 3);
  REQUIRE(ed.boundary == 0);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 3);
  REQUIRE(fn.right == 0);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 0);
  REQUIRE(fn.right == 3);
  // Normal
  norm = e->get_normal_vector(1, 3, 3); // x-fnode[3] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(-1.0).margin(1E-15));
  
  norm = e->get_normal_vector(1, 0, 3);      // x-fnode[0] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(-1.0).margin(1E-15));

  std::cout << "Element 1)" << "\n";

  // Element 1:
  e = mesh->elems[1];
  // snodes
  // sp 0
  n = e->transform(sd->snodes[0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 1
  n = e->transform(sd->snodes[1], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 2
  n = e->transform(sd->snodes[2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // sp 3
  n = e->transform(sd->snodes[3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  
  // x-fnodes
  // fp 0
  n = e->transform(sd->fnodes[0][0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 1
  n = e->transform(sd->fnodes[0][1], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 2
  n = e->transform(sd->fnodes[0][2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 3
  n = e->transform(sd->fnodes[0][3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 4
  n = e->transform(sd->fnodes[0][4], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 5
  n = e->transform(sd->fnodes[0][5], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));

  // y-fnodes
  // fp 0
  n = e->transform(sd->fnodes[1][0], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 1
  n = e->transform(sd->fnodes[1][1], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 2
  n = e->transform(sd->fnodes[1][2], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.6056624327).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 3
  n = e->transform(sd->fnodes[1][3], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.500).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 4
  n = e->transform(sd->fnodes[1][4], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+0.750).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));
  // fp 5
  n = e->transform(sd->fnodes[1][5], mesh->nodes);
  REQUIRE(n.coords[0] == Approx(+0.8943375673).margin(1E-15));
  REQUIRE(n.coords[1] == Approx(+1.000).margin(1E-15));
  REQUIRE(n.coords[2] == Approx(+0.000).margin(1E-15));

  // Neighboors
  // edge 0
  ed = e->edges[0];
  REQUIRE(ed.nodes[0] == 16);
  REQUIRE(ed.nodes[1] == 7);
  REQUIRE(ed.ghost == -1);
  REQUIRE(ed.left == 1);
  REQUIRE(ed.right == 2);
  REQUIRE(ed.lr_edge == 3);
  REQUIRE(ed.boundary == 0);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 0);
  REQUIRE(fn.right == 0);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 3);
  REQUIRE(fn.right == 3);
  // Normal
  norm = e->get_normal_vector(2, 0, 0); // y-fnode[0] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(-1.0).margin(1E-15));
  
  norm = e->get_normal_vector(2, 3, 0);      // y-fnode[3] at edge 0
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(-1.0).margin(1E-15));

  // edge 1
  ed = e->edges[1];
  REQUIRE(ed.nodes[0] == 7);
  REQUIRE(ed.nodes[1] == 2);
  REQUIRE(ed.ghost == 2);
  REQUIRE(ed.left == 1);
  REQUIRE(ed.right == -1);
  REQUIRE(ed.lr_edge == -1);
  REQUIRE(ed.boundary == 1);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 2);
  REQUIRE(fn.right == 2);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 5);
  REQUIRE(fn.right == 5);
  // Normal
  norm = e->get_normal_vector(1, 2, 1); // x-fnode[2] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));
  
  norm = e->get_normal_vector(1, 5, 1);      // x-fnode[5] at edge 1
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(+1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));

  // edge 2
  ed = e->edges[2];
  REQUIRE(ed.nodes[0] == 2);
  REQUIRE(ed.nodes[1] == 10);
  REQUIRE(ed.ghost == 3);
  REQUIRE(ed.left == 1);
  REQUIRE(ed.right == -1);
  REQUIRE(ed.lr_edge == -1);
  REQUIRE(ed.boundary == 1);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 5);
  REQUIRE(fn.right == 5);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 2);
  REQUIRE(fn.right == 2);
  // Normal
  norm = e->get_normal_vector(2, 5, 2); // y-fnode[5] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+1.0).margin(1E-15));
  
  norm = e->get_normal_vector(2, 2, 2);      // y-fnode[2] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(0.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+1.0).margin(1E-15));

  // edge 3
  ed = e->edges[3];
  REQUIRE(ed.nodes[0] == 10);
  REQUIRE(ed.nodes[1] == 16);
  REQUIRE(ed.ghost == -1);
  REQUIRE(ed.left == 1);
  REQUIRE(ed.right == 0);
  REQUIRE(ed.lr_edge == 0);
  REQUIRE(ed.boundary == 0);
  fn = ed.fnodes[0];
  REQUIRE(fn.local == 3);
  REQUIRE(fn.right == 3);
  fn = ed.fnodes[1];
  REQUIRE(fn.local == 0);
  REQUIRE(fn.right == 0);
  // Normal
  norm = e->get_normal_vector(1, 3, 3); // x-fnode[3] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(-1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));
  
  norm = e->get_normal_vector(1, 0, 3);      // x-fnode[0] at edge 2
  REQUIRE(sqrt(norm[0]*norm[0]+norm[1]*norm[1]) == Approx(+1.0).margin(1E-15)); // unit norm
  REQUIRE(norm[0] == Approx(-1.0).margin(1E-15));
  REQUIRE(norm[1] == Approx(+0.0).margin(1E-15));



  sd->to_vtk(mesh, 
            (cur_path.parent_path() / "results" / "pp_mesh_test_2_order.vtk").string());
}
