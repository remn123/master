#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <catch/catch.hpp>
#include <Dummies.h>
#include <Helpers.h>
#include <Poly.h>
#include <Mesh.h>
#include <SD.h>

using namespace Catch::literals;
namespace fs = std::filesystem;

// CONSTRUCTORS
// TEST_CASE("1: Test appending strings", "[posproc]")
// {

// std::vector<std::string> lines;

// lines.push_back("vou tomar cafeh ate o sono passar 1 !!");
// lines.push_back("vou tomar cafeh ate o sono passar 2 !!");
// lines.push_back("vou tomar cafeh ate o sono passar 3 !!");
// lines.push_back("vou tomar cafeh ate o sono passar 4 !!");
// lines.push_back("vou tomar cafeh ate o sono passar 5 !!");
// lines.push_back("");
// lines.push_back("sera que a linha em branco deu certo?");

// std::ofstream output;

// output.open("../results/test_posproc_lines.txt");
// for (auto &line : lines)
//   output << line << std::endl;
// output.close();
// }

TEST_CASE("1: Test string concatenation", "[posproc]")
{
  std::vector<std::string> lines;
  std::string result{"Node 10 has 4 vertices"};
  int n = 10;

  lines.push_back(std::string{"Node "} + std::to_string(n) + std::string{" has "} + std::to_string(4) + std::string{" vertices"});

  REQUIRE(lines[0] == result);
}

TEST_CASE("2: Test 1 n=2, vtk posproc", "[posproc]")
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "mesh_test_2.msh").string());

  int order = 3;
  auto sd = std::make_shared<SD<Euler>>(3, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  sd->setup(mesh);
  sd->to_vtk(mesh,
             (cur_path.parent_path() / "results" / "pp_mesh_test_2.vtk").string());
  // REQUIRE(lines[0] == result);
}

TEST_CASE("3: Test 2 n=2, vtk posproc", "[posproc]")
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "mesh_test.msh").string());

  int order = 3;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  sd->setup(mesh);
  sd->to_vtk(mesh,
             (cur_path.parent_path() / "results" / "pp_mesh_test.vtk").string());
  // REQUIRE(lines[0] == result);
}