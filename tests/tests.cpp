// 020-TestCase-2.cpp

// main() provided by Catch in file 020-TestCase-1.cpp.
#include <iostream>
#include <filesystem>
#include <string>

#include <catch/catch.hpp>
#include <Helpers.h>
#include <Mesh.h>

// #include "../Solver.h"

namespace fs = std::filesystem;

TEST_CASE("1: Check elements' id of a Mesh Object", "[multi-file:2]")
{
  fs::path cur_path = fs::current_path();

  Mesh mesh1{2};

  mesh1.read_gmsh((cur_path.parent_path() / "resources" / "mesh_test.msh").string());

  REQUIRE(mesh1.elems[0]->id == 0);
  REQUIRE(mesh1.elems[3]->id == 3);
}

TEST_CASE("2: Check nodes' id of an Element Object", "[mesh]")
{
  fs::path cur_path = fs::current_path();

  Mesh mesh1{2};

  mesh1.read_gmsh((cur_path.parent_path() / "resources" / "mesh_test.msh").string());

  long n0 = mesh1.elems[0]->nodes[0];
  long n1 = mesh1.elems[0]->nodes[1];
  long n2 = mesh1.elems[0]->nodes[2];
  long n3 = mesh1.elems[0]->nodes[3];

  SECTION("Check Element nodes component")
  {
    REQUIRE(n0 == 8);
    REQUIRE(n1 == 6);
    REQUIRE(n2 == 3);
    REQUIRE(n3 == 7);
  }

  SECTION("Check Mesh nodes id")
  {
    REQUIRE(mesh1.nodes[n0].id == 8);
    REQUIRE(mesh1.nodes[n1].id == 6);
    REQUIRE(mesh1.nodes[n2].id == 3);
    REQUIRE(mesh1.nodes[n3].id == 7);
  }
}

TEST_CASE("3: Check coordinates of a Node Object", "[mesh]")
{
  fs::path cur_path = fs::current_path();

  Mesh mesh1{2};

  mesh1.read_gmsh((cur_path.parent_path() / "resources" / "mesh_test.msh").string());

  long n = mesh1.elems[0]->nodes[0];

  REQUIRE(mesh1.nodes[n].coords[0] == 0.4983490244891706);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.4983490244891706);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  // for (auto& node: mesh1.elems[0]->nodes)
  //     std::cout << "Element[0]->Node[" << node << "]: ("
  //               << mesh1.nodes[node].coords[0] << ", "
  //               << mesh1.nodes[node].coords[1] << ", "
  //               << mesh1.nodes[node].coords[2] << ")\n";
}

TEST_CASE("4: Check edges of a Mesh Object", "[mesh2]")
{
  fs::path cur_path = fs::current_path();

  Mesh mesh1{2};

  mesh1.read_gmsh((cur_path.parent_path() / "resources" / "mesh_test.msh").string());

  // Check if element is a Quadrangle

  // Check the number of elements, nodes and edges
  REQUIRE(mesh1.elems.size() == 4);
  REQUIRE(mesh1.nodes.size() == 9);
  REQUIRE(mesh1.ghosts.size() == 8);
  REQUIRE(mesh1.Nel == 4);
  REQUIRE(mesh1.N == 9);
  REQUIRE(mesh1.Ngh == 8);
  REQUIRE(mesh1.Ned == 12);

  long n = 0;

  // Element 0: 8 6 3 7
  // Node 0
  n = mesh1.elems[0]->nodes[0];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.4983490244891706);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.4983490244891706);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 8);
  // Node 1
  n = mesh1.elems[0]->nodes[1];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.5000000000020595);
  REQUIRE(mesh1.nodes[n].coords[1] == 1.0);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 6);
  // Node 2
  n = mesh1.elems[0]->nodes[2];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.0);
  REQUIRE(mesh1.nodes[n].coords[1] == 1.0);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 3);
  // Node 3
  n = mesh1.elems[0]->nodes[3];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.0);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.5000000000020595);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 7);

  // Element 1: 8 5 2 6
  // Node 0
  n = mesh1.elems[1]->nodes[0];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.4983490244891706);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.4983490244891706);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 8);
  // Node 1
  n = mesh1.elems[1]->nodes[1];
  REQUIRE(mesh1.nodes[n].coords[0] == 1.0);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.4999999999986921);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 5);
  // Node 2
  n = mesh1.elems[1]->nodes[2];
  REQUIRE(mesh1.nodes[n].coords[0] == 1.0);
  REQUIRE(mesh1.nodes[n].coords[1] == 1.0);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 2);
  // Node 3
  n = mesh1.elems[1]->nodes[3];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.5000000000020595);
  REQUIRE(mesh1.nodes[n].coords[1] == 1.0);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 6);

  // Element 2: 8 4 1 5
  // Node 0
  n = mesh1.elems[2]->nodes[0];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.4983490244891706);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.4983490244891706);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 8);
  // Node 1
  n = mesh1.elems[2]->nodes[1];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.4999999999986921);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.0);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 4);
  // Node 2
  n = mesh1.elems[2]->nodes[2];
  REQUIRE(mesh1.nodes[n].coords[0] == 1.0);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.0);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 1);
  // Node 3
  n = mesh1.elems[2]->nodes[3];
  REQUIRE(mesh1.nodes[n].coords[0] == 1.0);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.4999999999986921);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 5);

  // Element 3: 7 0 4 8
  // Node 0
  n = mesh1.elems[3]->nodes[0];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.0);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.5000000000020595);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 7);
  // Node 1
  n = mesh1.elems[3]->nodes[1];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.0);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.0);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 0);
  // Node 2
  n = mesh1.elems[3]->nodes[2];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.4999999999986921);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.0);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 4);
  // Node 3
  n = mesh1.elems[3]->nodes[3];
  REQUIRE(mesh1.nodes[n].coords[0] == 0.4983490244891706);
  REQUIRE(mesh1.nodes[n].coords[1] == 0.4983490244891706);
  REQUIRE(mesh1.nodes[n].coords[2] == 0.0);
  REQUIRE(n == 8);

  // int local=0;
  // for (auto& e: mesh1.elems)
  // {
  //   std::cout << "Element[" << e->id << "]: \n";
  //   local = 0;
  //   for (auto& ed: e->edges)
  //   {
  //     std::cout << "    Edge(" << local << "/" << ed.id << "): "
  //               << ed.nodes[0] << " -> " << ed.nodes[1] << " | L = " << ed.left << " R = " << ed.right << "\n";
  //     local++;
  //   }
  // }
}

// 0.2758620689659038
// 0.7586206896550877
//TEST_CASE( "4: Test Pipeline", "[multi-file:2]" )
//{
//    fs::path cur_path = fs::current_path();

//    std::string infile1 = (cur_path.parent_path() / "data" / "mesh3.msh").string();
//    std::string infile2 = (cur_path.parent_path() / "data" / "mesh2.msh").string();
//    std::string outfile = (cur_path.parent_path() / "results" / "output.vtk").string();

/* Background Mesh */
//    dict pre_params1 = {
//        {"filename",  infile1},
//        {"dimension", 2},
//        {"static",    1},       // 0-Mesh, 1-Static_Mesh
//        {"verbose",   1}        // 0-None, 1-Print, 2-Log
//    };

/* Near-body Mesh */
//    dict pre_params2 = {
//        {"filename",  infile2},
//        {"dimension", 2},
//        {"static",    0},       // 0-Mesh, 1-Static_Mesh
//        {"verbose",   1}        // 0-None, 1-Print, 2-Log
//    };

//    dict sol_params = {
//        {"method",     2},   // 1-Jameson, 2-DG, 3-FR/CPR
//        {"equation",   1},   // 1-Euler, 2-Navier-Stokes (Haha)
//        {"turbulence", 0},   // 0-None, 1-RANS, 2-LES(HAHAHA)
//        {"limiter",    0},   // 0-None, 1-minmod, 2-superbee, 3-Albada
//        {"overset",    1},   // 0-None, 1-Face-Based, 2-Element-Based
//        {"acc_conv",   0},   // 0-None, 1-Multigrid, 2-GMRES (HAHAHAHAHAHA)
//        {"max_iter",   1E5}, // max number of iterations
//        {"verbose",    1}    // 0-None, 1-Print, 2-Log
//    };

//    dict pos_params = {
//        {"filename", outfile},
//        {"verbose",   1}
//    };

//    Pipeline test {std::vector<Pipe>{ Preproc{pre_params1},
//                                      Preproc{pre_params2},
//                                      Solver{sol_params},
//                                      Posproc{pos_params}}};

//    test.setup();
//    test.run();
//    test.save();

//}

// Compile: see 020-TestCase-1.cpp

// Expected compact output (all assertions):
//
// prompt> 020-TestCase --reporter compact --success
// 020-TestCase-2.cpp:13: failed: Factorial(0) == 1 for: 0 == 1
// 020-TestCase-2.cpp:17: passed: Factorial(1) == 1 for: 1 == 1
// 020-TestCase-2.cpp:18: passed: Factorial(2) == 2 for: 2 == 2
// 020-TestCase-2.cpp:19: passed: Factorial(3) == 6 for: 6 == 6
// 020-TestCase-2.cpp:20: passed: Factorial(10) == 3628800 for: 3628800 (0x375f00) == 3628800 (0x375f00)
// Failed 1 test case, failed 1 assertion.
