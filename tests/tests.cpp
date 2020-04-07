// 020-TestCase-2.cpp

// main() provided by Catch in file 020-TestCase-1.cpp.
#pragma once

#include <iostream>
#include <filesystem>
#include <string>

#include <catch/catch.hpp>
#include <Mesh.h>
#include <Helpers.h>
// #include "../Solver.h"

namespace fs = std::filesystem;


TEST_CASE( "1: Check elements' id of a Mesh Object", "[multi-file:2]" ) 
{
    fs::path cur_path = fs::current_path();
	
    Mesh mesh1 {2};

	mesh1.read_gmsh((cur_path.parent_path() / "resources" / "mesh1.msh").string());
	
    REQUIRE( mesh1.elems[0]->id == 0);
    REQUIRE( mesh1.elems[17]->id == 17);
}

TEST_CASE( "2: Check nodes' id of an Element Object", "[multi-file:2]" ) 
{
    fs::path cur_path = fs::current_path();
	
    Mesh mesh1 {2};

	mesh1.read_gmsh((cur_path.parent_path() / "resources" / "mesh1.msh").string());
	
    long n0 = mesh1.elems[0]->nodes[0];
    long n1 = mesh1.elems[0]->nodes[1];
    long n2 = mesh1.elems[0]->nodes[2];
    long n3 = mesh1.elems[0]->nodes[3];

    SECTION( "Check Element nodes component" ) {
        REQUIRE( n0 == 701);
        REQUIRE( n1 == 620);
        REQUIRE( n2 == 627);
        REQUIRE( n3 == 686);
    }

    SECTION( "Check Mesh nodes id" ) {
        REQUIRE( mesh1.nodes[n0].id == 701);
        REQUIRE( mesh1.nodes[n1].id == 620);
        REQUIRE( mesh1.nodes[n2].id == 627);
        REQUIRE( mesh1.nodes[n3].id == 686);
    }
}

TEST_CASE( "3: Check coordinates of a Node Object", "[multi-file:2]" ) 
{
    fs::path cur_path = fs::current_path();
	
    Mesh mesh1 {2};

	mesh1.read_gmsh((cur_path.parent_path() / "resources" / "mesh1.msh").string());
	
    long n = mesh1.elems[0]->nodes[0];

    REQUIRE( mesh1.nodes[n].coords[0] == 0.7241379310348262);
    REQUIRE( mesh1.nodes[n].coords[1] == 0.7241379310346457);
    REQUIRE( mesh1.nodes[n].coords[2] == 0.0);
    // for (auto& node: mesh1.elems[0]->nodes)
    //     std::cout << "Element[0]->Node[" << node << "]: (" 
    //               << mesh1.nodes[node].coords[0] << ", "
    //               << mesh1.nodes[node].coords[1] << ", "
    //               << mesh1.nodes[node].coords[2] << ")\n";
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
