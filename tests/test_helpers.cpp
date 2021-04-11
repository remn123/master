
// 020-TestCase-2.cpp

// main() provided by Catch in file 020-TestCase-1.cpp.
#include <iostream>
#include <sstream>

#include <catch/catch.hpp>
#include <Helpers.h>
#include <Poly.h>

using namespace Catch::literals;

// POLY
TEST_CASE("1: Test Poly - factorial Zero case", "[Poly]")
{
  GL pol{};
  REQUIRE( pol.factorial(0) == 1);
}

TEST_CASE("2: Test Poly - factorial of x", "[Poly]")
{
  GL pol{};
  REQUIRE( pol.factorial(5) == 120);
}


TEST_CASE("3: Test Poly - factorial memoization", "[Poly]")
{
  
  GL pol{};
    
  long x1 = pol.factorial(5);

  long x2 = pol.factorial(5);
  REQUIRE(x1 == x2);

}

TEST_CASE("4: Test Poly - Pn zero case", "[Poly]")
{
  GL pol{};

  double x = pol.Pn(0.0, 0.0, 0, 0.0); // Guass-Legendre
  REQUIRE(x == 1.0_a);
}

TEST_CASE("5: Test Poly - dPn zero case", "[Poly]")
{
  GL pol{};

  double x = pol.dPn(0.0, 0.0, 0, 0.0); // Guass-Legendre
  REQUIRE(x == 0.0_a);
}

TEST_CASE("6: Test Poly - Pn root case", "[Poly]")
{
  GL pol{};

  double root = pol.Pn(0.0, 0.0, 2, -0.5773502691896257); // Guass-Legendre
  REQUIRE(root == Approx(0.0).margin(1E-5));
}
// ---------------------------------------------------------------- //

// Chebyshev
TEST_CASE("1: Test setup","[chebyshev]")
{
  Chebyshev cheb {};
  cheb.setup(2);

  REQUIRE(cheb.get_node(0) == -0.707107_a);
  REQUIRE(cheb.get_node(1) == 0.707107_a);
}

// Gauss-Legendre
TEST_CASE("4: Test Legendre - get_node method", "[gauss-legendre]")
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(2);
  REQUIRE( Helpers<GL>::get_node(0) == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE( Helpers<GL>::get_node(1) == Approx(0.577350269189626).margin(1E-15));
  Helpers<GL>::delete_nodes();
}

TEST_CASE("5: Test Legendre - nodes length", "[gauss-legendre]")
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(2);
  REQUIRE( Helpers<GL>::get_nodes().size() == 2);
  Helpers<GL>::delete_nodes();
}

TEST_CASE("6: Test Legendre - get_nodes method", "[gauss-legendre]")
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(2);
  std::vector<double> nodes = Helpers<GL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE( nodes[1] == Approx(0.577350269189626).margin(1E-15));
  Helpers<GL>::delete_nodes();
}

TEST_CASE("7: Test Legendre - set_nodes", "[gauss-legendre]")
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(2);
  std::vector<double> nodes = Helpers<GL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-0.577350269189626).margin(1E-15));
  REQUIRE( nodes[1] == Approx(0.577350269189626).margin(1E-15));
  Helpers<GL>::delete_nodes();

  nodes.clear();

  Helpers<GL>::set_nodes(3);
  nodes = Helpers<GL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-0.7745966692414833770359).margin(1E-15));
  REQUIRE( nodes[1] == Approx(0.0).margin(1E-15));
  REQUIRE( nodes[2] == Approx(0.7745966692414833770359).margin(1E-15));
  Helpers<GL>::delete_nodes(); 

  nodes.clear();

  Helpers<GL>::set_nodes(4);
  nodes = Helpers<GL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-0.8611363115940525752239).margin(1E-15));
  REQUIRE( nodes[1] == Approx(-0.3399810435848562648026).margin(1E-15));
  REQUIRE( nodes[2] == Approx(0.3399810435848562648026).margin(1E-15));
  REQUIRE( nodes[3] == Approx(0.8611363115940525752239).margin(1E-15));
  Helpers<GL>::delete_nodes(); 

  nodes.clear();

  Helpers<GL>::set_nodes(10);
  nodes = Helpers<GL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-0.9739065285171717200779640120844520534).margin(1E-15));
  REQUIRE( nodes[1] == Approx(-0.8650633666889845107320966884234930485).margin(1E-15));
  REQUIRE( nodes[2] == Approx(-0.6794095682990244062343273651148735757).margin(1E-15));
  REQUIRE( nodes[3] == Approx(-0.4333953941292471907992659431657841622).margin(1E-15));
  REQUIRE( nodes[4] == Approx(-0.1488743389816312108848260011297199846).margin(1E-15));
  REQUIRE( nodes[5] == Approx(0.1488743389816312108848260011297199846).margin(1E-15));
  REQUIRE( nodes[6] == Approx(0.4333953941292471907992659431657841622).margin(1E-15));
  REQUIRE( nodes[7] == Approx(0.6794095682990244062343273651148735757).margin(1E-15));
  REQUIRE( nodes[8] == Approx(0.8650633666889845107320966884234930485).margin(1E-15));
  REQUIRE( nodes[9] == Approx(0.9739065285171717200779640120844520534).margin(1E-15));
  Helpers<GL>::delete_nodes(); 
}
// ----------------------------------------------------------------- //


// Gauss-Legendre-Lobatto
TEST_CASE("1: Test Gauss-Lobatto - get_node method", "[gauss-lobatto]")
{
  Helpers<GLL>::init();
  Helpers<GLL>::set_nodes(2);
  REQUIRE( Helpers<GLL>::get_node(0) == -1.0);
  REQUIRE( Helpers<GLL>::get_node(1) == 1.0);
  Helpers<GLL>::delete_nodes();
}

TEST_CASE("2: Test Gauss-Lobatto - nodes length", "[gauss-lobatto]")
{
  Helpers<GLL>::init();
  Helpers<GLL>::set_nodes(2);
  REQUIRE( Helpers<GLL>::get_nodes().size() == 2);
  Helpers<GLL>::delete_nodes();
}

TEST_CASE("3: Test Gauss-Lobatto - get_nodes method", "[gauss-lobatto]")
{
  Helpers<GLL>::init();
  Helpers<GLL>::set_nodes(3);
  std::vector<double> nodes = Helpers<GLL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-1.0).margin(1E-15));
  REQUIRE( nodes[1] == Approx(0.0).margin(1E-15));
  REQUIRE( nodes[2] == Approx(1.0).margin(1E-15));
  Helpers<GLL>::delete_nodes();
}

TEST_CASE("4: Test Lobatto - high orders", "[gauss-lobatto]")
{
  Helpers<GLL>::init();
  Helpers<GLL>::set_nodes(2);
  std::vector<double> nodes = Helpers<GLL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-1.0).margin(1E-15));
  REQUIRE( nodes[1] == Approx(1.0).margin(1E-15));
  Helpers<GLL>::delete_nodes();

  nodes.clear();

  Helpers<GLL>::set_nodes(3);
  nodes = Helpers<GLL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-1.0).margin(1E-15));
  REQUIRE( nodes[1] == Approx(0.0).margin(1E-15));
  REQUIRE( nodes[2] == Approx(1.0).margin(1E-15));
  Helpers<GLL>::delete_nodes(); 

  nodes.clear();

  Helpers<GLL>::set_nodes(4);
  nodes = Helpers<GLL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-1.0).margin(1E-15));
  REQUIRE( nodes[1] == Approx(-0.5773502691896257645092).margin(1E-15));
  REQUIRE( nodes[2] == Approx(0.5773502691896257645092).margin(1E-15));
  REQUIRE( nodes[3] == Approx(1.0).margin(1E-15));
  Helpers<GLL>::delete_nodes(); 

  nodes.clear();

  Helpers<GLL>::set_nodes(5);
  nodes = Helpers<GLL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-1.0).margin(1E-15));
  REQUIRE( nodes[1] == Approx(-0.7745966692414833770359).margin(1E-15));
  REQUIRE( nodes[2] == Approx(0.0).margin(1E-15));
  REQUIRE( nodes[3] == Approx(0.7745966692414833770359).margin(1E-15));
  REQUIRE( nodes[4] == Approx(1.0).margin(1E-15));
  Helpers<GLL>::delete_nodes();

  nodes.clear();

  Helpers<GLL>::set_nodes(10);
  nodes = Helpers<GLL>::get_nodes();
  REQUIRE( nodes[0] == Approx(-1.0).margin(1E-15));
  REQUIRE( nodes[1] == Approx(-0.9602898564975362316835608685694729904).margin(1E-15));
  REQUIRE( nodes[2] == Approx(-0.7966664774136267395915539364758304368).margin(1E-15));
  REQUIRE( nodes[3] == Approx(-0.5255324099163289858177390491892463490).margin(1E-15));
  REQUIRE( nodes[4] == Approx(-0.1834346424956498049394761423601839806).margin(1E-15));
  REQUIRE( nodes[5] == Approx(0.1834346424956498049394761423601839806).margin(1E-15));
  REQUIRE( nodes[6] == Approx(0.5255324099163289858177390491892463490).margin(1E-15));
  REQUIRE( nodes[7] == Approx(0.7966664774136267395915539364758304368).margin(1E-15));
  REQUIRE( nodes[8] == Approx(0.9602898564975362316835608685694729904).margin(1E-15));
  REQUIRE( nodes[9] == Approx(1.0).margin(1E-15));
  Helpers<GLL>::delete_nodes(); 
}
// ----------------------------------------------------------------- //

TEST_CASE("1: Test GL and GLL - get_nodes method", "[GL-GLL]")
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(2);
  std::vector<double> nodes_gl = Helpers<GL>::get_nodes();
  REQUIRE( nodes_gl[0] == Approx(-0.5773502691896257645092).margin(1E-15));
  REQUIRE( nodes_gl[1] == Approx(0.5773502691896257645092).margin(1E-15));
  Helpers<GL>::delete_nodes();
  
  Helpers<GLL>::init();
  Helpers<GLL>::set_nodes(3);
  std::vector<double> nodes_gll = Helpers<GLL>::get_nodes();
  REQUIRE( nodes_gll[0] == Approx(-1.0).margin(1E-15));
  REQUIRE( nodes_gll[1] == Approx(0.0).margin(1E-15));
  REQUIRE( nodes_gll[2] == Approx(1.0).margin(1E-15));
  Helpers<GLL>::delete_nodes();

  nodes_gl.clear();
  nodes_gll.clear();
}

// Lagrange
 TEST_CASE("1: Test Lagrange - get_node method", "[lagrange]")
 {
   // First I will create a Legendre polynomial 
   Helpers<GL>::init();
   Helpers<GL>::set_nodes(2);
   std::vector<double> nodes_gl = Helpers<GL>::get_nodes();

   // I create a Lagrange polynomial using the gauss-legendre roots as nodes
   Helpers<Lagrange>::init();
   Helpers<Lagrange>::set_nodes(nodes_gl);
   REQUIRE( Helpers<Lagrange>::Pn(0, nodes_gl[0]) == 1.0);
   REQUIRE( Helpers<Lagrange>::Pn(0, nodes_gl[1]) == 0.0);
   REQUIRE( Helpers<Lagrange>::Pn(1, nodes_gl[0]) == 0.0);
   REQUIRE( Helpers<Lagrange>::Pn(1, nodes_gl[1]) == 1.0);

   Helpers<Lagrange>::delete_nodes();
   Helpers<GL>::delete_nodes();
   nodes_gl.clear();
 }
// ----------------------------------------------------------------- //

// Lagrange
// TEST_CASE("1: Test Lagrange - Pn method", "[lagrange]")
// {
//   // 
//   // 1) Initialize Lagrange polynomial with a set of nodes
//   // 
//   int Ns = 3;
//   Helpers<GL>::init();
//   Helpers<GL>::set_nodes(Ns);
//   std::vector<double> nodes_gl = Helpers<GL>::get_nodes();
  
//   Helpers<Lagrange>::init();
//   Helpers<Lagrange>::set_nodes(nodes_gl); // initialize Lagrange polynomial with a set of nodes



//   for (auto i=0; i<= Ns; i++)
//   {
//     // All nodes are
//     L = Helpers<Lagrange>::Pn(i, j, nodes_gl[i]);
//     REQUIRE(  == Approx(n).margin(1E-15));
//   }
  
//   REQUIRE( Helpers<Lagrange>::get_node(1) == 1.0);
//   Helpers<Lagrange>::delete_nodes();
// }


//
//TEST_CASE("4: Test 2 Legendre Nodes", "[multi-file:2]")
//{
//    Helpers<GL>::init();
//    Helpers<GL>::set_nodes(2);
//    std::vector<double> nodes = Helpers<GL>::get_nodes();
//    REQUIRE( nodes[0] == -1.0/(sqrt(3.0)));
//    REQUIRE( nodes[1] == 1.0/(sqrt(3.0)));
//}
//
//TEST_CASE("5: Test 3 Legendre Nodes", "[multi-file:2]")
//{
//    Helpers<GL>::init();
//    Helpers<GL>::set_nodes(3);
//    std::vector<double> nodes = Helpers<GL>::get_nodes();
//    REQUIRE( nodes[0] == -sqrt(3.0/5.0));
//    REQUIRE( nodes[1] == 0.0);
//    REQUIRE( nodes[2] == sqrt(3.0/5.0));
//}
