
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

  REQUIRE(cheb.get_node(0) == 0.707107_a);
  REQUIRE(cheb.get_node(1) == -0.707107_a);
}


TEST_CASE("4: Test Legendre - get_node method", "[multi-file:2]")
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(2);
  REQUIRE( Helpers<GL>::get_node(0) == Approx(0.0).margin(1E-10));
}

//TEST_CASE("5: Test Poly - print_nodes method", "[multi-file:2]")
//{
//  GL::setup(5);
//  GL::print_nodes();
//}
//
TEST_CASE("6: Test Legendre - get_nodes method", "[multi-file:2]")
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(1);
  std::vector<double> nodes = Helpers<GL>::get_nodes();
  REQUIRE( nodes[0] == Approx(0.0).margin(1E-10));
}
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
