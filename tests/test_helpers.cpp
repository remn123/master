
// 020-TestCase-2.cpp

// main() provided by Catch in file 020-TestCase-1.cpp.
#include <iostream>
#include <sstream>

#include <catch/catch.hpp>
#include <Helpers.h>
#include <Poly.h>

TEST_CASE("1: Test Poly factorial Zero case", "[multi-file:3]")
{
  GL pol{};
  REQUIRE( pol.factorial(0) == 1);
}

TEST_CASE("2: Test Poly factorial of x", "[multi-file:3]")
{
  GL pol{};
  REQUIRE( pol.factorial(5) == 120);
}


TEST_CASE("3: Test Poly factorial memoization", "[multi-file:3]")
{
  
  GL pol{};
    
  long x1 = pol.factorial(5);

  long x2 = pol.factorial(5);
  REQUIRE(x1 == x2);

}


TEST_CASE("4: Test Legendre Nodes", "[multi-file:2]")
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(1);
  std::vector<double> nodes = Helpers<GL>::get_nodes();
  REQUIRE( nodes[0] == 0.0);
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
