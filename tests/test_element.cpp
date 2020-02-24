#include <iostream>
#include <sstream>

#include <catch/catch.hpp>
#include <Mesh.h>

using namespace Catch::literals;

// CONSTRUCTORS
TEST_CASE("1: Test Quadrangle - constructor", "[elems]")
{
  Quadrangle q = {std::vector<std::string>{"1 2 3 4"}};
  REQUIRE(q.id == 1);
  REQUIRE(q.boundary == 0);
  REQUIRE(q.fringe == 0);
  REQUIRE(q.J == 0.0);
  REQUIRE(q.nodes.size() == 4);
  REQUIRE(q.nodes[0] == 0);
  REQUIRE(q.nodes[1] == 1);
  REQUIRE(q.nodes[2] == 2);
  REQUIRE(q.nodes[3] == 3);
}

TEST_CASE("1: Test Quadrangle - constructor", "[elems]")
{
  Quadrangle q = {std::vector<std::string>{"1 2 3 4"}};
  REQUIRE(q.id == 1);
  REQUIRE(q.boundary == 0);
  REQUIRE(q.fringe == 0);
  REQUIRE(q.J == 0.0);
  REQUIRE(q.nodes.size() == 4);
  REQUIRE(q.nodes[0] == 0);
  REQUIRE(q.nodes[1] == 1);
  REQUIRE(q.nodes[2] == 2);
  REQUIRE(q.nodes[3] == 3);
}
