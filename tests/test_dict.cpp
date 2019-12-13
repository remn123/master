#include <iostream>
#include <sstream>

#include <catch/catch.hpp>
#include <Dict.h>

using namespace Catch::literals;


// Dict
TEST_CASE("1: Test dict - dict namespace", "[dict]")
{
  dict params = {
    {"path",      "/usr/dev/mestrado"},
    {"type",      "1"},
    {"weight",    "0.5"},
    {"dimension", "2"},
    {"static",    "1"},       // 0-Mesh, 1-Static_Mesh
    {"verbose",   "1"}        // 0-None, 1-Print, 2-Log
    };

  REQUIRE( params["path"] == "/usr/dev/mestrado");
  REQUIRE( params["type"] == "1");
  REQUIRE( params["weight"] == "0.5");
  REQUIRE( params["dimension"] == "2");
  REQUIRE( params["static"] == "1");
  REQUIRE( params["verbose"] == "1");
}
