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

TEST_CASE("2: Test dict - split_params", "[dict]")
{
  dict params = {
    {"d_path",      "/usr/dev/mestrado/results/log/cylinder.log"},
    {"t_max_iter",  "100000"},
    {"t_CFL",       "0.5"},
    {"s_order",     "2"},
    {"s_dimension", "2"},
    {"d_verbose",   "1"}        // 0-None, 1-Print, 2-Log
    };

  dict space_params   = split_params(params, "s_");
  dict time_params    = split_params(params, "t_");
  dict default_params = split_params(params, "d_");

  REQUIRE( default_params["path"] == "/usr/dev/mestrado/results/log/cylinder.log");
  REQUIRE( default_params["verbose"] == "1");

  REQUIRE( space_params["order"] == "2");
  REQUIRE( space_params["dimension"] == "2");

  REQUIRE( time_params["max_iter"] == "100000");
  REQUIRE( time_params["CFL"] == "0.5");
}
