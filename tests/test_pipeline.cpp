#include <iostream>
#include <experimental/filesystem>
#include <string>
#include <map>
#include <memory>
//#include <crtdbg.h>

#include <catch/catch.hpp>
#include <Dict.h>
#include <Mock.h>
#include <Pipeline.h>

using namespace Catch::literals;

typedef std::vector<std::shared_ptr<IPipe>> DAG;

// CONSTRUCTORS
TEST_CASE("1: Test Pipeline - ", "[pipeline]")
{

  std::shared_ptr<IPipe> preproc = std::make_shared<Pipe<Preproc<Mesh>>>();
  std::shared_ptr<IPipe> solver  = std::make_shared<Pipe<Solver<SD<Euler>, ExplicitEuler<>>>>();
  std::shared_ptr<IPipe> posproc = std::make_shared<Pipe<Posproc>>();

  Pipeline cylinder2D {DAG{preproc, solver, posproc}};

  std::vector<dict> dicts = {
    { // 1) Preproc dict
      {"path",      "/usr/dev/mestrado/resources/cylinder.msh"},
      {"type",      "1"},
      {"dimension", "2"},
      //{"static",    "0"},       // 0-Mesh, 1-Static_Mesh
      {"verbose",   "1"}        // 0-None, 1-Print, 2-Log
    },
    { // 2) Solver dict
      {"d_path",      "/usr/dev/mestrado/results/log/cylinder.log"},
      {"d_type",      "0"},
      {"t_max_iter", "100000"},
      {"t_CFL",       "0.5"},
      {"s_order",     "2"},
      {"s_dimension", "2"},
      {"d_verbose",   "1"}        // 0-None, 1-Print, 2-Log
    },
    { // 3) Posproc dict
      {"path",      "/usr/dev/mestrado/results/cylinder.vtk"},
      {"type",      "0"},
      {"verbose",   "1"}        // 0-None, 1-Print, 2-Log
    }
  };

  cylinder2D.setup(dicts);

  cylinder2D.start();
}
