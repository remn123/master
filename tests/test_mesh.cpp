
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>

#include <catch/catch.hpp>
#include <Dummies.h>
#include <Element.h>
#include <Mesh.h>
#include <Node.h>
#include <SD.h>

using namespace Catch::literals;
namespace fs = std::filesystem;

TEST_CASE("1: Test 1 High Order Mesh", "[mesh]")
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "ringleb_v0.msh").string());

}