#include <iostream>
#include <sstream>

#include <catch/catch.hpp>

using namespace Catch::literals;

double foo1(double x, double y)
{
  return x + 2.0 * y;
}

double foo2(double x, double y)
{
  return x + 5.0 * y;
}

double bar(double x, double y, double (*func)(double, double) = nullptr)
{
  if (func != nullptr)
    return func(x, y);
  return -1.0;
}

// Functional
TEST_CASE("1: Test functional", "[functional]")
{
  double z1 = bar(1.0, 2.0, foo1);
  double z2 = bar(1.0, 2.0, foo2);
  REQUIRE(z1 == 5.0);
  REQUIRE(z2 == 11.0);
}

// Functional
TEST_CASE("2: Test default functional", "[functional]")
{
  double z1 = bar(1.0, 2.0, foo1);
  double z2 = bar(1.0, 2.0);
  REQUIRE(z1 == 5.0);
  REQUIRE(z2 == -1.0);
}