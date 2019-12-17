#include <iostream>
#include <sstream>

#include <catch/catch.hpp>
#include <DVector.h>

using namespace Catch::literals;

// CONSTRUCTORS
TEST_CASE("1: Test DVector - default constructor", "[dvector]")
{
  DVector d {};
  REQUIRE(d.size() == 4);
  REQUIRE(d[0] == 0.0);
  REQUIRE(d[1] == 0.0);
  REQUIRE(d[2] == 0.0);
  REQUIRE(d[3] == 0.0);
}

TEST_CASE("2: Test DVector - Parameterized constructor", "[dvector]")
{
  DVector d {std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0}};
  REQUIRE(d.size() == 5);
  REQUIRE(d[0] == 1.0);
  REQUIRE(d[1] == 2.0);
  REQUIRE(d[2] == 3.0);
  REQUIRE(d[3] == 4.0);
  REQUIRE(d[4] == 5.0);
}

TEST_CASE("3: Test DVector - Copy constructor", "[dvector]")
{
  DVector d0 {std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0}};
  DVector d{d0};

  REQUIRE(&d != &d0); // different objects
  REQUIRE(d.size() == 5);
  REQUIRE(d[0] == 1.0);
  REQUIRE(d[1] == 2.0);
  REQUIRE(d[2] == 3.0);
  REQUIRE(d[3] == 4.0);
  REQUIRE(d[4] == 5.0);
}

// ASSIGNEMENTS
TEST_CASE("4: Test DVector - Scalar assignement", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 1.0);
  REQUIRE(d[1] == 2.0);
  REQUIRE(d[2] == 3.0);
  
  d = 1.0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 1.0);
  REQUIRE(d[1] == 1.0);
  REQUIRE(d[2] == 1.0);
}

TEST_CASE("5: Test DVector - DVector assignement", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{0.0, 1.0, 5.0}};
  d = d0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 0.0);
  REQUIRE(d[1] == 1.0);
  REQUIRE(d[2] == 5.0);
}

TEST_CASE("6: Test DVector - std::vector assignement", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  std::vector<double> vec {0.0, 1.0, 5.0};
  d = vec;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 0.0);
  REQUIRE(d[1] == 1.0);
  REQUIRE(d[2] == 5.0);
}

// ADDITION
TEST_CASE("7: Test DVector - Addition operator += scalar (double)", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  d += 2.0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 3.0);
  REQUIRE(d[1] == 4.0);
  REQUIRE(d[2] == 5.0);
}

TEST_CASE("8: Test DVector - Addition operator += scalar (int)", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  d += 2;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 3.0);
  REQUIRE(d[1] == 4.0);
  REQUIRE(d[2] == 5.0);
}

TEST_CASE("9: Test DVector - Addition operator += DVector", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d += d0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 6.0);
  REQUIRE(d[1] == 4.0);
  REQUIRE(d[2] == 4.0);
}

TEST_CASE("10: Test DVector - Addition operator + scalar (double) from right", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d = d0 + 1.0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 6.0);
  REQUIRE(d[1] == 3.0);
  REQUIRE(d[2] == 2.0);
}

TEST_CASE("11: Test DVector - Addition operator + scalar (double) from left", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d = 1.0 + d0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 6.0);
  REQUIRE(d[1] == 3.0);
  REQUIRE(d[2] == 2.0);

  // d0 must remain the same
  REQUIRE(d0.size() == 3);
  REQUIRE(d0[0] == 5.0);
  REQUIRE(d0[1] == 2.0);
  REQUIRE(d0[2] == 1.0);
}

TEST_CASE("12: Test DVector - Addition operator + scalar (int)", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d = d0 + 1;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 6.0);
  REQUIRE(d[1] == 3.0);
  REQUIRE(d[2] == 2.0);

  // d0 must remain the same
  REQUIRE(d0.size() == 3);
  REQUIRE(d0[0] == 5.0);
  REQUIRE(d0[1] == 2.0);
  REQUIRE(d0[2] == 1.0);
}

TEST_CASE("13: Test DVector - Addition operator + DVector", "[dvector]")
{
  DVector d0{std::vector<double>{1.0, 3.0, 7.0}};
  DVector d1{std::vector<double>{5.0, 2.0, 1.0}};
  DVector d{d0 + d1};
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 6.0);
  REQUIRE(d[1] == 5.0);
  REQUIRE(d[2] == 8.0);
}

// SUBTRACTION
TEST_CASE("14: Test DVector - Subtraction operator -= scalar (double)", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  d -= 2.0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == -1.0);
  REQUIRE(d[1] == 0.0);
  REQUIRE(d[2] == 1.0);
}

TEST_CASE("15: Test DVector - Subtraction operator -= scalar (int)", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  d -= 2;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == -1.0);
  REQUIRE(d[1] == 0.0);
  REQUIRE(d[2] == 1.0);
}

TEST_CASE("16: Test DVector - Subtraction operator -= DVector", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d -= d0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == -4.0);
  REQUIRE(d[1] == 0.0);
  REQUIRE(d[2] == 2.0);
}

TEST_CASE("17: Test DVector - Subtraction operator - scalar (double) from right", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d = d0 - 1.0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 4.0);
  REQUIRE(d[1] == 1.0);
  REQUIRE(d[2] == 0.0);

  // d0 must remain the same
  REQUIRE(d0.size() == 3);
  REQUIRE(d0[0] == 5.0);
  REQUIRE(d0[1] == 2.0);
  REQUIRE(d0[2] == 1.0);
}

TEST_CASE("18: Test DVector - Subtraction operator - scalar (double) from left", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d = 1.0 - d0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == -4.0);
  REQUIRE(d[1] == -1.0);
  REQUIRE(d[2] == 0.0);

  REQUIRE(d0.size() == 3);
  REQUIRE(d0[0] == 5.0);
  REQUIRE(d0[1] == 2.0);
  REQUIRE(d0[2] == 1.0);
}

TEST_CASE("19: Test DVector - Subtraction operator - scalar (int)", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d = d0 - 1;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 4.0);
  REQUIRE(d[1] == 1.0);
  REQUIRE(d[2] == 0.0);
}

TEST_CASE("20: Test DVector - Subtraction operator - DVector", "[dvector]")
{
  DVector d0{std::vector<double>{1.0, 3.0, 7.0}};
  DVector d1{std::vector<double>{5.0, 2.0, 1.0}};
  DVector d{d0-d1};
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == -4.0);
  REQUIRE(d[1] == 1.0);
  REQUIRE(d[2] == 6.0);


  // d0 and d1 must remain the same
  REQUIRE(d0.size() == 3);
  REQUIRE(d0[0] == 1.0);
  REQUIRE(d0[1] == 3.0);
  REQUIRE(d0[2] == 7.0);

  REQUIRE(d1.size() == 3);
  REQUIRE(d1[0] == 5.0);
  REQUIRE(d1[1] == 2.0);
  REQUIRE(d1[2] == 1.0);
}

// MULTIPLICATION
TEST_CASE("21: Test DVector - Multiplcation operator *= scalar (double)", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  d *= 2.0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 2.0);
  REQUIRE(d[1] == 4.0);
  REQUIRE(d[2] == 6.0);
}

TEST_CASE("22: Test DVector - Multiplication operator *= scalar (int)", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  d *= 2;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 2.0);
  REQUIRE(d[1] == 4.0);
  REQUIRE(d[2] == 6.0);
}

TEST_CASE("23: Test DVector - Multiplication operator *= DVector", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d *= d0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 5.0);
  REQUIRE(d[1] == 4.0);
  REQUIRE(d[2] == 3.0);
}

TEST_CASE("24: Test DVector - Multiplication operator * scalar (double) from right", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d = d0 * 2.0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 10.0);
  REQUIRE(d[1] == 4.0);
  REQUIRE(d[2] == 2.0);
}

TEST_CASE("25: Test DVector - Multiplication operator * scalar (double) from left", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d = 2.0 * d0;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 10.0);
  REQUIRE(d[1] == 4.0);
  REQUIRE(d[2] == 2.0);
}

TEST_CASE("26: Test DVector - Multiplication operator * scalar (int)", "[dvector]")
{
  DVector d{std::vector<double>{1.0, 2.0, 3.0}};
  DVector d0{std::vector<double>{5.0, 2.0, 1.0}};
  d = d0 * 3;
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 15.0);
  REQUIRE(d[1] == 6.0);
  REQUIRE(d[2] == 3.0);
}

TEST_CASE("27: Test DVector - Elementwise multiplication operator * DVector", "[dvector]")
{
  DVector d0{std::vector<double>{1.0, 3.0, 7.0}};
  DVector d1{std::vector<double>{5.0, 2.0, 1.0}};
  DVector d{d0*d1};
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == 5.0);
  REQUIRE(d[1] == 6.0);
  REQUIRE(d[2] == 7.0);
}

TEST_CASE("28: Test DVector - negative DVector", "[dvector]")
{
  DVector d0{std::vector<double>{0.0, 0.0, 0.0}};
  DVector d1{std::vector<double>{5.0, 2.0, 1.0}};
  DVector d{-d1};
  REQUIRE(d.size() == 3);
  REQUIRE(d[0] == -5.0);
  REQUIRE(d[1] == -2.0);
  REQUIRE(d[2] == -1.0);

  REQUIRE(d1.size() == 3);
  REQUIRE(d1[0] == 5.0);
  REQUIRE(d1[1] == 2.0);
  REQUIRE(d1[2] == 1.0);
}