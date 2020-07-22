#include <memory>

#include <catch/catch.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Time.h>

using namespace Catch::literals;
namespace ublas = boost::numeric::ublas;


TEST_CASE("1: Test UBlas - vector container", "[ublas]")
{
  ublas::vector<double> v (4);
  std::cout << "1: Test UBlas - vector container" << std::endl;

  for (size_t i = 0; i < v.size (); ++ i)
      v (i) = i;
  std::cout << v << std::endl;
  REQUIRE(v.size() == 4);
  REQUIRE(v[0] == 0.0);
  REQUIRE(v[1] == 1.0);
  REQUIRE(v[2] == 2.0);
  REQUIRE(v[3] == 3.0);
}


TEST_CASE("2: Test Ublas - unit_vector", "[ublas]")
{
  std::cout << "2: Test Ublas - unit_vector" << std::endl;
  for (int i = 0; i < 3; ++ i) 
  {
    ublas::unit_vector<double> v (3, i);
    std::cout << v << std::endl;
  }
}

TEST_CASE("3: Test Ublas - zero_vector", "[ublas]")
{
  std::cout << "3: Test Ublas - zero_vector" << std::endl;
  ublas::zero_vector<double> v (3);
  std::cout << v << std::endl;
}

TEST_CASE("4: Test Ublas - scalar_vector", "[ublas]")
{
  std::cout << "4: Test Ublas - scalar_vector" << std::endl;
  ublas::scalar_vector<double> v (3);
  std::cout << v << std::endl;
}


// TEST_CASE("2: Test UBlas - vector of vectors", "[ublas]")
// {
  

//   std::vector<ublas::vector<double>> Q;
//   Q(1.0, 2.0 , 7.0, -4.0);
//   (0.0, 2.0 , 8.0,  4.0);
//   (1.5, 1.0 , 9.0, -3.0);

//   REQUIRE(Q.size() == 3);
//   REQUIRE(Q[0].size() == 4);
//   REQUIRE(Q[0][0] == 1.0);
//   REQUIRE(Q[0][1] == 2.0);
//   REQUIRE(Q[0][2] == 7.0);
//   REQUIRE(Q[0][3] == -4.0);

//   REQUIRE(Q[1][0] == 0.0);
//   REQUIRE(Q[1][1] == 2.0);
//   REQUIRE(Q[1][2] == 8.0);
//   REQUIRE(Q[1][3] == 4.0);

//   REQUIRE(Q[2][0] == 1.5);
//   REQUIRE(Q[2][1] == 1.0);
//   REQUIRE(Q[2][2] == 9.0);
//   REQUIRE(Q[2][3] == -3.0);
// }