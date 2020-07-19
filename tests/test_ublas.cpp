#include <catch/catch.hpp>
#include <boost/numeric/ublas/vector.hpp>
// #include <boost/numeric/ublas/io.hpp>

using namespace Catch::literals;
namespace ublas = boost::numeric::ublas;


TEST_CASE("1: Test UBlas - vector container", "[ublas]")
{
  ublas::vector<double> v (4);
  for (size_t i = 0; i < v.size (); ++ i)
      v (i) = i;
      
  REQUIRE(v.size() == 4);
  REQUIRE(v[0] == 0.0);
  REQUIRE(v[1] == 1.0);
  REQUIRE(v[2] == 2.0);
  REQUIRE(v[3] == 3.0);
}
