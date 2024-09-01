#include <memory>

#include <catch/catch.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <Time.h>
#include <Numericals.h>

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

TEST_CASE("5: Test Ublas - vector norm", "[ublas]")
{
  std::cout << "5: Test Ublas - vector norm" << std::endl;
  ublas::vector<double> v (3);  
  // [0, 1, 2] -> norm(v) = srqt(0*0 + 1*1 + 2*2) = sqrt(5) ~ 2.2360679775
  for (size_t i = 0; i < v.size(); ++ i)
    v(i) = i;
   
  std::cout << v << std::endl;
  auto norm = ublas::norm_2(v);
  std::cout << norm << std::endl;
  REQUIRE(norm == Approx(2.2360679775).margin(1E-15));
}

TEST_CASE("6: Test Ublas - vector inner_prod", "[ublas]")
{
  std::cout << "6: Test Ublas - vector inner_prod" << std::endl;
  ublas::vector<double> v1 (3);
  ublas::vector<double> v2 (3);
  
  for (size_t i = 0; i < v1.size(); ++ i)
  {
    v1(i) = i;
    v2(i) = 2.0*i+5.0;
  }
  // v1*v2 = 0.0*(2.0*0.0+5.0) + 1.0*(2.0*1.0+5.0) + 2.0*(2.0*2.0+5.0) = 0.0 + 7.0 + 18.0 = 25.0
  std::cout << v1 << std::endl;
  std::cout << v2 << std::endl;
  auto inner = ublas::inner_prod(v1, v2);
  std::cout << inner << std::endl;
  REQUIRE(inner == 25.0);
}

TEST_CASE("7: Test Ublas - matrix", "[ublas]")
{
  std::cout << "7: Test Ublas - matrix" << std::endl;
  ublas::matrix<double> m (3,3);
  std::cout << m << std::endl;
  for (size_t i = 0; i < m.size1(); ++i)
    for (size_t j = 0; j < m.size2(); ++j)
      m(i, j) = 3 * i + j;
  /*
    m = [[0.0, 1.0, 2.0],
         [3.0, 4.0, 5.0],
         [6.0, 7.0, 8.0]] 
  */
  std::cout << m << std::endl;
  REQUIRE(m(0,0) == 0.0);
  REQUIRE(m(0,1) == 1.0);
  REQUIRE(m(0,2) == 2.0);

  REQUIRE(m(1,0) == 3.0);
  REQUIRE(m(1,1) == 4.0);
  REQUIRE(m(1,2) == 5.0);

  REQUIRE(m(2,0) == 6.0);
  REQUIRE(m(2,1) == 7.0);
  REQUIRE(m(2,2) == 8.0);
}

TEST_CASE("8: Test Ublas - matrix prod vector", "[ublas]")
{
  std::cout << "8: Test Ublas - matrix prod vector" << std::endl;
  ublas::vector<double> v(3);
  ublas::matrix<double> m(3,3);

  for (size_t i = 0; i < m.size1(); ++i)
    for (size_t j = 0; j < m.size2(); ++j)
      m(i, j) = 3 * i + j;
  
  for (size_t i = 0; i < m.size1(); ++i)
    v(i) = i;

  /*
    v = [ 0.0,
          1.0, 
          2.0]
    m = [[0.0, 1.0, 2.0],
         [3.0, 4.0, 5.0],
         [6.0, 7.0, 8.0]] 
    m(3,3)*v(3, 1) = (3, 1)
    [ 0.0*0.0 + 1.0*1.0 + 2.0*2.0 = 5.0
      0.0*3.0 + 1.0*4.0 + 2.0*5.0 = 14.0
      0.0*6.0 + 1.0*7.0 + 2.0*8.0 = 23.0
    ]
  */
  v = prod (m, v);
  REQUIRE(v(0) == 5.0);
  REQUIRE(v(1) == 14.0);
  REQUIRE(v(2) == 23.0);
}

TEST_CASE("9: Test Ublas - matrix inverse", "[ublas]")
{
  std::cout << "9: Test Ublas - matrix inverse" << std::endl;
  ublas::matrix<double> m(2,2);
  ublas::matrix<double> inv_m(2,2);
  m(0,0) = 1.0;
  m(0,1) = 2.0;
  m(1,0) = 3.0;
  m(1,1) = 4.0;
  /*
    m      = [[1.0, 2.0],
              [3.0, 4.0]]
    
    inv(m) = [[ -2.0,      1.0  ],
              [3.0/2.0, -1.0/2.0]]
  
   */
  bool inverted = invert_matrix(m, inv_m);
  REQUIRE(inverted == true);
  REQUIRE(inv_m(0,0) == Approx(-2.0).margin(1E-15));
  REQUIRE(inv_m(0,1) == Approx(1.0).margin(1E-15));
  REQUIRE(inv_m(1,0) == Approx(3.0/2.0).margin(1E-15));
  REQUIRE(inv_m(1,1) == Approx(-1.0/2.0).margin(1E-15));
}

TEST_CASE("10: Test Ublas - newton-raphson solver for system of non-linear equations", "[ublas]")
{
  std::cout << "10: Test Ublas - newton-raphson solver for system of non-linear equations" << std::endl;
  ublas::vector<double> guess(2);
  ublas::vector<double> solution(2);

  guess(0) = 1.5;
  guess(1) = 0.1;
  std::cout << guess << "\n";
  auto f = [&](const ublas::vector<double>& u, 
               ublas::vector<double>& res) -> bool
  {
    auto x = u(0);
    auto y = u(1);

    res(0) = x + x*y - 4.0;
    res(1) = x + y - 3.0;

    return true;
  };

  auto df = [&](const ublas::vector<double>& u, 
                ublas::matrix<double>& res) -> bool
  {
    auto x = u(0);
    auto y = u(1);

    res(0,0) = 1.0 + y;
    res(0,1) = x;
    res(1,0) = 1.0;
    res(1,1) = 1.0;

    return true;
  };

  bool solved = newton_raphson(guess, f, df, solution, 1E-15, 1E+5);
  REQUIRE(solved == true);
  REQUIRE(solution(0) == Approx(2.0).margin(1E-15));
  REQUIRE(solution(1) == Approx(1.0).margin(1E-15));
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