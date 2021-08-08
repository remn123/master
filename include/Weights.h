#pragma once

#include <DVector.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::matrix<double> MatrixDouble;
typedef boost::numeric::ublas::vector<double> uvector;

class Weights
{
private:
  std::size_t sp;
  std::size_t fp;
  MatrixDouble Mxfp;
  MatrixDouble Myfp;
  MatrixDouble Mxsp;
  MatrixDouble Mysp;

public:
  Weights(int, int);
  ~Weights();

  std::vector<DVector> prod(int , std::vector<DVector>&);
  std::vector<DVector> prod(int , std::vector<DVector>&, std::vector<double>&);
  void set(int, uint, uint, double);
private:
  uvector _prod(int, uvector);
  std::size_t _get_nrows(int);
};