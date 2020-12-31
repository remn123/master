#pragma once

#include <vector>
#include <memory>
// #include <boost/numeric/ublas/vector.hpp>

#include <Dummies.h>
#include <DVector.h>
#include <Element.h>
#include <Ghost.h>
#include <Mesh.h>
#include <Node.h>

// namespace ublas = boost::numeric::ublas;

template <typename Method>
class Time
{
public:
  double CFL;
  long MAX_ITER;
  int stages;
  int order;
  double dt;
  long iter;

  std::vector<DVector> Q;
  std::vector<DVector> res;
  DVector res_sum;
  DVector alpha;
  DVector beta;

  Time(double, long);
  Time(double, long, int, int, size_t);
  virtual ~Time();

  void update(std::shared_ptr<Mesh> &, void (*)(std::shared_ptr<Mesh> &));
  void loop(std::shared_ptr<Mesh> &mesh, void (*)(std::shared_ptr<Mesh> &));
  void read_solution(const std::shared_ptr<Mesh> &, size_t);
  void read_residue(const std::shared_ptr<Mesh> &, size_t);
  void write_solution(std::shared_ptr<Mesh> &, size_t);
  double c(int, int);
};

// U(0)   = U(n)
// U(1)   = U(0) + dt*(c10*L0)
// U(2)   = U(0) + dt*(c20*L0 + c21*L1)
// U(3)   = U(0) + dt*(c30*L0 + c31*L1 + c32*L2)
// U(4)   = U(0) + dt*(c40*L0 + c41*L1 + c42*L2 + c43*L3)
// U(5)   = U(0) + dt*(c50*L0 + c51*L1 + c52*L2 + c53*L3 + c54*L4)
// U(n+1) = U(5)