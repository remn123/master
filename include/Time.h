#pragma once

#include <vector>
#include <memory>
#include <boost/numeric/ublas/vector.hpp>

#include <Dummies.h>
#include <Element.h>
#include <Ghost.h>
#include <Mesh.h>
#include <Node.h>

namespace ublas = boost::numeric::ublas;

template <typename Method>
class Time
{
public:
  double CFL;
  long MAX_ITER;
  int stages;
  int order;
  double dt;
  
  ublas::vector<double> Q;
  ublas::vector<double> res_sum;
  ublas::vector<double> Qnew;


  Time(double, long);
  Time(double, long, int, int);
  virtual ~Time();

  void update(std::shared_ptr<Mesh>&, void*);
  void read_from_mesh(const std::shared_ptr<Mesh>&);
  void write_into_mesh(std::shared_ptr<Mesh>&);
  void get_residue(const std::shared_ptr<Mesh>&);
};

