#pragma once

#include <functional>
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
  double tinf;

  
  std::vector<DVector> Q;
  std::vector<DVector> res;
  DVector res_sum;
  DVector alpha;
  DVector beta;
  std::string log_path;

  Time(double, long);
  Time(double, long, int, int, size_t, int, double, double, const std::string);
  virtual ~Time();

  void update(
    std::shared_ptr<Mesh> &, 
    std::function<void(std::shared_ptr<Mesh> &)>,
    std::function<std::vector<double>(std::shared_ptr<Mesh> &, double)>
  );
  void update(std::shared_ptr<Static_Mesh> &,
              std::shared_ptr<Static_Mesh> &,
              std::function<void(std::shared_ptr<Mesh> &)>,
              std::function<std::vector<double>(std::shared_ptr<Mesh> &, double)>,
              std::function<void(std::shared_ptr<Static_Mesh> &, 
                                 const std::shared_ptr<Static_Mesh> &)>);
  void loop(
    std::shared_ptr<Mesh> &, 
    std::function<void(std::shared_ptr<Mesh> &)>,
    std::function<std::vector<double>(std::shared_ptr<Mesh> &, double)>,
    const std::string &,
    std::function<void(const std::shared_ptr<Mesh> &, const std::string &)>
  );
  void loop(
    std::shared_ptr<Static_Mesh> &,
    std::shared_ptr<Static_Mesh> &, 
    std::function<void(std::shared_ptr<Mesh> &)>,
    std::function<std::vector<double>(std::shared_ptr<Mesh> &, double)>,
    std::function<void(std::shared_ptr<Static_Mesh> &, 
                       const std::shared_ptr<Static_Mesh> &)>,
    const std::string &,
    const std::string &,
    std::function<void(const std::shared_ptr<Mesh> &, const std::string &)>
  );
  void read_solution(const std::shared_ptr<Mesh> &, size_t);
  void read_residue(const std::shared_ptr<Mesh> &, size_t);
  void write_solution(std::shared_ptr<Mesh> &, size_t);

  void read_solution(const std::shared_ptr<Static_Mesh> &,
                     const std::shared_ptr<Static_Mesh> &,
                     size_t);
  void read_residue(const std::shared_ptr<Static_Mesh> &,
                    const std::shared_ptr<Static_Mesh> &,
                    size_t);
  void write_solution(std::shared_ptr<Static_Mesh> &,
                      std::shared_ptr<Static_Mesh> &,
                      size_t);
  void save(const std::shared_ptr<Mesh> &,
            const std::string &,
            std::function<void(const std::shared_ptr<Mesh> &, const std::string &)>);
  void save_log(const std::string &);
  double c(int, int);
};

// U(0)   = U(n)
// U(1)   = U(0) + dt*(c10*L0)
// U(2)   = U(0) + dt*(c20*L0 + c21*L1)
// U(3)   = U(0) + dt*(c30*L0 + c31*L1 + c32*L2)
// U(4)   = U(0) + dt*(c40*L0 + c41*L1 + c42*L2 + c43*L3)
// U(5)   = U(0) + dt*(c50*L0 + c51*L1 + c52*L2 + c53*L3 + c54*L4)
// U(n+1) = U(5)