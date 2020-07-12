#pragma once

#include <vector>
#include <memory>

#include <Dummies.h>
#include <Element.h>
#include <Ghost.h>
#include <Mesh.h>
#include <Node.h>

template <typename Equation>
class SD
{
public:
  int order;
  int dimension;
  double M;
  double Re;
  double MU;
  double MUv;
  double GAMMA;
  double KAPPA;

  std::vector<std::vector<Node>> fnodes; // flux points (FP)
  std::vector<Node> snodes;              // solution points (SP)

private:
public:
  SD(int, int);
  SD(int, int, double, double, double);
  ~SD();

  void setup(std::shared_ptr<Mesh> &);
  void create_nodes(void);

  void initialize_properties(std::shared_ptr<Element> &,
                             const std::vector<Vertice> &);
  void initialize_properties(Ghost &);
  void update_edges(std::shared_ptr<Element> &,
                    std::vector<std::shared_ptr<Element>> &,
                    std::vector<Ghost> &);

  void boundary_condition(Ghost &,
                          const std::vector<std::shared_ptr<Element>> &);
  void interpolate_interface(Mesh &, std::shared_ptr<Element> &);
  void interpolate_sp2fp(std::shared_ptr<Element> &);
  //void calculate_fluxes(std::shared_ptr<Element>&);
  // void calculate_interface_fluxes(std::shared_ptr<Element>&,
  //                                 const std::vector<std::shared_ptr<Element>>&,
  //                                 const std::vector<Ghost>&);
  void calculate_fluxes_sp(std::shared_ptr<Element> &);
  void calculate_fluxes_fp(std::shared_ptr<Element> &);
  void calculate_fluxes_fp(Ghost &, const std::vector<std::shared_ptr<Element>> &);
  void riemann_solver(std::shared_ptr<Element> &, const std::vector<std::shared_ptr<Element>> &, const std::vector<Ghost> &);
  void interpolate_fp2sp(std::shared_ptr<Element> &);
  void residue(std::shared_ptr<Element> &);
  void solve(std::shared_ptr<Mesh> &);
  void to_vtk(const std::shared_ptr<Mesh> &, const std::string &);

private:
  void _init_dvec(std::vector<DVector> &, size_t);
  void _init_dvec(std::vector<std::vector<DVector>> &, size_t);
  long _n(const long, const long, const long);
};
