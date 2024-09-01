#pragma once

#include <functional>
#include <memory>
#include <vector>

#include <Dummies.h>
#include <Element.h>
#include <Field.h>
#include <Ghost.h>
#include <Mesh.h>
#include <Node.h>
#include <Weights.h>

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
  std::shared_ptr<Weights> weights;
  //std::vector<std::vector<std::vector<std::vector<double>>>> weights;
  //std::vector<std::vector<std::vector<std::vector<double>>>> d_weights;
  

private:
public:
  SD(int, int);
  SD(int, int, double, double, double);
  ~SD();

  void solve(std::shared_ptr<Mesh> &);
  void to_vtk(const std::shared_ptr<Mesh> &, const std::string &);

  void setup(std::shared_ptr<Mesh> &, 
             std::vector<double> (*)(const Node &) = FIELDS::DEFAULT_FIELD_MAPPING);
  void create_nodes(void);
  void create_weights(void);

  void initialize_element_properties(std::shared_ptr<Element> &,
                                     const std::vector<Vertice> &,
                                     std::vector<double> (*)(const Node &));
  void initialize_ghost_properties(Ghost &);
  void update_edges(std::shared_ptr<Element> &,
                    std::vector<std::shared_ptr<Element>> &,
                    std::vector<Ghost> &);

  void boundary_condition(Ghost &,
                          const std::vector<std::shared_ptr<Element>> &,
                          const std::vector<Vertice> &);
  DVector interpolate_solution_to_node(std::shared_ptr<Element> &, const Node &);
  DVector interpolate_solution_to_fp(std::shared_ptr<Element> &, const Node &, long, int);
  void interpolate_sp2fp(std::shared_ptr<Element> &);
  //void calculate_fluxes(std::shared_ptr<Element>&);
  // void calculate_interface_fluxes(std::shared_ptr<Element>&,
  //                                 const std::vector<std::shared_ptr<Element>>&,
  //                                 const std::vector<Ghost>&);
  void calculate_fluxes_sp(std::shared_ptr<Element> &);
  void calculate_fluxes_fp(std::shared_ptr<Element> &);
  void calculate_fluxes_gfp(Ghost &, const std::vector<std::shared_ptr<Element>> &);
  void riemann_solver(std::shared_ptr<Element> &, const std::vector<std::shared_ptr<Element>> &, const std::vector<Ghost> &);
  void interpolate_fp2sp(std::shared_ptr<Element> &);
  void residue(std::shared_ptr<Element> &);
  void update_fluxes(std::shared_ptr<Element> &, std::vector<std::shared_ptr<Element>> &);
  double get_min_dx(std::shared_ptr<Mesh> &);
  // For Static_Meshes:
  void update_overset(std::shared_ptr<Static_Mesh>&, std::shared_ptr<Static_Mesh>&);
  long mark_fringes_and_holes(std::shared_ptr<Static_Mesh>&, std::shared_ptr<Static_Mesh>&);
  void propagate_holes(std::shared_ptr<Static_Mesh>&, long);
  void update_background_neighboors(std::shared_ptr<Static_Mesh>&, std::shared_ptr<Static_Mesh>&);
  void communicate_data(std::shared_ptr<Static_Mesh>&, const std::shared_ptr<Static_Mesh>&);
  
  void interpolate_and_calculate_flux_sp(std::shared_ptr<Element> &);
  void reconstruct_flux_and_residue(std::shared_ptr<Element> &, const std::vector<std::shared_ptr<Element>> &, const std::vector<Ghost> &);
  Node& project_node(fNode&, Ghost &, std::shared_ptr<Static_Mesh>&, const std::shared_ptr<Static_Mesh>&);
  std::vector<double> get_property_error_at_time(std::shared_ptr<Mesh>&, double);

private:
  void _init_dvec(std::vector<DVector> &, size_t);
  void _init_dvec(std::vector<std::vector<DVector>> &, size_t);
  long _n(const long, const long, const long);
};
