#include <vector>
#include <memory>

#include <DVector.h>
#include <Mesh.h>


template <typename Equation>
class SD
{
public:
  int order;
  int dimension;
  
  std::vector<std::vector<Node>> fnodes; // flux points (FP)
  std::vector<Node> snodes;              // solution points (SP)
  
private:

public:
  SD(int, int);
  ~SD();

  void setup(Mesh&);
  void create_nodes(void);
  void initialize_properties(Mesh&);
  
  void boundary_condition(std::shared_ptr<Element>&);
  void calc_high_order_nodes(std::shared_ptr<Element>&);
  void interpolate_sp2fp(std::shared_ptr<Element>&);
  void riemann_solver(std::shared_ptr<Element>&);
  void interpolate_fp2sp(std::shared_ptr<Element>&);
  void residue(std::shared_ptr<Element>&);
  void solve(std::shared_ptr<Element>&);

private:
  void _init_dvec(std::vector<DVector>&, size_t);
  void _init_dvec(std::vector<std::vector<DVector>>&, size_t);
};
