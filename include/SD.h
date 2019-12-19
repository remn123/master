#include <vector>

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
  
  void boundary_condition(Element&);
  void calc_high_order_nodes(Element&);
  void interpolate_sp2fp(Element&);
  void riemann_solver(Element&);
  void interpolate_fp2sp(Element&);
  void residue(Element&);
  void solve(Element&);

private:
  void _init_dvec(std::vector<DVector>&, size_t);
  void _init_dvec(std::vector<std::vector<DVector>>&, size_t);
};
