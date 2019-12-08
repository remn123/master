#include <vector>

#include <Mesh.h>
#include <Solver.h>
#include <Helpers.h>


SD::SD(int order, int dimension)
{
  this->order = order;
  this->dim = dimension;

  // Solution Points will be located at Gauss-Legendre Nodes
  Helpers<GL>::init();
  Helpers<GL>::setup(order);
  std::vector<double> nodes_gl = Helpers<GL>::get_nodes();

  /*  
      Second Order SPs
        ________________
       |                |
       |  (2)      (4)  |
       |                |
       |                |
       |  (1)      (3)  |
       |________________|
    */
  for (auto n1 : nodes_gl) // x
  {
    for (auto n2 : nodes_gl) // y
    {
      this->snodes.push_back(Node{n1, n2, 0.0});
    } 
  }

  // Flux Points will be located at Gauss-Legendre-Lobatto Nodes 
  Helpers<GLL>::init();
  Helpers<GLL>::setup(order+1); // Flux must be one order higher
  std::vector<double> nodes_gll = Helpers<GLL>::get_nodes();

  std::vector<Node> vecx;
  std::vector<Node> vecy;

  vecx.resize(order*(order+1));
  vecy.resize(order*(order+1));

  if (dimension==2)
  {
    /*  
      Third Order FPs on x direction
        ________________
       |                |
      (4)     (5)      (6)
       |                |
       |                |
      (1)     (2)      (3)
       |________________|
    */
    for (auto n2 : nodes_gl) // y
    {
      for (auto n1 : nodes_gll) // x
      {
        vecx.push_back(Node{n1, n2, 0.0});
      } 
    }
    this->fnodes.push_back(vecx);

    /*  
      Third Order FPs on y direction

        __(3)_____(6)____
       |                |
       |                |
       |  (2)     (5)   |
       |                |
       |                |
       |__(1)_____(4)___|
    */
    for (auto n1 : nodes_gl) // x
    {
      for (auto n2 : nodes_gll) // y
      {
        vecy.push_back(Node{n1, n2, 0.0});
      } 
    }
    this->fnodes.push_back(vecy);
  }
  else if (dimension==3)
  {
    // NOT IMPLEMENTED FOR THIS LIFE SO SOON
    // std::vector<Node> vecz;  
    // vecz.resize();
    // for (auto n1 : nodes_gl) // x
    // {
    //   for (auto n2 : nodes_gll) // y
    //   {
    //     vecz.push_back(Node{n1, n2, 0.0});
    //   } 
    // }
    // this->fnodes.push_back(vecz);
  }
}

SD::~SD()
{

}

// 0)
void SD::setup(void)
{
  this->calc_sp_nodes();
  this->calc_fp_nodes();
}

// 1)
void SD::boundary_condition (Element& e)
{
  //pass
}

// 2)
void SD::interpolate_sp2fp (Element& e)
{
  //pass
}

// 3)
void SD::riemann_solver (Element& e)
{
  //pass
}

// 4)
void SD::interpolate_fp2sp (Element& e)
{
  //pass
}

// 5)
void SD::residue (Element& e)
{
  //pass
}

// 6)
void SD::solve (Element& e)
{
  this->interpolate_sp2fp(e);
  this->riemann_solver(e);
  this->interpolate_fp2sp(e);
  this->residue(e);
}
