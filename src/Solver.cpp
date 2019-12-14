#include <vector>

#include <Helpers.h>
#include <Mesh.h>
#include <Solver.h>


SD::SD(int order, int dimension)
{
  this->order = order;
  this->dim = dimension;
}

SD::~SD()
{

}

// 0)
void SD::setup(void)
{
  // Solution Points will be located at Gauss-Legendre Nodes
  Helpers<GL>::init();
  Helpers<GL>::setup(order);
  std::vector<double> nodes_gl = Helpers<GL>::get_nodes();

  Helpers<GL>::delete_nodes();

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
  
  Helpers<GLL>::delete_nodes();

  std::vector<Node> vecx;
  std::vector<Node> vecy;

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

// 1)
void SD::boundary_condition (Element& e)
{
  //pass
}

// 2)
void SD::interpolate_sp2fp (Element& e)
{
  Helpers<GL>::init();
  Helpers<GL>::setup(this->order);
  
  Helpers<Lagrange>::init();
  Helpers<Lagrange>::setup(Helpers<GL>::get_nodes());
  
  double csi=0.0, eta=0.0;
  double Lcsi, Leta;
  unsigned int i, j;
  unsigned int index, s_index, f_index;
  
  index = 0;
  for (auto& vec_lines : this->fnodes)
  {
    index++;

    f_index = 0;
    for (auto& node : vec_lines) // line nodes for a specific direction (x, y)
    {
      f_index++;
      // flux node coordinates
      csi = node.coords[0]; 
      eta = node.coords[1];

      // initialize flux nodes solution
      e.Qfp[index-1][f_index-1] = 0.0;

      s_index = 0;
      for (auto& n : this->snodes)
      {
        s_index++;
        // s_index = (j+1) + this->order*i;

        i = (int) s_index / this->order;
        j = s_index % this->order;

        Lcsi = Helpers<Lagrange>::Pn(i, csi);
        Leta = Helpers<Lagrange>::Pn(j, eta);
        e.Qfp[index-1][f_index-1] += (Lcsi*Leta)*e.Qsp[s_index-1];
      }
    }
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GL>::delete_nodes();
}

// 3)
void SD::riemann_solver (Element& e)
{
  //pass
}

// 4)
void SD::interpolate_fp2sp (Element& e)
{
  Helpers<GLL>::init();
  Helpers<GLL>::setup(this->order+1);
  
  Helpers<Lagrange>::init();
  Helpers<Lagrange>::setup(Helpers<GLL>::get_nodes());
  
  double csi=0.0, eta=0.0;
  //double Lcsi, Leta;
  double dLcsi, dLeta;
  unsigned int i, j;
  unsigned int index, s_index, f_index;
  
  s_index = 0;
  for (auto& node : this->snodes)
  {
    // for each solution point
    // I will calculate the lagrange polynomial at its position
    // interpolate the flux in (x/y) direction from fps
    s_index++;

    csi = node.coords[0];
    eta = node.coords[1];

    index = 0;
    for (auto& vec_lines : this->fnodes)
    {
      index++;

      //e.Fsp[index-1][s_index-1] = 0.0;
      e.dFsp[index-1][s_index-1] = 0.0;
  
      f_index = 0;
      for (auto& node : vec_lines)
      {
        f_index++;      
        i = (int) f_index / this->order;
        j = f_index % this->order;

        //Lcsi = Helpers<Lagrange>::Pn(i, csi);
        //Leta = Helpers<Lagrange>::Pn(j, eta);
        dLcsi = Helpers<Lagrange>::dPn(i, csi);
        dLeta = Helpers<Lagrange>::dPn(j, eta);
        
        // index-1 is related to the flux direction (x/0 or y/1)
        //e.Fsp[index-1][s_index-1] += Lcsi*Leta*e.Ffp[index-1][f_index-1];
        e.dFsp[index-1][s_index-1] += (dLcsi*dLeta)*e.Ffp[index-1][f_index-1];
      }
    }
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GLL>::delete_nodes();
}

// 5)
void SD::residue (Element& e)
{
  double csi=0.0, eta=0.0;
  //double Lcsi, Leta;
  double dLcsi, dLeta;
  unsigned int index, s_index;

  s_index = 0;
  for (auto& node : this->snodes)
  {
    s_index++;

    e.res[s_index-1] = 0.0;
    index = 0;
    for (auto& vec_lines : this->fnodes)
    {
      index++;
      e.res[s_index-1] += -e.dFcsp[index-1][s_index-1] + e.d2Fdsp[index-1][s_index-1];
    }
  }
  
}

// 6)
void SD::solve (Element& e)
{
  this->bondary_condition(e);
  this->interpolate_sp2fp(e);
  this->riemann_solver(e);
  this->interpolate_fp2sp(e);
  this->residue(e);
}
