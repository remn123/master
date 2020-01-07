#include <vector>
#include <memory>

#include <Dummies.h>
#include <Poly.h>
#include <Helpers.h>
#include <Mesh.h>
#include <SD.h>


// Constructor for Euler Equations
template <>
SD<Euler>::SD(int order, int dimension)
{
  this->order = order;
  this->dimension = dimension;
  std::cout << "Initializing Euler SD solver...\n";
}

// Constructor for Navier-Stokes Equations
template <>
SD<NavierStokes>::SD(int order, int dimension)
{
  this->order = order;
  this->dimension = dimension;
  std::cout << "Initializing Navier-Stokes SD solver...\n";
}

template <typename Equation>
SD<Equation>::~SD()
{

}

// 0.1)
template <typename Equation>
void SD<Equation>::create_nodes(void)
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

// Auxiliar function
template <typename Equation>
void SD<Equation>::_init_dvec(std::vector<DVector>& vec, size_t num_nodes)
{
  std::vector<double> init_vec(this->dimension+2, 1.0); 

  vec.clear();
  vec.resize(num_nodes);
  for (auto& dvec : vec)
  {
    dvec = init_vec; 
  }
}

template <typename Equation>
void SD<Equation>::_init_dvec(std::vector<std::vector<DVector>>& vec, size_t num_nodes)
{
  std::vector<double> init_vec(this->dimension+2, 1.0); 
  
  vec.clear();
  vec.resize(this->dimension);
  for (auto& dirvec : vec) // x, y, ...
  {
    dirvec.clear();
    dirvec.resize(num_nodes); // number of nodes
    for (auto& dvec : dirvec)
    {
      dvec = init_vec; 
    }
  }
}

// 0.2) Initialize property template
template <>
void SD<Euler>::initialize_properties(Mesh& mesh)
{
  for (auto& e : mesh.elems)
  {
    // Conservative Properties
    this->_init_dvec(e->Qsp,  this->snodes.size());
    this->_init_dvec(e->Qfp,  this->fnodes[0].size());

    // Convective Fluxes
    this->_init_dvec(e->Fcsp, this->snodes.size());
    this->_init_dvec(e->Fcfp, this->fnodes[0].size());

    // Gradients
    this->_init_dvec(e->dFcsp, this->snodes.size());
    this->_init_dvec(e->dFcfp, this->fnodes[0].size());
  }
}

template <>
void SD<NavierStokes>::initialize_properties(Mesh& mesh)
{
  for (auto& e : mesh.elems)
  {
    // Conservative Properties
    this->_init_dvec(e->Qsp,  this->snodes.size());
    this->_init_dvec(e->Qfp,  this->fnodes[0].size());

    // Convective Fluxes
    this->_init_dvec(e->Fcsp, this->snodes.size());
    this->_init_dvec(e->Fcfp, this->fnodes[0].size());

    // Diffusive Fluxes
    this->_init_dvec(e->Fdsp, this->snodes.size());
    this->_init_dvec(e->Fdfp, this->fnodes[0].size());

    // Gradients
    // Conservative Properties
    this->_init_dvec(e->dQsp,  this->snodes.size());
    this->_init_dvec(e->dQfp,  this->fnodes[0].size());

    // Convective Fluxes
    this->_init_dvec(e->dFcsp, this->snodes.size());
    this->_init_dvec(e->dFcfp, this->fnodes[0].size());

    // Diffusive Fluxes
    this->_init_dvec(e->dFdsp, this->snodes.size());
    this->_init_dvec(e->dFdfp, this->fnodes[0].size());

    // Residue
    this->_init_dvec(e->res, this->snodes.size());
  }
}

// 0)
template <typename Equation>
void SD<Equation>::setup(Mesh& mesh)
{
  /* 
    On setup method, all solution and flux nodes
    are allocated and calculated in create_nodes.

    Then initialize_properties will allocate and
    populate the mesh with an initial value
  */
  this->create_nodes();
  this->initialize_properties(mesh);
}


// 1) BOUNDARY CONDITIONS
/* 
  Types of Boundary Conditions:
    - Dirichlet BC
    - Neumann BC

  Boundary Conditions:
    - Solid Inviscid Wall (Euler):
        Un - Uwall,n = (U-Uwall).n = 0
        Uwall,t = some value (SLIP - there's no boundary layer here)

    - Solid Viscous Wall (NavierStokes):
        Un - Uwall,n = (U-Uwall).n = 0
        Uwall,t = 0 (NO-SLIP - there's a boundary layer)

    - Inlet:
        Q = constant at specific place

    - Outlet (non-reflexive):
        Q_L = Q_R

    - Periodic:
        Q_inflow = Q_outflow

*/
template <typename Equation>
void SD<Equation>::boundary_condition (std::shared_ptr<Element>& e)
{
  
}

// 2) INTERPOLATE FROM SOLUTION POINTS TO FLUX POINTS
template <typename Equation>
void SD<Equation>::interpolate_sp2fp (std::shared_ptr<Element>& e)
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

// 3) RIEMANN SOLVER
template <>
void SD<Euler>::riemann_solver (std::shared_ptr<Element>& e)
{
  //pass
}

template <>
void SD<NavierStokes>::riemann_solver (std::shared_ptr<Element>& e)
{
  //pass
}

// 4) INTERPOLATE FROM FLUX POINTS TO SOLUTION POINTS
// Euler
template<>
void SD<Euler>::interpolate_fp2sp (std::shared_ptr<Element>& e)
{
  Helpers<GLL>::init();
  Helpers<GLL>::set_nodes(this->order+1);
  
  Helpers<Lagrange>::init();
  Helpers<Lagrange>::set_nodes(Helpers<GLL>::get_nodes());
  
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
      e.dFcsp[index-1][s_index-1] = 0.0;
  
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
        e.dFcsp[index-1][s_index-1] += (dLcsi*dLeta)*e.Fcfp[index-1][f_index-1];
      }
    }
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GLL>::delete_nodes();
}

template<>
void SD<NavierStokes>::interpolate_fp2sp (std::shared_ptr<Element>& e)
{
  Helpers<GLL>::init();
  Helpers<GLL>::set_nodes(this->order+1);
  
  Helpers<Lagrange>::init();
  Helpers<Lagrange>::set_nodes(Helpers<GLL>::get_nodes());
  
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
      e.dFcsp[index-1][s_index-1] = 0.0;
      e.dFdsp[index-1][s_index-1] = 0.0;
  
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
        e.dFcsp[index-1][s_index-1] += (dLcsi*dLeta)*e.Fcfp[index-1][f_index-1];
        e.dFdsp[index-1][s_index-1] += (dLcsi*dLeta)*e.Fdfp[index-1][f_index-1];
      }
    }
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GLL>::delete_nodes();
}

// 5) RESIDUE 
// Euler
template <>
void SD<Euler> ::residue(std::shared_ptr<Element>& e)
{
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
      e.res[s_index-1] += -e.dFcsp[index-1][s_index-1];
    }
  }
}

// Navier-Stokes
template <>
void SD<NavierStokes>::residue (std::shared_ptr<Element>& e)
{
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
      e.res[s_index-1] += -e.dFcsp[index-1][s_index-1] + e.dFdsp[index-1][s_index-1];
    }
  }
}

// 6) SOLVE
template <typename Equation>
void SD<Equation>::solve (std::shared_ptr<Element>& e)
{
  this->bondary_condition(e);
  this->interpolate_sp2fp(e);
  this->calculate_fluxes(e);
  this->riemann_solver(e);
  this->interpolate_fp2sp(e);
  this->residue(e);
}



/*
  1) Setup (all element in Mesh)
     1.1) Calculate solution and fluxes points
     1.2) Initialize solution and fluxes

  2) Solver Loop (for each element in Mesh)
     2.1) Apply Boundary Conditions on interfaces (flux points)
     2.2) Interpolate Q and dQ from SPs to FPs
     2.3) Calculate Fluxes at internal FPs
     2.4) Calculate Fluxes at faces FPs (riemann solver and BR2)
     2.5) Interpolate Fluxes derivatives from FPs to SPs
     2.6) Calculate Residue
  
  3) Iterate in time
     3.1) Calculate residue norm
     3.2) Check if it's already converged
     3.3) (if not) Apply time iteration then go to (2)
*/

