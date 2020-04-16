#include <vector>
#include <memory>

#include <Dummies.h>
#include <Poly.h>
#include <Helpers.h>
#include <Mesh.h>
#include <SD.h>

// explicit instantiations
template class SD<Euler>;
template class SD<NavierStokes>;


// Constructor for Euler Equations
template <>
SD<Euler>::SD(int order, int dimension)
{
  this->order = order;
  this->dimension = dimension;
  this->MU = 0.0;
  this->MUv = 0.0;
  std::size_t snodes_size = order*order;
  std::size_t fnodes_size = (order+1)*order;
  std::size_t dim = dimension;

  this->snodes = std::vector<Node>{snodes_size};
  this->fnodes = std::vector<std::vector<Node>>{dim, std::vector<Node>{fnodes_size}};    

  std::cout << "Initializing Euler SD solver...\n";
}

// Constructor for Navier-Stokes Equations
template <>
SD<NavierStokes>::SD(int order, int dimension, double mu, double Mach, double Reynolds)
{
  this->order = order;
  this->dimension = dimension;
  this->M = Mach;
  this->Re = Reynolds;
  this->MU = mu;
  this->MUv = (-2.0/3.0)*mu;
  std::size_t snodes_size = order*order;
  std::size_t fnodes_size = (order+1)*order;
  std::size_t dim = dimension;

  this->snodes = std::vector<Node>{snodes_size};
  this->fnodes = std::vector<std::vector<Node>>{dim, std::vector<Node>{fnodes_size}};
  
  std::cout << "Initializing Navier-Stokes SD solver...\n";
}

template <typename Equation>
SD<Equation>::~SD()
{

}

template <>
SD<Euler>::~SD()
{

}

// 0.1)
template <typename Equation>
void SD<Equation>::create_nodes(void)
{
   // Solution Points will be located at Gauss-Legendre Nodes
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(order);
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
  std::size_t s_index=0;
  for (auto n1 : nodes_gl) // x
  {
    for (auto n2 : nodes_gl) // y
    { 
      this->snodes[s_index] = Node{n1, n2, 0.0};
      s_index++;
    } 
  }

  // Flux Points will be located at Gauss-Legendre-Lobatto Nodes 
  Helpers<GLL>::init();
  Helpers<GLL>::set_nodes(order+1); // Flux must be one order higher
  std::vector<double> nodes_gll = Helpers<GLL>::get_nodes();
  
  Helpers<GLL>::delete_nodes();

  std::vector<Node> vecx;
  std::vector<Node> vecy;

  std::size_t f_index=0;

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
        this->fnodes[0][f_index] = Node{n1, n2, 0.0};
        f_index++;
      } 
    }

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
    f_index=0;
    for (auto n1 : nodes_gl) // x
    {
      for (auto n2 : nodes_gll) // y
      {
        this->fnodes[1][f_index] = Node{n1, n2, 0.0};
        f_index++;
      } 
    }
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
// INITIALIZE GHOST PROPERTIES
template <>
void SD<Euler>::initialize_properties(Ghost& g)
{
  // PHYSICAL
  // Conservative Properties
  this->_init_dvec(g.physical->Qsp,  this->snodes.size());
  this->_init_dvec(g.physical->Qfp,  this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(g.physical->Fcsp, this->snodes.size());
  this->_init_dvec(g.physical->Fcfp, this->fnodes[0].size());

  // Gradients
  // Convective Fluxes
  this->_init_dvec(g.physical->dFcsp, this->snodes.size());
  this->_init_dvec(g.physical->dFcfp, this->fnodes[0].size());

  // COMPUTATIONAL
  // Conservative Properties
  this->_init_dvec(g.computational->Qsp, this->snodes.size());
  this->_init_dvec(g.computational->Qfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(g.computational->Fcsp, this->snodes.size());
  this->_init_dvec(g.computational->Fcfp, this->fnodes[0].size());

  // Gradients
  // Convective Fluxes
  this->_init_dvec(g.computational->dFcsp, this->snodes.size());
  this->_init_dvec(g.computational->dFcfp, this->fnodes[0].size());


  // Residue
  this->_init_dvec(g.computational->res, this->snodes.size());
}

template <>
void SD<NavierStokes>::initialize_properties(Ghost& g)
{
  // PHYSICAL
  // Conservative Properties
  this->_init_dvec(g.physical->Qsp,  this->snodes.size());
  this->_init_dvec(g.physical->Qfp,  this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(g.physical->Fcsp, this->snodes.size());
  this->_init_dvec(g.physical->Fcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(g.physical->Fdsp, this->snodes.size());
  this->_init_dvec(g.physical->Fdfp, this->fnodes[0].size());

  // Gradients
  // Conservative Properties
  this->_init_dvec(g.physical->dQsp,  this->snodes.size());
  this->_init_dvec(g.physical->dQfp,  this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(g.physical->dFcsp, this->snodes.size());
  this->_init_dvec(g.physical->dFcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(g.physical->dFdsp, this->snodes.size());
  this->_init_dvec(g.physical->dFdfp, this->fnodes[0].size());

  // COMPUTATIONAL
  // Conservative Properties
  this->_init_dvec(g.computational->Qsp, this->snodes.size());
  this->_init_dvec(g.computational->Qfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(g.computational->Fcsp, this->snodes.size());
  this->_init_dvec(g.computational->Fcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(g.computational->Fdsp, this->snodes.size());
  this->_init_dvec(g.computational->Fdfp, this->fnodes[0].size());

  // Gradients
  // Conservative Properties
  this->_init_dvec(g.computational->dQsp, this->snodes.size());
  this->_init_dvec(g.computational->dQfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(g.computational->dFcsp, this->snodes.size());
  this->_init_dvec(g.computational->dFcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(g.computational->dFdsp, this->snodes.size());
  this->_init_dvec(g.computational->dFdfp, this->fnodes[0].size());

  // Residue
  this->_init_dvec(g.computational->res, this->snodes.size());
}


// INITIALIZE ELEMENTS PROPERTIES
template <>
void SD<Euler>::initialize_properties(std::shared_ptr<Element>& e, const std::vector<Vertice>& enodes)
{
  // PHYSICAL
  // Conservative Properties
  this->_init_dvec(e->physical->Qsp,  this->snodes.size());
  this->_init_dvec(e->physical->Qfp,  this->fnodes[0].size());
  
  // Convective Fluxes
  this->_init_dvec(e->physical->Fcsp, this->snodes.size());
  this->_init_dvec(e->physical->Fcfp, this->fnodes[0].size());

  // Gradients
  // Convective Fluxes
  this->_init_dvec(e->physical->dFcsp, this->snodes.size());
  this->_init_dvec(e->physical->dFcfp, this->fnodes[0].size());

  for (auto& ed: e->edges)
  {
    this->_init_dvec(ed.physical->Qfp,   this->order);
    this->_init_dvec(ed.physical->Fcfp,  this->order);
    this->_init_dvec(ed.physical->dFcfp, this->order);
  }

  // Metrics 
  e->allocate_jacobian(this->order);
  e->calculate_jacobian(this->snodes, this->fnodes, enodes);

  // COMPUTATIONAL
  // Conservative Properties
  this->_init_dvec(e->computational->Qsp, this->snodes.size());
  this->_init_dvec(e->computational->Qfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(e->computational->Fcsp, this->snodes.size());
  this->_init_dvec(e->computational->Fcfp, this->fnodes[0].size());

  // Gradients
  // Convective Fluxes
  this->_init_dvec(e->computational->dFcsp, this->snodes.size());
  this->_init_dvec(e->computational->dFcfp, this->fnodes[0].size());


  // Residue
  this->_init_dvec(e->computational->res, this->snodes.size());
}

template <>
void SD<NavierStokes>::initialize_properties(std::shared_ptr<Element>& e, const std::vector<Vertice>& enodes)
{
  // PHYSICAL
  // Conservative Properties
  this->_init_dvec(e->physical->Qsp,  this->snodes.size());
  this->_init_dvec(e->physical->Qfp,  this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(e->physical->Fcsp, this->snodes.size());
  this->_init_dvec(e->physical->Fcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(e->physical->Fdsp, this->snodes.size());
  this->_init_dvec(e->physical->Fdfp, this->fnodes[0].size());

  // Gradients
  // Conservative Properties
  this->_init_dvec(e->physical->dQsp,  this->snodes.size());
  this->_init_dvec(e->physical->dQfp,  this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(e->physical->dFcsp, this->snodes.size());
  this->_init_dvec(e->physical->dFcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(e->physical->dFdsp, this->snodes.size());
  this->_init_dvec(e->physical->dFdfp, this->fnodes[0].size());

  // Metrics 
  e->allocate_jacobian(this->order);
  e->calculate_jacobian(this->snodes, this->fnodes, enodes);

  for (auto& ed: e->edges)
  {
    this->_init_dvec(ed.physical->Qfp,   this->order);
    this->_init_dvec(ed.physical->dQfp,  this->order);
    this->_init_dvec(ed.physical->Fcfp,  this->order);
    this->_init_dvec(ed.physical->dFcfp, this->order);
    this->_init_dvec(ed.physical->Fdfp,  this->order);
    this->_init_dvec(ed.physical->dFdfp, this->order);
  }

  // COMPUTATIONAL
  // Conservative Properties
  this->_init_dvec(e->computational->Qsp, this->snodes.size());
  this->_init_dvec(e->computational->Qfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(e->computational->Fcsp, this->snodes.size());
  this->_init_dvec(e->computational->Fcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(e->computational->Fdsp, this->snodes.size());
  this->_init_dvec(e->computational->Fdfp, this->fnodes[0].size());

  // Gradients
  // Conservative Properties
  this->_init_dvec(e->computational->dQsp, this->snodes.size());
  this->_init_dvec(e->computational->dQfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(e->computational->dFcsp, this->snodes.size());
  this->_init_dvec(e->computational->dFcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(e->computational->dFdsp, this->snodes.size());
  this->_init_dvec(e->computational->dFdfp, this->fnodes[0].size());

  // Residue
  this->_init_dvec(e->computational->res, this->snodes.size());
}

// 0)
template <typename Equation>
void SD<Equation>::setup(std::shared_ptr<Mesh>& mesh)
{
  /* 
    On setup method, all solution and flux nodes
    are allocated and calculated in create_nodes.

    Then initialize_properties will allocate and
    populate the mesh with an initial value
  */
  this->create_nodes();
  for (auto& e : mesh->elems)
    this->initialize_properties(e, mesh->nodes);
  for (auto& g : mesh->ghosts)
    this->initialize_properties(g);
}


// 1) BOUNDARY CONDITIONS
/* 
  Types of Boundary Conditions:
    - Dirichlet BC
    - Neumann BC

  Boundary Conditions:
    - Solid Inviscid Wall (Euler) [TYPE 0]:
        Un - Uwall,n = (U-Uwall).n = 0
        Uwall,t = some value (SLIP - there's no boundary layer here)

    - Solid Viscous Wall (NavierStokes) [TYPE 4]:
        Un - Uwall,n = (U-Uwall).n = 0
        Uwall,t = 0 (NO-SLIP - there's a boundary layer)

    - Inlet [TYPE 1]:
        Q = constant at specific place

    - Outlet (non-reflexive) [TYPE 2]:
        Q_L = Q_R

    - Periodic [TYPE 3]:
        Q_inflow = Q_outflow

*/
template <typename Equation>
void SD<Equation>::boundary_condition (Ghost& g)
{
  // switch (g.type)
  // {
  //   //- Solid Inviscid Wall (Euler) [TYPE 0]:
  //   case 0:
  //     //  Un - Uwall,n = (U-Uwall).n = 0
  //     g.physical->Qfp[0]
  //     //  Uwall,t = some value (SLIP - there's no boundary layer here)


  // }
}

// Interpolation from:
template <typename Equation>
void SD<Equation>::interpolate_interface (Mesh& mesh, std::shared_ptr<Element>& e)
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(this->order);
  
  Helpers<Lagrange>::init();
  Helpers<Lagrange>::set_nodes(Helpers<GL>::get_nodes());

  double x=0.0, y=0.0; 
  double csi=0.0, eta=0.0;
  double Lcsi, Leta;
  unsigned int i, j;
  unsigned int index, s_index, f_index;
  
  // node coordinates
  
  for (auto& ed: e->edges)
  {
    for ()
    {

    }
    x = node.coords[0]; 
    y = node.coords[1];
  }


  // initialize flux nodes solution
  e->computational->Qfp[f_index-1] = 0.0;
  e->physical->Qfp[f_index-1] = 0.0;

  s_index = 0;
  for (auto& n : this->snodes)
  {
    s_index++;
    // s_index = (j+1) + this->order*i;

    i = (int) s_index / this->order;
    j = s_index % this->order;

    Lcsi = Helpers<Lagrange>::Pn(i, csi);
    Leta = Helpers<Lagrange>::Pn(j, eta);
    e->computational->Qfp[f_index-1] += ((Lcsi*Leta)*e->computational->Qsp[s_index-1]);
    e->physical->Qfp[f_index-1] += ((Lcsi*Leta)*e->physical->Qsp[s_index-1]);
  } 


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
      e->computational->Qfp[f_index-1] = 0.0;
      e->physical->Qfp[f_index-1] = 0.0;

      s_index = 0;
      for (auto& n : this->snodes)
      {
        s_index++;
        // s_index = (j+1) + this->order*i;

        i = (int) s_index / this->order;
        j = s_index % this->order;

        Lcsi = Helpers<Lagrange>::Pn(i, csi);
        Leta = Helpers<Lagrange>::Pn(j, eta);
        e->computational->Qfp[f_index-1] += ((Lcsi*Leta)*e->computational->Qsp[s_index-1]);
        e->physical->Qfp[f_index-1] += ((Lcsi*Leta)*e->physical->Qsp[s_index-1]);
      }
    }
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GL>::delete_nodes();
}

// ---------------------- //



// 2) INTERPOLATE FROM SOLUTION POINTS TO FLUX POINTS
template <typename Equation>
void SD<Equation>::interpolate_sp2fp (std::shared_ptr<Element>& e)
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(this->order);
  
  Helpers<Lagrange>::init();
  Helpers<Lagrange>::set_nodes(Helpers<GL>::get_nodes());
  
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
      e->computational->Qfp[f_index-1] = 0.0;
      e->physical->Qfp[f_index-1] = 0.0;

      s_index = 0;
      for (auto& n : this->snodes)
      {
        s_index++;
        // s_index = (j+1) + this->order*i;

        i = (int) s_index / this->order;
        j = s_index % this->order;

        Lcsi = Helpers<Lagrange>::Pn(i, csi);
        Leta = Helpers<Lagrange>::Pn(j, eta);
        e->computational->Qfp[f_index-1] += ((Lcsi*Leta)*e->computational->Qsp[s_index-1]);
        e->physical->Qfp[f_index-1] += ((Lcsi*Leta)*e->physical->Qsp[s_index-1]);
      }
    }
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GL>::delete_nodes();
}

// 3) Calculate Fluxes
// 3.1) Calculate Fluxes at SPs
// Euler
template <>
void SD<Euler>::calculate_fluxes_sp (std::shared_ptr<Element>& e)
{
  double q1, q2, q3, q4, q5;

  unsigned int s_index;
  
  s_index = 0;
  // for each solution point, calculate the flux vector
  for (auto& node : this->snodes)
  {
    s_index++;

    // these Q must be the physical solution, so I divide the computational solution Q by |J|
    q1 = e->physical->Qsp[s_index-1][0]/e->J; 
    q2 = e->physical->Qsp[s_index-1][1]/e->J;
    q3 = e->physical->Qsp[s_index-1][2]/e->J;
    q4 = e->physical->Qsp[s_index-1][3]/e->J;

    if (this->dimension+2 == 4)
    {
      /*
        Ji = [a b; = [dcsi_dx  dcsi_dy; =  [ (1,1)   (1,2);
              c d]    deta_dx  deta_dy]      (2,1)   (2,2)]

        Fc = [a b; c d].[Fxc; Fyc] = [a*Fxc + b*Fyc; c*Fxc + d*Fyc]
        Fc = [dcsi_dx*Fxc + dcsi_dy*Fyc; 
              deta_dx*Fxc + deta_dy*Fyc]

        Fcsic = dcsi_dx*Fxc + dcsi_dy*Fyc
        Fetac = deta_dx*Fxc + deta_dy*Fyc

        Fcsic = e->Fc[0]
        Fetac = e->Fc[1]
      */
      e->computational->Fcsp[0][s_index-1] = {e->Ji[1][1][s_index]*q2 + e->Ji[1][2][s_index]*q3, 
                                              e->Ji[1][1][s_index]*(q2*q2/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1)) + e->Ji[1][2][s_index]*q3*q2/q1,
                                              e->Ji[1][1][s_index]*(q2*q3/q1) + e->Ji[1][2][s_index]*(q3*q3/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1)),
                                              e->Ji[1][1][s_index]*((q2/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1)) + e->Ji[1][2][s_index]*((q3/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1))};
                              
      e->computational->Fcsp[1][s_index-1] = {e->Ji[2][1][s_index]*q2 + e->Ji[2][2][s_index]*q3, 
                                              e->Ji[2][1][s_index] * (q2 * q2 / q1 + (this->GAMMA - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + e->Ji[2][2][s_index] * q3 * q2 / q1,
                                              e->Ji[2][1][s_index] * (q2 * q3 / q1) + e->Ji[2][2][s_index] * (q3 * q3 / q1 + (this->GAMMA - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
                                              e->Ji[2][1][s_index]*((q2/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1)) + e->Ji[2][2][s_index]*((q3/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1))};
      
      // e->Fcsp[0][s_index-1] = {q2, 
      //                          q2*q2/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1),
      //                          q2*q3/q1,
      //                          (q2/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1)};
      
      // e->Fcsp[1][s_index-1] = {q3, 
      //                          q3*q2/q1,
      //                          q3*q3/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1),
      //                          (q3/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1)};
    }
    else if (this->dimension+2 == 5)
    {
      q5 = e->physical->Qsp[s_index-1][4];

       /*
        Ji = [a b c; = [dcsi_dx  dcsi_dy dcsi_dz; =  [ (1,1)  (1,2)  (1,3);
              d e f;    deta_dx  deta_dy deta_dz;      (2,1)  (2,2)  (2,3);
              g h i]    dzet_dx  dzet_dy dzet_dz;]     (3,1)  (3,2)  (3,3) ]

        Fc = [a b c; d e f; g h i].[Fxc; Fyc; Fzc] = [a*Fxc + b*Fyc + c*Fzc; d*Fxc + e*Fyc + f*Fzc; g*Fxc + h*Fyc + i*Fzc]
        Fc = [dcsi_dx*Fxc + dcsi_dy*Fyc + dcsi_dz*Fzc; 
              deta_dx*Fxc + deta_dy*Fyc + deta_dz*Fzc; 
              dzet_dx*Fxc + dzet_dy*Fyc + dzet_dz*Fzc ] 

        Fcsic = dcsi_dx*Fxc + dcsi_dy*Fyc + dcsi_dz*Fzc
        Fetac = deta_dx*Fxc + deta_dy*Fyc + deta_dz*Fzc
        Fzetc = dzet_dx*Fxc + dzet_dy*Fyc + dzet_dz*Fzc

        Fcsic = e->Fc[0]
        Fetac = e->Fc[1]
        Fzetc = e->Fc[2]
      */

      e->computational->Fcsp[0][s_index-1] = {e->Ji[1][1][s_index]*(q2) + e->Ji[1][2][s_index]*(q3) + e->Ji[1][3][s_index]*(q4),
                                              e->Ji[1][1][s_index]*(q2*q2/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[1][2][s_index]*(q3*q2/q1) + e->Ji[1][3][s_index]*(q4*q2/q1),
                                              e->Ji[1][1][s_index]*(q2*q3/q1) + e->Ji[1][2][s_index]*(q3*q3/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[1][3][s_index]*(q4*q3/q1),
                                              e->Ji[1][1][s_index]*(q2*q4/q1) + e->Ji[1][2][s_index]*(q3*q4/q1) + e->Ji[1][3][s_index]*(q4*q4/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)),
                                              e->Ji[1][1][s_index]*((q2/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[1][2][s_index]*((q3/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[1][3][s_index]*((q4/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1))};

      e->computational->Fcsp[1][s_index-1] = {e->Ji[2][1][s_index]*(q2) + e->Ji[2][2][s_index]*(q3) + e->Ji[2][3][s_index]*(q4),
                                              e->Ji[2][1][s_index]*(q2*q2/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[2][2][s_index]*(q3*q2/q1) + e->Ji[2][3][s_index]*(q4*q2/q1),
                                              e->Ji[2][1][s_index]*(q2*q3/q1) + e->Ji[2][2][s_index]*(q3*q3/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[2][3][s_index]*(q4*q3/q1),
                                              e->Ji[2][1][s_index]*(q2*q4/q1) + e->Ji[2][2][s_index]*(q3*q4/q1) + e->Ji[2][3][s_index]*(q4*q4/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)),
                                              e->Ji[2][1][s_index]*((q2/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[2][2][s_index]*((q3/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[2][3][s_index]*((q4/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1))};

      e->computational->Fcsp[2][s_index-1] = {e->Ji[3][1][s_index]* (q2)+e->Ji[3][2][s_index]* (q3)+e->Ji[3][3][s_index]* (q4),
                                              e->Ji[3][1][s_index]* (q2 * q2 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1)) + e->Ji[3][2][s_index]* (q3 * q2 / q1) + e->Ji[3][3][s_index]* (q4 * q2 / q1),
                                              e->Ji[3][1][s_index]* (q2 * q3 / q1) + e->Ji[3][2][s_index]* (q3 * q3 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1)) + e->Ji[3][3][s_index]* (q4 * q3 / q1),
                                              e->Ji[3][1][s_index]* (q2 * q4 / q1) + e->Ji[3][2][s_index]* (q3 * q4 / q1) + e->Ji[3][3][s_index]* (q4 * q4 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1)),
                                              e->Ji[3][1][s_index]* ((q2 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1)) + e->Ji[3][2][s_index]* ((q3 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1)) + e->Ji[3][3][s_index]*((q4/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1))};
      

      // e->Fcsp[0][s_index-1] = {q2,
      //                          q2*q2/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
      //                          q2*q3/q1,
      //                          q2*q4/q1,
      //                          (q2/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};

      // e->Fcsp[1][s_index-1] = {q3,
      //                          q3*q2/q1,
      //                          q3*q3/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
      //                          q3*q4/q1,
      //                          (q3/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};

      
      // e->Fcsp[2][s_index-1] = {q4,
      //                          q4*q2/q1,
      //                          q4*q3/q1,
      //                          q4*q4/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
      //                          (q4/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};
    }
  }
}

// Navier-Stokes
template <>
void SD<NavierStokes>::calculate_fluxes_sp (std::shared_ptr<Element>& e)
{
  double q1, q2, q3, q4, q5;
  double dq1_dx, dq2_dx, dq3_dx, dq4_dx, dq5_dx;
  double dq1_dy, dq2_dy, dq3_dy, dq4_dy, dq5_dy;
  double dq1_dz, dq2_dz, dq3_dz, dq4_dz, dq5_dz;
  double dT_dx, dT_dy, dT_dz;
  double du_dx, du_dy, du_dz;
  double dv_dx, dv_dy, dv_dz;
  double dw_dx, dw_dy, dw_dz;
  double txx, txy, txz;
  double tyx, tyy, tyz;
  double tzx, tzy, tzz;

  unsigned int s_index;
  
  s_index = 0;
  // for each solution point, calculate the flux vectors
  for (auto& node : this->snodes)
  {
    s_index++;
    // these Q must be the physical solution, so I divide the computational solution Q by |J|
    q1 = e->physical->Qsp[s_index - 1][0] / e->J;
    q2 = e->physical->Qsp[s_index - 1][1] / e->J;
    q3 = e->physical->Qsp[s_index - 1][2] / e->J;
    q4 = e->physical->Qsp[s_index - 1][3] / e->J;

    dq1_dx = e->physical->dQsp[0][s_index - 1][0];
    dq2_dx = e->physical->dQsp[0][s_index - 1][1];
    dq3_dx = e->physical->dQsp[0][s_index - 1][2];
    dq4_dx = e->physical->dQsp[0][s_index - 1][3];

    dq1_dy = e->physical->dQsp[1][s_index - 1][0];
    dq2_dy = e->physical->dQsp[1][s_index - 1][1];
    dq3_dy = e->physical->dQsp[1][s_index - 1][2];
    dq4_dy = e->physical->dQsp[1][s_index-1][3];

    if (this->dimension+2 == 4)
    {
      du_dx = (dq2_dx/q1 - q2*dq1_dx/(q1*q1));
      du_dy = (dq2_dy/q1 - q2*dq1_dy/(q1*q1));
    
      dv_dx = (dq3_dx/q1 - q2*dq1_dx/(q1*q1));
      dv_dy = (dq3_dy/q1 - q2*dq1_dx/(q1*q1));

      txx = 2.0*this->MU*du_dx + this->MUv*(du_dx + dv_dy);
      tyy = 2.0*this->MU*dv_dy + this->MUv*(du_dx + dv_dy);
      txy = this->MU*(dv_dx + du_dy);
      tyx = txy;

      /*
        Ji = [a b; = [dcsi_dx  dcsi_dy; =  [ (1,1)   (1,2);
              c d]    deta_dx  deta_dy]      (2,1)   (2,2)]

        Fc = [a b; c d].[Fxc; Fyc] = [a*Fxc + b*Fyc; c*Fxc + d*Fyc]
        Fc = [dcsi_dx*Fxc + dcsi_dy*Fyc; 
              deta_dx*Fxc + deta_dy*Fyc]

        Fcsic = dcsi_dx*Fxc + dcsi_dy*Fyc
        Fetac = deta_dx*Fxc + deta_dy*Fyc

        Fcsic = e->Fc[0]
        Fetac = e->Fc[1]
      */

      // Convective Flux
      e->computational->Fcsp[0][s_index-1] = {e->Ji[1][1][s_index]*q2,
                                              e->Ji[1][1][s_index]*(q2*q2/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1)),
                                              e->Ji[1][1][s_index]*(q2*q3/q1),
                                              e->Ji[1][1][s_index]*((q2/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1))};
      
      e->computational->Fcsp[0][s_index-1] = {q2,
                               q2*q2/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1),
                               q2*q3/q1,
                               (q2/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1)};
      
      e->computational->Fcsp[1][s_index-1] = {q3,
                               q3*q2/q1,
                               q3*q3/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1),
                               (q3/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1)};
     
      // Diffusive Flux
      e->computational->Fdsp[0][s_index-1] = {0.0,
                               txx,
                               tyx,
                               txx*q2/q1 + tyx*q3/q1 + this->KAPPA*dT_dx};
      
      e->computational->Fdsp[1][s_index-1] = {0.0,
                               txy,
                               tyy,
                               txy*q2/q1 + tyy*q3/q1 + this->KAPPA*dT_dy};
    }
    else if (this->dimension+2 == 5)
    {
      q5  = e->physical->Qsp[s_index-1][4];
      
      dq5_dx = e->physical->dQsp[0][s_index-1][4];
      
      dq5_dy = e->physical->dQsp[1][s_index - 1][4];

      dq1_dz = e->physical->dQsp[2][s_index - 1][0];
      dq2_dz = e->physical->dQsp[2][s_index - 1][1];
      dq3_dz = e->physical->dQsp[2][s_index - 1][2];
      dq4_dz = e->physical->dQsp[2][s_index - 1][3];
      dq5_dz = e->physical->dQsp[2][s_index-1][4];
    
      dw_dx = (dq4_dx/q1 - q2*dq1_dx/(q1*q1));
      dw_dy = (dq4_dy/q1 - q2*dq1_dy/(q1*q1));

      du_dz = (dq2_dz/q1 - q2*dq1_dz/(q1*q1));
      dv_dz = (dq3_dz/q1 - q2*dq1_dz/(q1*q1));
      dw_dz = (dq4_dz/q1 - q2*dq1_dz/(q1*q1));
      
      
      txx = 2.0*this->MU*du_dx + this->MUv*(du_dx + dv_dy + dw_dz);
      tyy = 2.0*this->MU*dv_dy + this->MUv*(du_dx + dv_dy + dw_dz);
      tzz = 2.0*this->MU*dw_dz + this->MUv*(du_dx + dv_dy + dw_dz);
      
      txy = this->MU*(dv_dx + du_dy);;
      tyx = txy;
      
      txz = this->MU*(dw_dx + du_dz);;
      tzx = txy;
      
      tzy = this->MU*(dv_dz + dw_dy);;
      tyz = tzy;

      // Convective Flux
      e->computational->Fcsp[0][s_index - 1] = { q2,
                               q2 * q2 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1),
                               q2 * q3 / q1,
                               q2 * q4 / q1,
                               (q2 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1) };

      e->computational->Fcsp[1][s_index - 1] = { q3,
                               q3 * q2 / q1,
                               q3 * q3 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1),
                               q3 * q4 / q1,
                               (q3 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1) };

      e->computational->Fcsp[2][s_index - 1] = { q4,
                               q4 * q2 / q1,
                               q4 * q3 / q1,
                               q4 * q4 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1),
                               (q4 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1) };

      // Diffussive Flux
      e->computational->Fdsp[0][s_index - 1] = { 0.0,
                               txx,
                               tyx,
                               tzx,
                               txx * q2 / q1 + tyx * q3 / q1 + tzx * q4 / q1 + this->KAPPA * dT_dx };

      e->computational->Fdsp[1][s_index - 1] = { 0.0,
                               txy,
                               tyy,
                               tzy,
                               txy * q2 / q1 + tyy * q3 / q1 + tzz * q4 / q1 + this->KAPPA * dT_dy };

      e->computational->Fdsp[2][s_index-1] = {0.0,
                               txz,
                               tyz,
                               tzz,
                               txz*q2/q1 + tyz*q3/q1 + tzz*q4/q1 + this->KAPPA*dT_dz};
    }
  }
}

// 3.2) Calculate Fluxes at Internal FPs
// Euler
template <>
void SD<Euler>::calculate_fluxes_fp (std::shared_ptr<Element>& e)
{
  double q1, q2, q3, q4, q5;
  unsigned int s_index;
  
  s_index = 0;
  // for each solution point, calculate the flux vector
  for (auto& node : this->snodes)
  {
    s_index++;

    q1 = e->physical->Qsp[s_index-1][0];
    q2 = e->physical->Qsp[s_index-1][1];
    q3 = e->physical->Qsp[s_index-1][2];
    q4 = e->physical->Qsp[s_index-1][3];

    if (this->dimension+2 == 4)
    {
      e->computational->Fcsp[0][s_index-1] = {q2, 
                               q2*q2/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1),
                               q2*q3/q1,
                               (q2/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1)};
      
      e->computational->Fcsp[1][s_index-1] = {q3,
                               q3*q2/q1,
                               q3*q3/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1),
                               (q3/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1)};
    }
    else if (this->dimension+2 == 5)
    {
      q5 = e->physical->Qsp[s_index-1][4];

      e->computational->Fcsp[0][s_index-1] = {q2,
                               q2*q2/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
                               q2*q3/q1,
                               q2*q4/q1,
                               (q2/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};

      e->computational->Fcsp[1][s_index-1] = {q3,
                               q3*q2/q1,
                               q3*q3/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
                               q3*q4/q1,
                               (q3/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};

      
      e->computational->Fcsp[2][s_index-1] = {q4,
                               q4*q2/q1,
                               q4*q3/q1,
                               q4*q4/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
                               (q4/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};
    }
  }
}

// Navier-Stokes
template <>
void SD<NavierStokes>::calculate_fluxes_fp (std::shared_ptr<Element>& e)
{
  double q1, q2, q3, q4, q5;
  double dq1_dx, dq2_dx, dq3_dx, dq4_dx, dq5_dx;
  double dq1_dy, dq2_dy, dq3_dy, dq4_dy, dq5_dy;
  double dq1_dz, dq2_dz, dq3_dz, dq4_dz, dq5_dz;
  double dT_dx, dT_dy, dT_dz;
  double du_dx, du_dy, du_dz;
  double dv_dx, dv_dy, dv_dz;
  double dw_dx, dw_dy, dw_dz;
  double txx, txy, txz;
  double tyx, tyy, tyz;
  double tzx, tzy, tzz;

  unsigned int f_index;
  
  f_index = 0;
  // for each flux point, calculate the flux vectors
  for (auto& node : this->fnodes[0])
  {
    f_index++;

    q1  = e->physical->Qfp[f_index-1][0];
    q2  = e->physical->Qfp[f_index-1][1];
    q3  = e->physical->Qfp[f_index-1][2];
    q4  = e->physical->Qfp[f_index-1][3];

    dq1_dx = e->physical->dQfp[0][f_index-1][0];
    dq2_dx = e->physical->dQfp[0][f_index-1][1];
    dq3_dx = e->physical->dQfp[0][f_index-1][2];
    dq4_dx = e->physical->dQfp[0][f_index-1][3];

    dq1_dy = e->physical->dQfp[1][f_index-1][0];
    dq2_dy = e->physical->dQfp[1][f_index-1][1];
    dq3_dy = e->physical->dQfp[1][f_index-1][2];
    dq4_dy = e->physical->dQfp[1][f_index-1][3];

    if (this->dimension+2 == 4)
    {
      du_dx = (dq2_dx/q1 - q2*dq1_dx/(q1*q1));
      du_dy = (dq2_dy/q1 - q2*dq1_dy/(q1*q1));
    
      dv_dx = (dq3_dx/q1 - q2*dq1_dx/(q1*q1));
      dv_dy = (dq3_dy/q1 - q2*dq1_dx/(q1*q1));

      txx = 2.0*this->MU*du_dx + this->MUv*(du_dx + dv_dy);
      tyy = 2.0*this->MU*dv_dy + this->MUv*(du_dx + dv_dy);
      txy = this->MU*(dv_dx + du_dy);
      tyx = txy;

      // Convective Flux
      e->computational->Fcfp[0][f_index-1] = {q2,
                               q2*q2/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1),
                               q2*q3/q1,
                               (q2/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1)};
      
      e->computational->Fcfp[1][f_index-1] = {q3,
                               q3*q2/q1,
                               q3*q3/q1 + (this->GAMMA-1.0)*(q4 - 0.5*(q2*q2 + q3*q3)/q1),
                               (q3/q1)*(this->GAMMA*q4 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3)/q1)};
     
      // Diffusive Flux
      e->computational->Fdfp[0][f_index-1] = {0.0,
                               txx,
                               tyx,
                               txx*q2/q1 + tyx*q3/q1 + this->KAPPA*dT_dx};
      
      e->computational->Fdfp[1][f_index-1] = {0.0,
                               txy,
                               tyy,
                               txy*q2/q1 + tyy*q3/q1 + this->KAPPA*dT_dy};
    }
    else if (this->dimension+2 == 5)
    {
      q5  = e->physical->Qfp[f_index-1][4];
      
      dq5_dx = e->physical->dQfp[0][f_index-1][4];
      
      dq5_dy = e->physical->dQfp[1][f_index-1][4];
      
      dq1_dz = e->physical->dQfp[2][f_index-1][0];
      dq2_dz = e->physical->dQfp[2][f_index-1][1];
      dq3_dz = e->physical->dQfp[2][f_index-1][2];
      dq4_dz = e->physical->dQfp[2][f_index-1][3];
      dq5_dz = e->physical->dQfp[2][f_index-1][4];
    
      dw_dx = (dq4_dx/q1 - q2*dq1_dx/(q1*q1));
      dw_dy = (dq4_dy/q1 - q2*dq1_dy/(q1*q1));

      du_dz = (dq2_dz/q1 - q2*dq1_dz/(q1*q1));
      dv_dz = (dq3_dz/q1 - q2*dq1_dz/(q1*q1));
      dw_dz = (dq4_dz/q1 - q2*dq1_dz/(q1*q1));
      
      
      txx = 2.0*this->MU*du_dx + this->MUv*(du_dx + dv_dy + dw_dz);
      tyy = 2.0*this->MU*dv_dy + this->MUv*(du_dx + dv_dy + dw_dz);
      tzz = 2.0*this->MU*dw_dz + this->MUv*(du_dx + dv_dy + dw_dz);
      
      txy = this->MU*(dv_dx + du_dy);;
      tyx = txy;
      
      txz = this->MU*(dw_dx + du_dz);;
      tzx = txy;
      
      tzy = this->MU*(dv_dz + dw_dy);;
      tyz = tzy;

      // Convective Flux
      e->computational->Fcfp[0][f_index-1] = {q2,
                               q2*q2/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
                               q2*q3/q1,
                               q2*q4/q1,
                               (q2/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};
      
      e->computational->Fcfp[1][f_index-1] = {q3,
                               q3*q2/q1,
                               q3*q3/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
                               q3*q4/q1,
                               (q3/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};
      
      e->computational->Fcfp[2][f_index-1] = {q4,
                               q4*q2/q1,
                               q4*q3/q1,
                               q4*q4/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
                               (q4/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};
     
      // Diffussive Flux
      e->computational->Fdfp[0][f_index-1] = {0.0,
                               txx,
                               tyx,
                               tzx,
                               txx*q2/q1 + tyx*q3/q1 + tzx*q4/q1 + this->KAPPA*dT_dx};

      e->computational->Fdfp[1][f_index-1] = {0.0,
                               txy,
                               tyy,
                               tzy,
                               txy*q2/q1 + tyy*q3/q1 + tzz*q4/q1 + this->KAPPA*dT_dy};

      e->computational->Fdfp[2][f_index-1] = {0.0,
                               txz,
                               tyz,
                               tzz,
                               txz*q2/q1 + tyz*q3/q1 + tzz*q4/q1 + this->KAPPA*dT_dz};
    }
  }
}


// 4) RIEMANN SOLVER
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

// 5) INTERPOLATE FROM FLUX POINTS TO SOLUTION POINTS
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

      e->computational->dFcsp[index-1][s_index-1] = 0.0;
  
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
        e->computational->dFcsp[index-1][s_index-1] += ((dLcsi*dLeta)*e->computational->Fcfp[index-1][f_index-1]);
        e->physical->dFcsp[index-1][s_index-1] += ((dLcsi*dLeta)*e->physical->Fcfp[index-1][f_index-1]);
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
      e->computational->dFcsp[index-1][s_index-1] = 0.0;
      e->computational->dFdsp[index-1][s_index-1] = 0.0;
  
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
        e->computational->dFcsp[index-1][s_index-1] += ((dLcsi*dLeta)*e->computational->Fcfp[index-1][f_index-1]);
        e->physical->dFdsp[index-1][s_index-1] += ((dLcsi*dLeta)*e->physical->Fdfp[index-1][f_index-1]);
      }
    }
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GLL>::delete_nodes();
}

// 6) RESIDUE 
// Euler
template <>
void SD<Euler> ::residue(std::shared_ptr<Element>& e)
{
  unsigned int index, s_index;

  s_index = 0;
  for (auto& node : this->snodes)
  {
    s_index++;

    e->computational->res[s_index-1] = 0.0;
    index = 0;
    for (auto& vec_lines : this->fnodes)
    {
      index++;
      e->computational->res[s_index-1] += (-e->computational->dFcsp[index-1][s_index-1]);
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

    e->computational->res[s_index-1] = 0.0;
    index = 0;
    for (auto& vec_lines : this->fnodes)
    {
      index++;
      e->computational->res[s_index-1] += (-e->computational->dFcsp[index-1][s_index-1] + (this->M/this->Re)*e->computational->dFdsp[index-1][s_index-1]);
    }
  }
}

// 7) SOLVE
template <typename Equation>
void SD<Equation>::solve (std::shared_ptr<Mesh>& mesh)
{
  
  // Step 0)
  for (auto& g : mesh->ghosts)
  {
    this->boundary_condition(g, mesh->elems);
  }
  
  // Step 1)
  for (auto& e : mesh->elems)
  {
    this->interpolate_sp2fp(e);
    //this->calculate_fluxes_sp(e);
    this->calculate_fluxes_fp(e);
  }
  
  // Step 2)
  for (auto& e : mesh->elems)
  {
    this->calculate_interface_fluxes(e, mesh->elems, mesh->ghosts);
    this->riemann_solver(e);
    this->interpolate_fp2sp(e);
    this->residue(e);
  }
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

