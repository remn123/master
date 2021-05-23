#pragma once

#include <algorithm>
#include <any>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <sstream>
#include <thread>
#include <vector>
#include <unordered_map>
#include <unordered_set>


#include <Dummies.h>
#include <Helpers.h>
#include <Mesh.h>
#include <Node.h>
#include <Poly.h>
#include <SD.h>

// explicit instances
template class SD<Euler>;
template class SD<NavierStokes>;

// Constructor for Euler Equations
template <>
SD<Euler>::SD(int order, int dimension)
{
  this->order = order;
  this->dimension = dimension;
  this->GAMMA = 1.4;
  this->MU = 0.0;
  this->MUv = 0.0;
  std::size_t snodes_size = order * order;
  std::size_t fnodes_size = (order + 1) * order;
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
  this->GAMMA = 1.4;
  this->M = Mach;
  this->Re = Reynolds;
  this->MU = mu;
  this->MUv = (-2.0 / 3.0) * mu;
  std::size_t snodes_size = order * order;
  std::size_t fnodes_size = (order + 1) * order;
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
  std::size_t s_index = 0;
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
  Helpers<GLL>::set_nodes(order + 1); // Flux must be one order higher
  std::vector<double> nodes_gll = Helpers<GLL>::get_nodes();

  Helpers<GLL>::delete_nodes();

  std::vector<Node> vecx;
  std::vector<Node> vecy;

  std::size_t f_index = 0;

  if (dimension == 2)
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
    f_index = 0;
    for (auto n1 : nodes_gl) // x
    {
      for (auto n2 : nodes_gll) // y
      {
        this->fnodes[1][f_index] = Node{n1, n2, 0.0};
        f_index++;
      }
    }
  }
  else if (dimension == 3)
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

template <typename Equation>
void SD<Equation>::create_weights(void)
{
  // this->weights.clear();
  // this->weights.resize(num_nodes);
}


// Auxiliar function
template <typename Equation>
void SD<Equation>::_init_dvec(std::vector<DVector> &vec, size_t num_nodes)
{
  std::vector<double> init_vec(this->dimension + 2, 0.0);

  vec.clear();
  vec.resize(num_nodes);
  for (auto &dvec : vec)
  {
    dvec = init_vec;
  }
}

template <typename Equation>
void SD<Equation>::_init_dvec(std::vector<std::vector<DVector>> &vec, size_t num_nodes)
{
  std::vector<double> init_vec(this->dimension + 2, 0.0);

  vec.clear();
  vec.resize(this->dimension);
  for (auto &dirvec : vec) // x, y, ...
  {
    dirvec.clear();
    dirvec.resize(num_nodes); // number of nodes
    for (auto &dvec : dirvec)
    {
      dvec = init_vec;
    }
  }
}

// 0.2) Initialize property template
// INITIALIZE GHOST PROPERTIES
template <>
void SD<Euler>::initialize_ghost_properties(Ghost &g)
{
  // PHYSICAL
  // Conservative Properties
  this->_init_dvec(g.physical->Qsp, this->snodes.size());
  this->_init_dvec(g.physical->Qfp, this->fnodes[0].size());

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
void SD<NavierStokes>::initialize_ghost_properties(Ghost &g)
{
  // PHYSICAL
  // Conservative Properties
  this->_init_dvec(g.physical->Qsp, this->snodes.size());
  this->_init_dvec(g.physical->Qfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(g.physical->Fcsp, this->snodes.size());
  this->_init_dvec(g.physical->Fcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(g.physical->Fdsp, this->snodes.size());
  this->_init_dvec(g.physical->Fdfp, this->fnodes[0].size());

  // Gradients
  // Conservative Properties
  this->_init_dvec(g.physical->dQsp, this->snodes.size());
  this->_init_dvec(g.physical->dQfp, this->fnodes[0].size());

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
void SD<Euler>::initialize_element_properties(std::shared_ptr<Element> &e,
                                              const std::vector<Vertice> &enodes,
                                              std::vector<double> (*field)(const Node &))
{
  // Metrics
  //e->allocate_jacobian(this->order);
  //e->calculate_jacobian(this->snodes, this->fnodes, enodes);

  e->physical->Qsp.clear();
  e->physical->Qsp.resize(this->snodes.size());
  e->computational->Qsp.clear();
  e->computational->Qsp.resize(this->snodes.size());

  // PHYSICAL
  Node n;
  int index = 0;
  for (auto& sp : this->snodes)
  {
    n = e->transform(sp, enodes);
    auto vec = field(n);
    e->physical->Qsp[index] = vec;
    e->computational->Qsp[index] = e->J[0][index]*e->physical->Qsp[index];
    index++;
  }

  // Conservative Properties
  //this->_init_dvec(e->physical->Qsp, this->snodes.size());
  this->_init_dvec(e->physical->Qfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(e->physical->Fcsp, this->snodes.size());
  this->_init_dvec(e->physical->Fcfp, this->fnodes[0].size());

  // Gradients
  // Convective Fluxes
  this->_init_dvec(e->physical->dFcsp, this->snodes.size());
  this->_init_dvec(e->physical->dFcfp, this->fnodes[0].size());

  std::size_t global_fn = 0;
  std::size_t local_face = 0;
  std::size_t dir = 0, indice = 0, swap = 0;
  std::vector<fNode> fn;

  
  for (auto &ed : e->edges)
  {
    this->_init_dvec(ed.physical->Qfp, this->order);
    this->_init_dvec(ed.physical->Fcfp, this->order);
    this->_init_dvec(ed.physical->dFcfp, this->order);

    this->_init_dvec(ed.computational->Qfp, this->order);
    this->_init_dvec(ed.computational->Fcfp, this->order);
    this->_init_dvec(ed.computational->dFcfp, this->order);

    fn.reserve(this->order);
    if (local_face == 0 || local_face == 2)
      dir = 1;
    else
      dir = 0;

    if (local_face == 0 || local_face == 3)
      indice = 0;
    else
      indice = 1;
    
    if (local_face == 0 || local_face == 1)
      swap = 0;
    else
      swap = 1;


    for (size_t idx = 0; idx < this->order; idx++)
    {
      global_fn = this->order * indice + (this->order + 1) * (idx + swap*(this->order-2*idx-1));
      //std::cout << " e[" << e->id << "]->edge[" << ed.id << "]: f.local = " << global_fn << " - f.id = " << idx << "\n";
      // Project node from computational to physical space
      n = e->transform(this->fnodes[dir][global_fn], enodes); // need copy constructor
      fn.push_back(fNode{idx, global_fn, n});
    }

    ed.fnodes = fn;

    local_face++;
    fn.clear();
  }

  // COMPUTATIONAL
  // Conservative Properties
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
void SD<NavierStokes>::initialize_element_properties(std::shared_ptr<Element> &e, 
                                                     const std::vector<Vertice> &enodes,
                                                     std::vector<double> (*field)(const Node &))
{
  // Metrics
  e->allocate_jacobian(this->order);
  e->calculate_jacobian(this->snodes, this->fnodes, enodes);

  e->physical->Qsp.clear();
  e->physical->Qsp.resize(this->snodes.size());

  // PHYSICAL
  Node n;
  int index = 0;
  for (auto sp : this->snodes)
  {
    n = e->transform(sp, enodes);
    auto vec = field(n);
    e->physical->Qsp[index] = vec;

    index++;
  }

  // PHYSICAL
  // Conservative Properties
  //this->_init_dvec(e->physical->Qsp, this->snodes.size());
  this->_init_dvec(e->physical->Qfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(e->physical->Fcsp, this->snodes.size());
  this->_init_dvec(e->physical->Fcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(e->physical->Fdsp, this->snodes.size());
  this->_init_dvec(e->physical->Fdfp, this->fnodes[0].size());

  // Gradients
  // Conservative Properties
  this->_init_dvec(e->physical->dQsp, this->snodes.size());
  this->_init_dvec(e->physical->dQfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(e->physical->dFcsp, this->snodes.size());
  this->_init_dvec(e->physical->dFcfp, this->fnodes[0].size());

  // Diffusive Fluxes
  this->_init_dvec(e->physical->dFdsp, this->snodes.size());
  this->_init_dvec(e->physical->dFdfp, this->fnodes[0].size());

  std::size_t global_fn = 0;
  std::size_t local_face = 0;
  std::size_t dir = 0, indice = 0, swap =0;
  std::vector<fNode> fn;

  for (auto &ed : e->edges)
  {
    this->_init_dvec(ed.physical->Qfp, this->order);
    this->_init_dvec(ed.physical->dQfp, this->order);
    this->_init_dvec(ed.physical->Fcfp, this->order);
    this->_init_dvec(ed.physical->dFcfp, this->order);
    this->_init_dvec(ed.physical->Fdfp, this->order);
    this->_init_dvec(ed.physical->dFdfp, this->order);

    this->_init_dvec(ed.computational->Qfp, this->order);
    this->_init_dvec(ed.computational->dQfp, this->order);
    this->_init_dvec(ed.computational->Fcfp, this->order);
    this->_init_dvec(ed.computational->dFcfp, this->order);
    this->_init_dvec(ed.computational->Fdfp, this->order);
    this->_init_dvec(ed.computational->dFdfp, this->order);

    fn.reserve(this->order);
    if (local_face == 0 || local_face == 2)
      dir = 1;
    else
      dir = 0;

    if (local_face == 0 || local_face == 3)
      indice = 0;
    else
      indice = 1;

    if (local_face == 0 || local_face == 1)
      swap = 0;
    else
      swap = 1;


    for (size_t idx = 0; idx < this->order; idx++)
    {
      global_fn = this->order * indice + (this->order + 1) * (idx + swap*(this->order-2*idx-1));
      // Project node from computational to physical space
      n = e->transform(this->fnodes[dir][global_fn], enodes); // need copy constructor
      fn.push_back(fNode{idx, global_fn, n});
    }

    ed.fnodes = fn;

    local_face++;
    fn.clear();
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

// UPDATE EDGES
template <typename Equation>
void SD<Equation>::update_edges(std::shared_ptr<Element> &e, 
                                std::vector<std::shared_ptr<Element>> &elems, 
                                std::vector<Ghost> &ghosts)
{
  long neighbor = 0, ghost = 0;
  long local_ed1 = 0, local_ed2 = 0;
  std::size_t global_fn = 0;
  std::size_t dir = 0;
  std::vector<fNode> fn;

  for (auto &ed1 : e->edges)
  {
    neighbor = ed1.right;
    ghost = ed1.ghost;
    if (neighbor > -1) // cell
    {
      // for each fnode in edge ed
      for (auto &fn1 : ed1.fnodes)
      {
        if (fn1.right != -1)
          break; // already updated

        local_ed2 = 0;
        // search through all neighbor's edges
        for (auto &ed2 : elems[neighbor]->edges)
        {
          //std::cout << "Element " << e->id << ", Edge " << ed1.id << ", Neighbor " << neighbor << ", Edge " << ed2.id << std::endl;
          //std::cout << "Edge " << ed1.id << " -> nodes: {" << ed1.nodes[0] << ", " << ed1.nodes[1] << "}" << std::endl;
          //std::cout << "Edge " << ed2.id << " -> nodes: {" << ed2.nodes[0] << ", " << ed2.nodes[1] << "}" << std::endl;
          // and look for the flux node which has the "same" physical coordinates
          for (auto &fn2 : ed2.fnodes)
          {
            //std::cout << "fn1 " << fn1.id << ", local " << fn1.local << ", right " << fn1.right << ", x " << fn1.coords[0] << ", y " << fn1.coords[1] << std::endl;
            //std::cout << "fn2 " << fn2.id << ", local " << fn2.local << ", right " << fn2.right << ", x " << fn2.coords[0] << ", y " << fn2.coords[1] << std::endl;

            if (std::abs(fn1.coords[0] - fn2.coords[0]) <= 1E-7 &&
                std::abs(fn1.coords[1] - fn2.coords[1]) <= 1E-7 &&
                std::abs(fn1.coords[2] - fn2.coords[2]) <= 1E-7)
            {
              //std::cout << "Found it!" << std::endl;
              fn1.right = fn2.local;
              fn2.right = fn1.local;
              break;
            }
          }

          if (fn1.right != -1)
          {
            ed1.lr_edge = local_ed2;
            ed2.lr_edge = local_ed1;
            break; // just updated
          }
          local_ed2++;
        }
        if (fn1.right == -1)
        {
          std::cout << "Not Found!" << std::endl;
          throw "Flux Node " +
              std::to_string(fn1.id) +
              " at Edge " +
              std::to_string(ed1.id) +
              " not found in Element " +
              std::to_string(neighbor) +
              "\n";
        }
      }
    }
    else // ghost cells
    {
      fn.reserve(this->order);
      // for each fnode in edge ed
      for (auto &fn1 : ed1.fnodes)
      {
        if (fn1.right != -1)
          break; // already updated

        // Appending fNode copy from edge1 at element 1 into ghost fnodes
        fn.push_back(fNode{fn1.id, fn1.local, fn1.local, fn1.coords});

        // Update communication
        fn1.right = fn1.local; // I mirror the edge order into the ghost edge

        if (fn1.right != -1)
        {
          ed1.lr_edge = -1;
          if (ghosts[ghost].lr_edge != local_ed1)
          {
            throw "Ghost communication mismatch. Ghost " +
                std::to_string(ghost) +
                " has local right edge " +
                std::to_string(ghosts[ghost].lr_edge) +
                " but " +
                std::to_string(local_ed1) +
                " was found! \n";
          }
        }
        else
        {
          throw "Flux Node " +
              std::to_string(fn1.id) +
              " at Edge " +
              std::to_string(ed1.id) +
              " not found in Ghost " +
              std::to_string(ghost) +
              "\n";
        }
      }
      ghosts[ghost].fnodes = fn;
      fn.clear();
    }
    local_ed1++;
  }
}

// 0)
template <typename Equation>
void SD<Equation>::setup(
  std::shared_ptr<Mesh> &mesh, 
  std::vector<double> (*field)(const Node &)
)
{
  /* 
    On setup method, all solution and flux nodes
    are allocated and calculated in create_nodes.

    Then initialize_properties will allocate and
    populate the mesh with an initial condition
  */

  //Ghost::Qbnds = bnd_map;


  this->create_nodes();

  //this->calculate_weights();

  for (auto &e : mesh->elems) {
    // Metrics
    e->allocate_jacobian(this->order);
    e->calculate_jacobian(this->snodes, this->fnodes, mesh->nodes);
  }

  std::cout << "Initializing Element Properties! \n";
  // std::vector<std::thread> threads;
  // for (auto &e : mesh->elems) {
  //   threads.push_back(std::thread(&SD<Equation>::initialize_element_properties, this, std::ref(e), std::ref(mesh->nodes), std::ref(field)));
  // }
 
  // for (auto &th : threads) {
  //   th.join();
  // }

  
  
  for (auto &e : mesh->elems)
    this->initialize_element_properties(e, mesh->nodes, field);
  std::cout << "Element Properties initialized! \n";

  for (auto &g : mesh->ghosts)
    this->initialize_ghost_properties(g);

  // Update edge communication
  for (auto &e : mesh->elems)
    this->update_edges(e, mesh->elems, mesh->ghosts);

  std::cout << "Setup completed\n";
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
template <>
void SD<Euler>::boundary_condition(
  Ghost &g, 
  const std::vector<std::shared_ptr<Element>> &elems,
  const std::vector<Vertice> &enodes
)
{ 
  // std::cout << "Ghost Cell ID = " << g.id << "; Type = " << g.type << "\n";
  // ghost face direction (0: csi/x and 1: eta/y)
  int dir = (g.lr_edge == 0 || g.lr_edge == 2) ? 1 : 0;
  auto gamma = this->GAMMA;
  switch (g.type)
  {
  //- Solid Inviscid Wall (Euler) [TYPE 0]:
  case PhysicalEnum::WALL:
    // Un - Uwall,n = (U-Uwall).n = 0
    // Qn
    // std::cout << "PhysicalEnum::WALL\n";
    // std::cin.get();
    for (auto &fn : g.fnodes)
    { 
      // Normal unit vector in Physical Space
      auto normal = elems[g.elm_id]->get_normal_vector(dir+1, fn.right, g.lr_edge);
      auto nx = normal[0];
      auto ny = normal[1];
      // auto signN = ((g.lr_edge==0||g.lr_edge==3)) ? -1 : 1;
      // auto ds_dx = this->Ji[dir+1][fn.local][2 * dir];
      // auto ds_dy = this->Ji[dir+1][fn.local][2 * dir + 1];
      auto J = elems[g.elm_id]->J[dir+1][fn.right];
      
      auto rho = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][0];
      auto u = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][1]/rho;
      auto v = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][2]/rho;
      auto E = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][3];
      
      g.computational->Qfp[dir][fn.local][0] = J*rho;
      // g.computational->Qfp[dir][fn.local][1] = J*rho*(-u*(nx + ny));
      // g.computational->Qfp[dir][fn.local][2] = J*rho*(v*(nx - ny));
      g.computational->Qfp[dir][fn.local][1] = J*rho*(u-2.0*(u*nx + v*ny)*nx);
      g.computational->Qfp[dir][fn.local][2] = J*rho*(v-2.0*(u*nx + v*ny)*ny);
      g.computational->Qfp[dir][fn.local][3] = J*E; // E
    }
    break;
  // Inlet BC or Far Field
  case PhysicalEnum::INLET:
    for (auto &fn : g.fnodes)
    {
      //std::cin.get();
      auto normal = elems[g.elm_id]->get_normal_vector(dir+1, fn.right, g.lr_edge);
      auto nx = normal[0];
      auto ny = normal[1];
      auto J = elems[g.elm_id]->J[dir+1][fn.right];

      // Ghost's Computational Properties
      auto rho = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][0];
      auto u = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][1] / rho;
      auto v = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][2] / rho;
      auto E = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][3];
      auto p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
      auto ek = rho*(0.5*(u*u+v*v));
      auto a = std::sqrt(gamma*p/rho);
     
      auto vec = std::any_cast<std::vector<double>(*)(const Node&)> (Ghost::analytical_solution)(fn);
      DVector Qbnd = DVector(vec);

      //auto un = std::abs((Qbnd[1]/Qbnd[0])*nx) + std::abs((Qbnd[2]/Qbnd[0])*ny);
      auto un = std::abs(u*nx + v*ny);
      if (un < a) // Subsonic Inlet
      {
        // auto c_star = 2.0*gamma*((gamma-1.0)/(gamma+1.0))*Ghost::Cv*Ghost::T;
        // auto Tt = Ghost::T*(1.0 - ((gamma-1.0)/(gamma+1.0))*(u*u+v*v)/std::pow(c_star, 2.0));
        // auto Pt = Ghost::p*std::pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(u*u+v*v)/std::pow(c_star, 2.0)), gamma/(gamma-1.0));

        // // Tt = Tt     / temp_ref;
        // // Pt = Pt     / ene_ref;

        // //! Compute squared critical speed of sound.
        // auto asq   = gamma * F_R * Tt;
        // auto acrsq = 2.0*asq / (gamma+1.0);

        // //! Compute the static pressure and other variables from stagnation properties.
        // auto auxcte1 = (gamma-1.0) / (gamma+1.0);
        // auto auxcte2 = gamma / (gamma-1.0);

        // auto aux1 = 1.0 - (auxcte1 * (u*u + v*v)/acrsq);
        // auto eint = F_Cv * Tt * aux1;
        // p    = Pt * (std::pow(aux1, auxcte2));
        // rho  = p / ( (gamma-1.0)*eint );
        u = 2.0*Ghost::Mach - u;
        v = -v;
        rho = 2.0-rho;
        E  = p/(gamma - 1.0) + 0.5 * rho * (u * u + v * v);
        g.computational->Qfp[dir][fn.local][0] = J*rho;
        g.computational->Qfp[dir][fn.local][1] = J*rho*u;
        g.computational->Qfp[dir][fn.local][2] = J*rho*v;
        g.computational->Qfp[dir][fn.local][3] = J*E;
        //std::cout << "Subsonic Inlet \n";
        // auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
        // auto uface = 2.0*(Qbnd[1]/Qbnd[0])-(Qint[1]/Qint[0]);
        // auto vface = 2.0*(Qbnd[2]/Qbnd[0])-(Qint[2]/Qint[0]);
        // auto uface = Qint[1]/Qint[0];
        // auto vface = Qint[2]/Qint[0];
        // ek = (0.5*(uface*uface+vface*vface));

        // g.computational->Qfp[dir][fn.local][0] = J*Qbnd[0];
        // g.computational->Qfp[dir][fn.local][1] = J*Qbnd[0]*uface;
        // g.computational->Qfp[dir][fn.local][2] = J*Qbnd[0]*vface;
        // g.computational->Qfp[dir][fn.local][3] = J*(p/(gamma-1.0) + Qbnd[0]*ek);
        // g.computational->Qfp[dir][fn.local][0] = J*Qbnd[0];
        // g.computational->Qfp[dir][fn.local][1] = J*Qbnd[0]*uface;
        // g.computational->Qfp[dir][fn.local][2] = J*Qbnd[0]*vface;
        // g.computational->Qfp[dir][fn.local][3] = J*Qbnd[3];
        

      }
      else // Supersonic Inlet
      {
        //  aisen = ( Pt / press )**(gamm1/F_gamma)

        // pe   = press
        // te   = Tt / aisen
        // rhoe = pe / (F_R*te)
        // a    = sqrt( F_gamma * F_R * te )

        // amach = sqrt( 2.0_rdp*(aisen - 1.0_rdp)/gamm1 )
        // vel   = amach*a
        // ue    = vel*udir
        // ve    = vel*vdir

        // ! Update the ghost array.
        // sdmesh%face_list(fID)%rfp_pointer(i)%p%u%d  = rho
        // sdmesh%face_list(fID)%rfp_pointer(i)%p%u%mx = rho * ue
        // sdmesh%face_list(fID)%rfp_pointer(i)%p%u%my = rho * ve
        // sdmesh%face_list(fID)%rfp_pointer(i)%p%u%e  = pe/gamm1 + rhoe*( 0.5_rdp*( ue*ue + ve*ve ) )
        //std::cout << "Supersonic Inlet \n";
        // Inplace
        std::transform(Qbnd.begin(), Qbnd.end(), Qbnd.begin(), [&J](auto &c) { return c * J; });

        //g.computational->Qfp[dir][fn.local] = Qbnd;
        auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
        // Ghost's Computational Properties
        g.computational->Qfp[dir][fn.local] = (2.0 * Qbnd) - Qint;
      }
    }
    break;
  // Outlet BC (Neumann BC)
  case PhysicalEnum::OUTLET:
    for (auto &fn : g.fnodes)
    {
      //std::cin.get();
      auto normal = elems[g.elm_id]->get_normal_vector(dir+1, fn.right, g.lr_edge);
      auto nx = normal[0];
      auto ny = normal[1];
      auto J = elems[g.elm_id]->J[dir+1][fn.right];

      // Ghost's Computational Properties
      auto rho = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][0];
      auto u = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][1] / rho;
      auto v = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][2] / rho;
      auto E = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][3];
      auto p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
      auto ek = rho*(0.5*(u*u+v*v));
      auto a = std::sqrt(gamma*p/rho);

      auto un = std::abs(u*nx + v*ny);
      if (un >= a) // Supersonic Outlet
      {
        //std::cout << "Supersonic Outlet \n";
        auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
        g.computational->Qfp[dir][fn.local] = Qint;
      }
      else // Subsonic Outlet
      {
        //std::cout << "Subsonic Outlet \n";
        auto pbc = 2.0*Ghost::p - p;
        //auto pbc = Ghost::p;

        g.computational->Qfp[dir][fn.local][0] = J*rho;
        g.computational->Qfp[dir][fn.local][1] = J*rho*u;
        g.computational->Qfp[dir][fn.local][2] = J*rho*v;
        g.computational->Qfp[dir][fn.local][3] = J*((pbc/(gamma-1.0)) + ek);
      }
    }
    break;
  case PhysicalEnum::RINGLEB_WALL:
    // Un - Uwall,n = (U-Uwall).n = 0
    // Qn
    //std::cout << "PhysicalEnum::RINGLEB_WALL\n";
    for (auto &fn : g.fnodes)
    { 
      auto J = elems[g.elm_id]->J[dir+1][fn.right];
      auto vec = std::any_cast<std::vector<double>(*)(const Node&)> (Ghost::analytical_solution)(fn);
      DVector Qbnd = DVector(vec);
      std::transform(Qbnd.begin(), Qbnd.end(), Qbnd.begin(), [&J](auto &c) { return c * J; });
      auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
      g.computational->Qfp[dir][fn.local] = (2.0 * Qbnd) - Qint;
    }
    break;
  // Supersonic Inlet BC or Far Field
  case PhysicalEnum::RINGLEB_INLET:
    // std::cout << "PhysicalEnum::SUPERSONIC_INLET\n";
    for (auto &fn : g.fnodes)
    {
      auto normal = elems[g.elm_id]->get_normal_vector(dir+1, fn.right, g.lr_edge);
      auto nx = normal[0];
      auto ny = normal[1];
      auto J = elems[g.elm_id]->J[dir+1][fn.right];

      // Ghost's Computational Properties
      auto rho = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][0];
      auto u = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][1] / rho;
      auto v = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][2] / rho;
      auto E = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][3];
      auto p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
      auto ek = rho*(0.5*(u*u+v*v));
      auto a = std::sqrt(gamma*p/rho);
      auto vec = std::any_cast<std::vector<double>(*)(const Node&)> (Ghost::analytical_solution)(fn);
      DVector Qbnd = DVector(vec);

      auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
      //auto un = std::abs((Qbnd[1]/Qbnd[0])*nx) + std::abs((Qbnd[2]/Qbnd[0])*ny);
      //auto un = std::sqrt(u*u + v*v);
      //auto un = std::abs(u*nx + v*ny);
      auto un = u*nx + v*ny;
      if (un < a) // Subsonic Inlet
      {
        //std::cout << "Ringleb Subsonic Inlet \n";
        rho = 2.0*Qbnd[0] - rho;
        auto uface = 2.0*(Qbnd[1]/Qbnd[0])-u;
        auto vface = 2.0*(Qbnd[2]/Qbnd[0])-v;
        // auto uface = Qint[1]/Qint[0];
        // auto vface = Qint[2]/Qint[0];
        
        ek = (0.5*(uface*uface+vface*vface));
        

        g.computational->Qfp[dir][fn.local][0] = J*rho;
        g.computational->Qfp[dir][fn.local][1] = J*rho*uface;
        g.computational->Qfp[dir][fn.local][2] = J*rho*vface;
        //g.computational->Qfp[dir][fn.local][3] = J*Qbnd[3];
        g.computational->Qfp[dir][fn.local][3] = J*(p/(gamma-1.0) + rho*ek);
      }
      else // Supersonic Inlet
      {
        // auto ureal = Qbnd[1]/Qbnd[0];
        // auto vreal = Qbnd[2]/Qbnd[0];
        // ek = 0.5*( ureal*ureal + vreal*vreal );

        // auto aux  = 1.0 - (((gamma-1.0)/(gamma+1.0)))*(ureal*ureal + vreal*vreal));
        // auto eint = ((gamma+1.0)/(2.0*gamma*(gamma-1.0))) * aux;
        // p = ((gamma+1.0)/(2.0*gamma))*(std::pow(aux,(gamma/(gamma-1.0))));
        // rho  = p / ( (gamma-1.0)*eint );

        // fpr%u%d  = urho
        // fpr%u%mx = urho * uface
        // fpr%u%my = urho * vface
        // fpr%u%e  = ue
        //std::cout << "Ringleb Supersonic Inlet \n";
        //auto J = elems[g.elm_id]->J[dir + 1][fn.right]; // J[0] for snodes, J[1] for csi-fnodes and J[2] for eta-fnodes
        // Inplace
        //std::transform(Qbnd.begin(), Qbnd.end(), Qbnd.begin(), [&J](auto &c) { return c * J; });

        std::transform(Qbnd.begin(), Qbnd.end(), Qbnd.begin(), [&J](auto &c) { return c * J; });

        //g.computational->Qfp[dir][fn.local] = Qbnd;
        //auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
        // Ghost's Computational Properties
        g.computational->Qfp[dir][fn.local] = (2.0 * Qbnd) - Qint;
        //auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
        // Ghost's Computational Properties
        //g.computational->Qfp[dir][fn.local] = (2.0 * Qbnd) - Qint;
      }

      // auto vec = std::any_cast<std::vector<double>(*)(const Node&)> (Ghost::analytical_solution)(fn);
      // DVector Qbnd = DVector(vec);
     
      // auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
      // auto J = elems[g.elm_id]->J[dir + 1][fn.right]; // J[0] for snodes, J[1] for csi-fnodes and J[2] for eta-fnodes
      // // Inplace
      // std::transform(Qbnd.begin(), Qbnd.end(), Qbnd.begin(), [&J](auto &c) { return c * J; });
      // // Ghost's Computational Properties
      // g.computational->Qfp[dir][fn.local] = (2.0 * Qbnd) - Qint;
    }
    break;
  // Supersonic Outlet BC (Neumann BC)
  case PhysicalEnum::RINGLEB_OUTLET:
    // std::cout << "PhysicalEnum::SUPERSONIC_OUTLET\n";
    for (auto &fn : g.fnodes)
    {
      auto normal = elems[g.elm_id]->get_normal_vector(dir+1, fn.right, g.lr_edge);
      auto nx = normal[0];
      auto ny = normal[1];
      auto J = elems[g.elm_id]->J[dir+1][fn.right];

      // Ghost's Computational Properties
      auto rho = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][0];
      auto u = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][1] / rho;
      auto v = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][2] / rho;
      auto E = (1.0/J)*elems[g.elm_id]->computational->Qfp[dir][fn.right][3];
      auto p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
      auto ek = rho*(0.5*(u*u+v*v));
      auto a = std::sqrt(gamma*p/rho);

      auto vec = std::any_cast<std::vector<double>(*)(const Node&)> (Ghost::analytical_solution)(fn);
      DVector Qbnd = DVector(vec);

      //auto un = std::abs(u*nx + v*ny);
      auto un = u*nx + v*ny;
      //auto un = std::sqrt(u*u + v*v);
      if (un >= a) // Supersonic Outlet
      {
        //std::cout << "Ringleb Supersonic Outlet \n";
        auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
        g.computational->Qfp[dir][fn.local] = Qint;
        // Ghost's Computational Properties
        // Inplace
        // std::transform(Qbnd.begin(), Qbnd.end(), Qbnd.begin(), [&J](auto &c) { return c * J; });
        // g.computational->Qfp[dir][fn.local] = (2.0 * Qbnd) - Qint;
      }
      else // Subsonic Outlet
      {
        //std::cout << "Ringleb Subsonic Outlet \n";
        //auto J = elems[g.elm_id]->J[dir + 1][fn.right]; // J[0] for snodes, J[1] for csi-fnodes and J[2] for eta-fnodes
        // Inplace
        //std::transform(Qbnd.begin(), Qbnd.end(), Qbnd.begin(), [&J](auto &c) { return c * J; });
        
        auto pbc = (gamma - 1.0) * (Qbnd[3] - 0.5 * (Qbnd[1] * Qbnd[1] + Qbnd[2] * Qbnd[2])/Qbnd[0]);
        //auto pbc = 2.0*((gamma - 1.0) * (Qbnd[3] - 0.5 * (Qbnd[1] * Qbnd[1] + Qbnd[2] * Qbnd[2])/Qbnd[0])) - p;
        
        g.computational->Qfp[dir][fn.local][0] = J*rho;
        g.computational->Qfp[dir][fn.local][1] = J*rho*u;
        g.computational->Qfp[dir][fn.local][2] = J*rho*v;
        g.computational->Qfp[dir][fn.local][3] = J*(((2.0*pbc - p)/(gamma-1.0)) + ek);
        //g.computational->Qfp[dir][fn.local][3] = J*((2.0*pbc - p)/this->GAMMA + ek);
      }
    }
    break;
  }
}

template <>
void SD<NavierStokes>::boundary_condition(
  Ghost &g, 
  const std::vector<std::shared_ptr<Element>> &elems,
  const std::vector<Vertice> &enodes
)
{
  // std::cout << "Ghost Cell ID = " << g.id << "\n";
  // ghost face direction (0: csi/x and 1: eta/y)
  int dir = (g.lr_edge == 0 || g.lr_edge == 2) ? 1 : 0;
  switch (g.type)
  {
  //- Solid Inviscid Wall (Euler) [TYPE 0]:
  case 0:
    // Un - Uwall,n = (U-Uwall).n = 0
    // Qn
    for (auto &fn : g.fnodes)
    {
      // Ghost's Computational Properties
      g.computational->Qfp[dir][fn.local][0] = elems[g.elm_id]->computational->Qfp[dir][fn.right][0];  // rho
      g.computational->Qfp[dir][fn.local][1] = -elems[g.elm_id]->computational->Qfp[dir][fn.right][1]; // rho*u
      g.computational->Qfp[dir][fn.local][2] = -elems[g.elm_id]->computational->Qfp[dir][fn.right][2]; // rho*v
      g.computational->Qfp[dir][fn.local][3] = elems[g.elm_id]->computational->Qfp[dir][fn.right][3];  // E
    }
    break;
  // Supersonic Inlet BC or Far Field (Dirichlet BC)
  case 1:
    for (auto &fn : g.fnodes)
    {
      auto type = (int)g.type;
      //auto n = elems[g.elm_id]->transform(fn, enodes);
      auto field = std::any_cast<std::vector<double>(*)(const Node&)> (Ghost::Qbnds.find(type)->second);
      auto Qbnd = DVector(field(fn));
      auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
      auto dQint = elems[g.elm_id]->computational->dQfp[dir][fn.right];
      auto J = elems[g.elm_id]->J[dir + 1][fn.right]; // J[0] for snodes, J[1] for csi-fnodes and J[2] for eta-fnodes
      // Inplace
      std::transform(Qbnd.begin(), Qbnd.end(), Qbnd.begin(), [&J](auto &c) { return c / J; });
      // Ghost's Physical Properties
      g.computational->Qfp[dir][fn.local] = 2.0 * Qbnd - Qint;
      g.computational->dQfp[dir][fn.local] = dQint;
    }
    break;
  // Supersonic Outlet BC (Neumann BC)
  case 2:
    for (auto &fn : g.fnodes)
    {
      auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
      auto dQint = elems[g.elm_id]->computational->dQfp[dir][fn.right];
      // auto dQbnd = Ghost::dQbnds.search(g.group);

      g.computational->Qfp[dir][fn.local] = Qint;
      // g.computational->dQfp[fn.local] = dQint + 2.0*();
    }
    break;
  }
}

// Interpolate solution to a node in computational domain:
template <typename Equation>
DVector SD<Equation>::interpolate_solution_to_node(std::shared_ptr<Element> &e, const Node &n)
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(this->order);

  Helpers<Lagrange>::init();
  Helpers<Lagrange>::set_nodes(Helpers<GL>::get_nodes());

  double csi = n.coords[0], eta = n.coords[1];
  if (std::abs(csi)>1.0||std::abs(eta)>1.0) {
    std::cout << "(csi, eta) = (" << csi << ", " << eta << ")\n";
    std::cin.get();
  }
  double Lcsi, Leta;
  unsigned int i, j;
  unsigned int s_index;

  DVector solution{};

  s_index = 0;
  for (auto &sp : this->snodes)
  {
    i = (int)s_index / this->order;
    j = s_index % this->order;

    Lcsi = Helpers<Lagrange>::Pn(i, csi);
    Leta = Helpers<Lagrange>::Pn(j, eta);

    solution += ((Lcsi * Leta) * e->computational->Qsp[s_index]);

    s_index++;
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GL>::delete_nodes();

  return solution;
}

template <typename Equation>
DVector SD<Equation>::interpolate_solution_to_fp(std::shared_ptr<Element> &e, const Node &n, long f_index, int dir)
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(this->order);

  Helpers<Lagrange>::init();
  Helpers<Lagrange>::set_nodes(Helpers<GL>::get_nodes());

  double csi = n.coords[0], eta = n.coords[1];
  double Lcsi, Leta;
  unsigned int i, j, fp_i, fp_j;
  unsigned int s_index;

  DVector solution{};

  if (dir==1) 
  {
    fp_i = (unsigned int)f_index / (this->order+1);
    fp_j = (unsigned int)f_index % (this->order+1);
  }
  else
  {
    fp_j = (unsigned int)f_index / (this->order+1);
    fp_i = (unsigned int)f_index % (this->order+1);
  }

  s_index = 0;
  for (auto &sp : this->snodes)
  {
    i = (int)s_index / this->order;
    j = s_index % this->order;

    if (dir==0) 
    {
      if (j==fp_j)
      {
        Lcsi = Helpers<Lagrange>::Pn(i, csi);
        //Leta = Helpers<Lagrange>::Pn(j, eta);

        solution += (Lcsi * e->computational->Qsp[s_index]);
      }
    }
    else 
    {
      if (i==fp_i)
      {
        //Lcsi = Helpers<Lagrange>::Pn(i, csi);
        Leta = Helpers<Lagrange>::Pn(j, eta);

        solution += (Leta * e->computational->Qsp[s_index]);
      }
    }
    
    s_index++;
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GL>::delete_nodes();

  return solution;
}

// ---------------------- //

// 2) INTERPOLATE FROM SOLUTION POINTS TO FLUX POINTS
template <typename Equation>
void SD<Equation>::interpolate_sp2fp(std::shared_ptr<Element> &e)
{
  unsigned int index, f_index;
  index = 0;
  for (auto &vec_lines : this->fnodes)
  {
    f_index = 0;
    for (auto &node : vec_lines) // line nodes for a specific direction (x, y)
    {
      e->computational->Qfp[index][f_index] = 0.0;
      auto vec = this->interpolate_solution_to_fp(e, this->fnodes[index][f_index], f_index, index);
      e->computational->Qfp[index][f_index] = vec;

      f_index++;
    }
    index++;
  }
}

// 3) Calculate Fluxes
// 3.1) Calculate Fluxes at SPs
// Euler
template <>
void SD<Euler>::calculate_fluxes_sp(std::shared_ptr<Element> &e)
{
  double q1, q2, q3, q4, q5;
  double dcsi_dx, dcsi_dy, deta_dx, deta_dy;
  unsigned int s_index;
  
  auto gamma = this->GAMMA;

  s_index = 0;
  // for each solution point, calculate the flux vector
  for (auto &node : this->snodes)
  {
    e->physical->Qsp[s_index] = (1.0/e->J[0][s_index])*e->computational->Qsp[s_index];

    q1 = e->physical->Qsp[s_index][0];
    q2 = e->physical->Qsp[s_index][1];
    q3 = e->physical->Qsp[s_index][2];
    q4 = e->physical->Qsp[s_index][3];

    if (this->dimension + 2 == 4)
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
      
      e->computational->Fcsp[0][s_index - 1] = {dcsi_dx * q2 + dcsi_dy * q3,
                                                dcsi_dx * (q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + dcsi_dy * q3 * q2 / q1,
                                                dcsi_dx * (q2 * q3 / q1) + dcsi_dy * (q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
                                                dcsi_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + dcsi_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};

      e->computational->Fcsp[1][s_index - 1] = {deta_dx * q2 + deta_dy * q3,
                                                deta_dx * (q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + deta_dy * q3 * q2 / q1,
                                                deta_dx * (q2 * q3 / q1) + deta_dy * (q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
                                                deta_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + deta_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};

    }
    else if (this->dimension + 2 == 5)
    {
      //q5 = e->computational->Qsp[s_index][4];

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
    }
    s_index++;
  }
}

// Navier-Stokes
template <>
void SD<NavierStokes>::calculate_fluxes_sp(std::shared_ptr<Element> &e)
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
  for (auto &node : this->snodes)
  {
    s_index++;
    // these Q must be the physical solution, so I divide the computational solution Q by |J|
    q1 = e->physical->Qsp[s_index - 1][0] / e->J[0][s_index - 1];
    q2 = e->physical->Qsp[s_index - 1][1] / e->J[0][s_index - 1];
    q3 = e->physical->Qsp[s_index - 1][2] / e->J[0][s_index - 1];
    q4 = e->physical->Qsp[s_index - 1][3] / e->J[0][s_index - 1];

    dq1_dx = e->physical->dQsp[0][s_index - 1][0];
    dq2_dx = e->physical->dQsp[0][s_index - 1][1];
    dq3_dx = e->physical->dQsp[0][s_index - 1][2];
    dq4_dx = e->physical->dQsp[0][s_index - 1][3];

    dq1_dy = e->physical->dQsp[1][s_index - 1][0];
    dq2_dy = e->physical->dQsp[1][s_index - 1][1];
    dq3_dy = e->physical->dQsp[1][s_index - 1][2];
    dq4_dy = e->physical->dQsp[1][s_index - 1][3];

    if (this->dimension + 2 == 4)
    {
      du_dx = (dq2_dx / q1 - q2 * dq1_dx / (q1 * q1));
      du_dy = (dq2_dy / q1 - q2 * dq1_dy / (q1 * q1));

      dv_dx = (dq3_dx / q1 - q2 * dq1_dx / (q1 * q1));
      dv_dy = (dq3_dy / q1 - q2 * dq1_dx / (q1 * q1));

      txx = 2.0 * this->MU * du_dx + this->MUv * (du_dx + dv_dy);
      tyy = 2.0 * this->MU * dv_dy + this->MUv * (du_dx + dv_dy);
      txy = this->MU * (dv_dx + du_dy);
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
      e->computational->Fcsp[0][s_index - 1] = {e->Ji[1][1][s_index] * q2,
                                                e->Ji[1][1][s_index] * (q2 * q2 / q1 + (this->GAMMA - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
                                                e->Ji[1][1][s_index] * (q2 * q3 / q1),
                                                e->Ji[1][1][s_index] * ((q2 / q1) * (this->GAMMA * q4 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3) / q1))};

      e->computational->Fcsp[0][s_index - 1] = {q2,
                                                q2 * q2 / q1 + (this->GAMMA - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
                                                q2 * q3 / q1,
                                                (q2 / q1) * (this->GAMMA * q4 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3) / q1)};

      e->computational->Fcsp[1][s_index - 1] = {q3,
                                                q3 * q2 / q1,
                                                q3 * q3 / q1 + (this->GAMMA - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
                                                (q3 / q1) * (this->GAMMA * q4 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3) / q1)};

      // Diffusive Flux
      e->computational->Fdsp[0][s_index - 1] = {0.0,
                                                txx,
                                                tyx,
                                                txx * q2 / q1 + tyx * q3 / q1 + this->KAPPA * dT_dx};

      e->computational->Fdsp[1][s_index - 1] = {0.0,
                                                txy,
                                                tyy,
                                                txy * q2 / q1 + tyy * q3 / q1 + this->KAPPA * dT_dy};
    }
    else if (this->dimension + 2 == 5)
    {
      q5 = e->physical->Qsp[s_index - 1][4];

      dq5_dx = e->physical->dQsp[0][s_index - 1][4];

      dq5_dy = e->physical->dQsp[1][s_index - 1][4];

      dq1_dz = e->physical->dQsp[2][s_index - 1][0];
      dq2_dz = e->physical->dQsp[2][s_index - 1][1];
      dq3_dz = e->physical->dQsp[2][s_index - 1][2];
      dq4_dz = e->physical->dQsp[2][s_index - 1][3];
      dq5_dz = e->physical->dQsp[2][s_index - 1][4];

      dw_dx = (dq4_dx / q1 - q2 * dq1_dx / (q1 * q1));
      dw_dy = (dq4_dy / q1 - q2 * dq1_dy / (q1 * q1));

      du_dz = (dq2_dz / q1 - q2 * dq1_dz / (q1 * q1));
      dv_dz = (dq3_dz / q1 - q2 * dq1_dz / (q1 * q1));
      dw_dz = (dq4_dz / q1 - q2 * dq1_dz / (q1 * q1));

      txx = 2.0 * this->MU * du_dx + this->MUv * (du_dx + dv_dy + dw_dz);
      tyy = 2.0 * this->MU * dv_dy + this->MUv * (du_dx + dv_dy + dw_dz);
      tzz = 2.0 * this->MU * dw_dz + this->MUv * (du_dx + dv_dy + dw_dz);

      txy = this->MU * (dv_dx + du_dy);
      ;
      tyx = txy;

      txz = this->MU * (dw_dx + du_dz);
      ;
      tzx = txy;

      tzy = this->MU * (dv_dz + dw_dy);
      ;
      tyz = tzy;

      // Convective Flux
      e->computational->Fcsp[0][s_index - 1] = {q2,
                                                q2 * q2 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1),
                                                q2 * q3 / q1,
                                                q2 * q4 / q1,
                                                (q2 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1)};

      e->computational->Fcsp[1][s_index - 1] = {q3,
                                                q3 * q2 / q1,
                                                q3 * q3 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1),
                                                q3 * q4 / q1,
                                                (q3 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1)};

      e->computational->Fcsp[2][s_index - 1] = {q4,
                                                q4 * q2 / q1,
                                                q4 * q3 / q1,
                                                q4 * q4 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1),
                                                (q4 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1)};

      // Diffussive Flux
      e->computational->Fdsp[0][s_index - 1] = {0.0,
                                                txx,
                                                tyx,
                                                tzx,
                                                txx * q2 / q1 + tyx * q3 / q1 + tzx * q4 / q1 + this->KAPPA * dT_dx};

      e->computational->Fdsp[1][s_index - 1] = {0.0,
                                                txy,
                                                tyy,
                                                tzy,
                                                txy * q2 / q1 + tyy * q3 / q1 + tzz * q4 / q1 + this->KAPPA * dT_dy};

      e->computational->Fdsp[2][s_index - 1] = {0.0,
                                                txz,
                                                tyz,
                                                tzz,
                                                txz * q2 / q1 + tyz * q3 / q1 + tzz * q4 / q1 + this->KAPPA * dT_dz};
    }
  }
}

// 3.2) Calculate Fluxes at Internal FPs
// Euler
template <>
void SD<Euler>::calculate_fluxes_fp(std::shared_ptr<Element> &e)
{
  double q1, q2, q3, q4, q5, p;
  double ds_dx, ds_dy;
  //double dcsi_dx, dcsi_dy, deta_dx, deta_dy;
  unsigned int f_index = 0, index = 0;
  auto gamma = this->GAMMA;

  for (auto &vec_lines : this->fnodes)
  {
    f_index = 0;
    for (auto &node : vec_lines) // line nodes for a specific direction (x, y)
    {
      auto J = e->J[index+1][f_index];
      e->physical->Qfp[index][f_index] =  (1.0/J)*e->computational->Qfp[index][f_index];

      q1 = e->physical->Qfp[index][f_index][0];
      q2 = e->physical->Qfp[index][f_index][1];
      q3 = e->physical->Qfp[index][f_index][2];
      q4 = e->physical->Qfp[index][f_index][3];
      p = (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1);
      if (this->dimension + 2 == 4)
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
        ds_dx = e->Ji[index+1][f_index][2 * (index)];
        ds_dy = e->Ji[index+1][f_index][2 * (index) + 1];
        
        e->computational->Fcfp[index][f_index] = {ds_dx * q2 + ds_dy * q3,
                                                  ds_dx * (q2 * q2 / q1 + p) + ds_dy * q3 * q2 / q1,
                                                  ds_dx * (q2 * q3 / q1) + ds_dy * (q3 * q3 / q1 + p),
                                                  ds_dx * ((q2 / q1) * (q4 + p)) + ds_dy * ((q3 / q1) * (q4 + p))};
        
        // e->computational->Fcfp[index][f_index] = {ds_dx * q2 + ds_dy * q3,
        //                                           ds_dx * (q2 * q2 / q1 + p) + ds_dy * q3 * q2 / q1,
        //                                           ds_dx * (q2 * q3 / q1) + ds_dy * (q3 * q3 / q1 + p),
        //                                           ds_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + ds_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};


      }
      else if (this->dimension + 2 == 5)
      {
       
      }
      f_index++;
    }
    index++;
  }
 
}

// Navier-Stokes
template <>
void SD<NavierStokes>::calculate_fluxes_fp(std::shared_ptr<Element> &e)
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

  unsigned int f_index = 0, index = 0;

  for (auto &vec_lines : this->fnodes)
  {
    index++;

    f_index = 0;
    for (auto &node : vec_lines) // line nodes for a specific direction (x, y)
    {
      f_index++;

      q1 = e->physical->Qfp[index - 1][f_index - 1][0];
      q2 = e->physical->Qfp[index - 1][f_index - 1][1];
      q3 = e->physical->Qfp[index - 1][f_index - 1][2];
      q4 = e->physical->Qfp[index - 1][f_index - 1][3];

      dq1_dx = e->physical->dQfp[0][f_index - 1][0];
      dq2_dx = e->physical->dQfp[0][f_index - 1][1];
      dq3_dx = e->physical->dQfp[0][f_index - 1][2];
      dq4_dx = e->physical->dQfp[0][f_index - 1][3];

      dq1_dy = e->physical->dQfp[1][f_index - 1][0];
      dq2_dy = e->physical->dQfp[1][f_index - 1][1];
      dq3_dy = e->physical->dQfp[1][f_index - 1][2];
      dq4_dy = e->physical->dQfp[1][f_index - 1][3];

      if (this->dimension + 2 == 4)
      {
        du_dx = (dq2_dx / q1 - q2 * dq1_dx / (q1 * q1));
        du_dy = (dq2_dy / q1 - q2 * dq1_dy / (q1 * q1));

        dv_dx = (dq3_dx / q1 - q2 * dq1_dx / (q1 * q1));
        dv_dy = (dq3_dy / q1 - q2 * dq1_dx / (q1 * q1));

        txx = 2.0 * this->MU * du_dx + this->MUv * (du_dx + dv_dy);
        tyy = 2.0 * this->MU * dv_dy + this->MUv * (du_dx + dv_dy);
        txy = this->MU * (dv_dx + du_dy);
        tyx = txy;

        // Convective Flux
        e->computational->Fcfp[0][f_index - 1] = {q2,
                                                  q2 * q2 / q1 + (this->GAMMA - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
                                                  q2 * q3 / q1,
                                                  (q2 / q1) * (this->GAMMA * q4 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3) / q1)};

        e->computational->Fcfp[1][f_index - 1] = {q3,
                                                  q3 * q2 / q1,
                                                  q3 * q3 / q1 + (this->GAMMA - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
                                                  (q3 / q1) * (this->GAMMA * q4 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3) / q1)};

        // Diffusive Flux
        e->computational->Fdfp[0][f_index - 1] = {0.0,
                                                  txx,
                                                  tyx,
                                                  txx * q2 / q1 + tyx * q3 / q1 + this->KAPPA * dT_dx};

        e->computational->Fdfp[1][f_index - 1] = {0.0,
                                                  txy,
                                                  tyy,
                                                  txy * q2 / q1 + tyy * q3 / q1 + this->KAPPA * dT_dy};
      }
      else if (this->dimension + 2 == 5)
      {
        q5 = e->physical->Qfp[index - 1][f_index - 1][4];

        dq5_dx = e->physical->dQfp[0][f_index - 1][4];

        dq5_dy = e->physical->dQfp[1][f_index - 1][4];

        dq1_dz = e->physical->dQfp[2][f_index - 1][0];
        dq2_dz = e->physical->dQfp[2][f_index - 1][1];
        dq3_dz = e->physical->dQfp[2][f_index - 1][2];
        dq4_dz = e->physical->dQfp[2][f_index - 1][3];
        dq5_dz = e->physical->dQfp[2][f_index - 1][4];

        dw_dx = (dq4_dx / q1 - q2 * dq1_dx / (q1 * q1));
        dw_dy = (dq4_dy / q1 - q2 * dq1_dy / (q1 * q1));

        du_dz = (dq2_dz / q1 - q2 * dq1_dz / (q1 * q1));
        dv_dz = (dq3_dz / q1 - q2 * dq1_dz / (q1 * q1));
        dw_dz = (dq4_dz / q1 - q2 * dq1_dz / (q1 * q1));

        txx = 2.0 * this->MU * du_dx + this->MUv * (du_dx + dv_dy + dw_dz);
        tyy = 2.0 * this->MU * dv_dy + this->MUv * (du_dx + dv_dy + dw_dz);
        tzz = 2.0 * this->MU * dw_dz + this->MUv * (du_dx + dv_dy + dw_dz);

        txy = this->MU * (dv_dx + du_dy);
        ;
        tyx = txy;

        txz = this->MU * (dw_dx + du_dz);
        ;
        tzx = txy;

        tzy = this->MU * (dv_dz + dw_dy);
        ;
        tyz = tzy;

        // Convective Flux
        e->computational->Fcfp[0][f_index - 1] = {q2,
                                                  q2 * q2 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1),
                                                  q2 * q3 / q1,
                                                  q2 * q4 / q1,
                                                  (q2 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1)};

        e->computational->Fcfp[1][f_index - 1] = {q3,
                                                  q3 * q2 / q1,
                                                  q3 * q3 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1),
                                                  q3 * q4 / q1,
                                                  (q3 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1)};

        e->computational->Fcfp[2][f_index - 1] = {q4,
                                                  q4 * q2 / q1,
                                                  q4 * q3 / q1,
                                                  q4 * q4 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1),
                                                  (q4 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1)};

        // Diffussive Flux
        e->computational->Fdfp[0][f_index - 1] = {0.0,
                                                  txx,
                                                  tyx,
                                                  tzx,
                                                  txx * q2 / q1 + tyx * q3 / q1 + tzx * q4 / q1 + this->KAPPA * dT_dx};

        e->computational->Fdfp[1][f_index - 1] = {0.0,
                                                  txy,
                                                  tyy,
                                                  tzy,
                                                  txy * q2 / q1 + tyy * q3 / q1 + tzz * q4 / q1 + this->KAPPA * dT_dy};

        e->computational->Fdfp[2][f_index - 1] = {0.0,
                                                  txz,
                                                  tyz,
                                                  tzz,
                                                  txz * q2 / q1 + tyz * q3 / q1 + tzz * q4 / q1 + this->KAPPA * dT_dz};
      }
    }
  }
}

// Calculate Fluxes at Ghost Flux Points
// Euler
template <>
void SD<Euler>::calculate_fluxes_fp(Ghost &g, const std::vector<std::shared_ptr<Element>> &elems)
{
  double q1, q2, q3, q4, q5, p;
  double ds_dx, ds_dy;
  unsigned int f_index = 0;
  auto gamma = this->GAMMA;
  
  // ghost face direction (0: csi/x and 1: eta/y)
  int dir = (g.lr_edge == 0 || g.lr_edge == 2) ? 1 : 0;

  // for each ghost flux point, calculate the flux vector
  for (auto &fn : g.fnodes)
  {
    f_index = fn.local;
    auto J = elems[g.elm_id]->J[dir + 1][fn.right];
    g.physical->Qfp[dir][f_index] = g.computational->Qfp[dir][f_index] * (1.0 / J);

    q1 = g.physical->Qfp[dir][f_index][0];
    q2 = g.physical->Qfp[dir][f_index][1];
    q3 = g.physical->Qfp[dir][f_index][2];
    q4 = g.physical->Qfp[dir][f_index][3];
    p = (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1);
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

    ds_dx = elems[g.elm_id]->Ji[dir+1][fn.right][2 * (dir)];
    ds_dy = elems[g.elm_id]->Ji[dir+1][fn.right][2 * (dir) + 1];

    g.computational->Fcfp[dir][f_index] = {ds_dx * q2 + ds_dy * q3,
                                           ds_dx * (q2 * q2 / q1 + p) + ds_dy * q3 * q2 / q1,
                                           ds_dx * (q2 * q3 / q1) + ds_dy * (q3 * q3 / q1 + p),
                                           ds_dx * ((q2 / q1) * (q4 + p)) + ds_dy * ((q3 / q1) * (q4 + p))};

    // g.computational->Fcfp[dir][f_index] = {ds_dx * q2 + ds_dy * q3,
    //                                        ds_dx * (q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + ds_dy * q3 * q2 / q1,
    //                                        ds_dx * (q2 * q3 / q1) + ds_dy * (q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
    //                                        ds_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + ds_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};
    
  }
}

// Euler
template <>
void SD<NavierStokes>::calculate_fluxes_fp(Ghost &g, const std::vector<std::shared_ptr<Element>> &elems)
{
  double q1, q2, q3, q4, q5;
  double dcsi_dx, dcsi_dy, deta_dx, deta_dy;
  unsigned int f_index = 0;

  // ghost face direction (0: csi/x and 1: eta/y)
  int dir = (g.lr_edge == 0 || g.lr_edge == 2) ? 1 : 0;

  // for each solution point, calculate the flux vector
  for (auto &fn : g.fnodes)
  {
    f_index = fn.local;
    auto J = elems[g.elm_id]->J[dir + 1][fn.right];

    g.physical->Qfp[dir][f_index - 1] = g.computational->Qfp[dir][f_index - 1] * (1.0 / J);

    q1 = g.physical->Qfp[dir][f_index - 1][0];
    q2 = g.physical->Qfp[dir][f_index - 1][1];
    q3 = g.physical->Qfp[dir][f_index - 1][2];
    q4 = g.physical->Qfp[dir][f_index - 1][3];
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
    dcsi_dx = elems[g.elm_id]->Ji[0][f_index - 1][0];
    dcsi_dy = elems[g.elm_id]->Ji[0][f_index - 1][1];
    deta_dx = elems[g.elm_id]->Ji[0][f_index - 1][2];
    deta_dy = elems[g.elm_id]->Ji[0][f_index - 1][3];
    auto gamma = this->GAMMA;
    g.computational->Fcfp[0][f_index - 1] = {dcsi_dx * q2 + dcsi_dy * q3,
                                             dcsi_dx * (q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + dcsi_dy * q3 * q2 / q1,
                                             dcsi_dx * (q2 * q3 / q1) + dcsi_dy * (q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
                                             dcsi_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + dcsi_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};

    g.computational->Fcfp[1][f_index - 1] = {deta_dx * q2 + deta_dy * q3,
                                             deta_dx * (q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + deta_dy * q3 * q2 / q1,
                                             deta_dx * (q2 * q3 / q1) + deta_dy * (q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
                                             deta_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + deta_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};

    g.physical->Fcfp[0][f_index - 1] = {q2,
                                        q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
                                        q2 * q3 / q1,
                                        (q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)};

    g.physical->Fcfp[1][f_index - 1] = {q3,
                                        q3 * q2 / q1,
                                        q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
                                        (q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)};
  }
}

// 4) RIEMANN SOLVER
template <>
void SD<Euler>::riemann_solver(std::shared_ptr<Element> &e, const std::vector<std::shared_ptr<Element>> &elems, const std::vector<Ghost> &ghosts)
{
  double rhoL, uL, vL, EL, pL, hL;
  double rhoR, uR, vR, ER, pR, hR;
  double q1, q2, q3, q4;
  double J_L, J_R, nxL, nxR, nyL, nyR, ds_dxL, ds_dyL;
  double ds_dx, ds_dy;
  double rho_ROE, ux_ROE, uy_ROE, h_ROE, ek_ROE, un_ROE, c_ROE;
  double nx, ny;
  double l1, l2, l3, l4;
  double alp1, alp2, alp3, alp4;
  DVector central_term, upwind, fxL, fyL, fxR, fyR;
  int dirL=0, dirR=0;
  int ll_edge=0;
  int signN = 1, signT = 1;
  auto gamma = this->GAMMA;
  std::vector<double> normal;
  std::vector<std::vector<double>> T;

  for (auto &ed : e->edges)
  {
    dirL = (ll_edge == 0 || ll_edge == 2) ? 1 : 0; // 0: x, 1: y
    dirR = (ed.lr_edge == 0 || ed.lr_edge == 2) ? 1 : 0; // 0: x, 1: y
    auto signN = ((ll_edge==0||ll_edge==3)) ? -1 : 1;
    
    for (auto &fn : ed.fnodes)
    {
      J_L = e->J[dirL+1][fn.local];
      ds_dxL = e->Ji[dirL+1][fn.local][2 * dirL];
      ds_dyL = e->Ji[dirL+1][fn.local][2 * dirL + 1];
      
      // Physical Primitive Properties for edge's left side
      rhoL = (1.0/J_L)*e->computational->Qfp[dirL][fn.local][0];
      uL = (1.0/J_L)*e->computational->Qfp[dirL][fn.local][1] / rhoL;
      vL = (1.0/J_L)*e->computational->Qfp[dirL][fn.local][2] / rhoL;
      EL = (1.0/J_L)*e->computational->Qfp[dirL][fn.local][3];
      pL = (gamma - 1.0) * (EL - 0.5 * rhoL * (uL * uL + vL * vL));
      hL = (EL + pL) / rhoL;

      // Physical Conserved Properties
      q1 = rhoL;
      q2 = rhoL*uL;
      q3 = rhoL*vL;
      q4 = EL;


      // Physical Fcfp_x
      fxL = DVector{{rhoL*uL,
                     rhoL*uL*uL + pL,
                     rhoL*uL*vL,
                     uL*(EL + pL)}};

      // fxL = DVector{{q2,
      //                q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
      //                q2 * q3 / q1,
      //                (q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)}};
      
      // Physical Fcfp_y
      fyL = DVector{{rhoL*vL,
                     rhoL*uL*vL,
                     rhoL*vL*vL + pL,
                     vL*(EL + pL)}};
      // fyL = DVector{{q3,
      //                q3 * q2 / q1,
      //                q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
      //                (q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)}};


      // Normal unit vector in Physical Space
      normal = e->get_normal_vector(dirL+1, fn.local, ll_edge);
      nx = normal[0];
      ny = normal[1];

      if (ed.right != -1) // if not a boundary
      {
        J_R = elems[ed.right]->J[dirR+1][fn.right];
        // Physical Primitive Properties for edge's right side
        rhoR = (1.0/J_R)*elems[ed.right]->computational->Qfp[dirR][fn.right][0];
        uR = (1.0/J_R)*elems[ed.right]->computational->Qfp[dirR][fn.right][1] / rhoR;
        vR = (1.0/J_R)*elems[ed.right]->computational->Qfp[dirR][fn.right][2] / rhoR;
        ER = (1.0/J_R)*elems[ed.right]->computational->Qfp[dirR][fn.right][3];
      }
      else
      {
        J_R = J_L; // mirror
        dirR = (ghosts[ed.ghost].lr_edge == 0 || ghosts[ed.ghost].lr_edge == 2) ? 1 : 0; // 0: x, 1: y
        // Physical Primitive Properties for edge's right side (when it's a ghost cell)
        rhoR = (1.0/J_R)*ghosts[ed.ghost].computational->Qfp[dirR][fn.right][0];
        uR = (1.0/J_R)*ghosts[ed.ghost].computational->Qfp[dirR][fn.right][1] / rhoR;
        vR = (1.0/J_R)*ghosts[ed.ghost].computational->Qfp[dirR][fn.right][2] / rhoR;
        ER = (1.0/J_R)*ghosts[ed.ghost].computational->Qfp[dirR][fn.right][3];

        // std::cout << "dirR = " << dirR << "\n";
        // std::cout << "ghosts[ed.ghost].computational->Qfp[dirR][fn.right][0] = " << ghosts[ed.ghost].computational->Qfp[dirR][fn.right][0] << "\n";
        // std::cout << "ghosts[ed.ghost].computational->Qfp[dirR][fn.right][1] = " << ghosts[ed.ghost].computational->Qfp[dirR][fn.right][1] << "\n";
        // std::cout << "ghosts[ed.ghost].computational->Qfp[dirR][fn.right][2] = " << ghosts[ed.ghost].computational->Qfp[dirR][fn.right][2] << "\n";
        // std::cout << "ghosts[ed.ghost].computational->Qfp[dirR][fn.right][3] = " << ghosts[ed.ghost].computational->Qfp[dirR][fn.right][3] << "\n";

      }
      
      pR = (gamma - 1.0) * (ER - 0.5 * rhoR * (uR * uR + vR * vR));
      hR = (ER + pR) / rhoR;

      // Physical Conserved Properties
      q1 = rhoR;
      q2 = rhoR*uR;
      q3 = rhoR*vR;
      q4 = ER;

      // Physical Fcfp_x
      fxR = DVector{{rhoR*uR,
                     rhoR*uR*uR + pR,
                     rhoR*uR*vR,
                     uR*(ER + pR)}};
      // fxR = DVector{{q2,
      //                q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
      //                q2 * q3 / q1,
      //                (q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)}};

      // Physical Fcfp_y
      fyR = DVector{{rhoR*vR,
                     rhoR*uR*vR,
                     rhoR*vR*vR + pR,
                     vR*(ER + pR)}};
      // fyR = DVector{{q3,
      //                q3 * q2 / q1,
      //                q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
      //                (q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)}};

      // Roe-averages
      rho_ROE = std::sqrt(rhoL)*std::sqrt(rhoR);
      ux_ROE = (uL * std::sqrt(rhoL) + uR * std::sqrt(rhoR)) / (std::sqrt(rhoL)+std::sqrt(rhoR));
      uy_ROE = (vL * std::sqrt(rhoL) + vR * std::sqrt(rhoR)) / (std::sqrt(rhoL)+std::sqrt(rhoR));
      h_ROE = (hL * std::sqrt(rhoL) + hR * std::sqrt(rhoR)) / (std::sqrt(rhoL)+std::sqrt(rhoR));
      ek_ROE = (ux_ROE * ux_ROE + uy_ROE * uy_ROE) / 2.0;
      un_ROE = ux_ROE*nx + uy_ROE*ny;
      //ut_ROE = (dirL == 1) ? uy_ROE : ux_ROE;
      c_ROE = std::sqrt((gamma - 1.0) * (h_ROE - ek_ROE));

      // Roe lambdas
      l1 = un_ROE;
      l2 = un_ROE;
      l3 = un_ROE + c_ROE;
      l4 = un_ROE - c_ROE;

      // Harten's Entropy Correction
      //auto delta = c_ROE/10.0;
      if (std::abs(l3) < 0.2) l3 = (l3*l3)/0.4+0.1;
      if (std::abs(l4) < 0.2) l4 = (l4*l4)/0.4+0.1;
      // if (std::abs(l3) > delta) l3 = std::abs(l3);
      // else l3 = (l3*l3 + delta*delta)/(2.0*delta);
      // if (std::abs(l4) > delta) l4 = std::abs(l4);
      // else l4 = (l4*l4 + delta*delta)/(2.0*delta);

      /*
      T = {{  1.0,             0.0,              1.0,                1.0         },
           { ux_ROE,           ny,         ux_ROE+c_ROE*nx,     ux_ROE-c_ROE*nx  },
           { uy_ROE,          -nx,         uy_ROE+c_ROE*ny,     uy_ROE-c_ROE*ny  },
           { ek_ROE, ux_ROE*ny-uy_ROE*nx, h_ROE+c_ROE*un_ROE, h_ROE-c_ROE*un_ROE }};
      */

      // T = {{  1.0,             0.0,              1.0,                1.0         },
      //      { ux_ROE,           ny,         ux_ROE+c_ROE*nx,     ux_ROE-c_ROE*nx  },
      //      { uy_ROE,          -nx,         uy_ROE+c_ROE*ny,     uy_ROE-c_ROE*ny  },
      //      { ek_ROE, ux_ROE*ny-uy_ROE*nx, h_ROE+c_ROE*un_ROE, h_ROE-c_ROE*un_ROE }};
      
      T = {{  1.0,             0.0,              1.0,                1.0         },
           { ux_ROE,          -ny,         ux_ROE+c_ROE*nx,     ux_ROE-c_ROE*nx  },
           { uy_ROE,           nx,         uy_ROE+c_ROE*ny,     uy_ROE-c_ROE*ny  },
           { ek_ROE, -ux_ROE*ny+uy_ROE*nx, h_ROE+c_ROE*un_ROE, h_ROE-c_ROE*un_ROE }};

      // alp1 = (rhoR - rhoL) - (pR - pL) / (c_ROE * c_ROE);
      // alp2 = rho_ROE * ((uR*ny - vR*nx) - (uL*ny - vL*nx));
      // alp3 = ((pR - pL) + rho_ROE * c_ROE * ((uR*nx + vR*ny) - (uL*nx + vL*ny))) / (2.0 * c_ROE * c_ROE);
      // alp4 = ((pR - pL) - rho_ROE * c_ROE * ((uR*nx + vR*ny) - (uL*nx + vL*ny))) / (2.0 * c_ROE * c_ROE);

      alp1 = (rhoR - rhoL) - (pR - pL) / (c_ROE * c_ROE);
      alp2 = rho_ROE * ((-uR*ny + vR*nx) - (-uL*ny + vL*nx));
      alp3 = ((pR - pL) + rho_ROE * c_ROE * ((uR*nx + vR*ny) - (uL*nx + vL*ny))) / (2.0 * c_ROE * c_ROE);
      alp4 = ((pR - pL) - rho_ROE * c_ROE * ((uR*nx + vR*ny) - (uL*nx + vL*ny))) / (2.0 * c_ROE * c_ROE);


      // if (dirL == 1) // y
      // {
      //   nx = 0.0;
      //   ny = 1.0;
      //   T = {{1.0, 0.0, 1.0, 1.0},
      //        {ux_ROE, 1.0, ux_ROE, ux_ROE},
      //        {uy_ROE, 0.0, uy_ROE + c_ROE, uy_ROE - c_ROE},
      //        {ek_ROE, ux_ROE, h_ROE + c_ROE * un_ROE, h_ROE - c_ROE * un_ROE}};
      //   alp2 = rho_ROE * (uR - uL);
      //   alp3 = (pR - pL + rho_ROE * c_ROE * (vR - vL)) / (2.0 * c_ROE * c_ROE);
      //   alp4 = (pR - pL - rho_ROE * c_ROE * (vR - vL)) / (2.0 * c_ROE * c_ROE);
      // }
      // else // x
      // {
      //   nx = 1.0;
      //   ny = 0.0;
      //   T = {{1.0, 0.0, 1.0, 1.0},
      //        {ux_ROE, 0.0, ux_ROE + c_ROE, ux_ROE - c_ROE},
      //        {uy_ROE, -1.0, uy_ROE, uy_ROE},
      //        {ek_ROE, -uy_ROE, h_ROE + c_ROE * un_ROE, h_ROE - c_ROE * un_ROE}};
      //   alp2 = rho_ROE * (vR - vL);
      //   alp3 = (pR - pL + rho_ROE * c_ROE * (uR - uL)) / (2.0 * c_ROE * c_ROE);
      //   alp4 = (pR - pL - rho_ROE * c_ROE * (uR - uL)) / (2.0 * c_ROE * c_ROE);
      // }
      // T = {{  1.0,             0.0,              1.0,                1.0         },
      //      { ux_ROE,           ny,         ux_ROE+c_ROE*nx,     ux_ROE-c_ROE*nx  },
      //      { uy_ROE,          -nx,         uy_ROE+c_ROE*ny,     uy_ROE-c_ROE*ny  },
      //      { ek_ROE, ux_ROE*ny-uy_ROE*nx, h_ROE+c_ROE*un_ROE, h_ROE-c_ROE*un_ROE }};

      // Flux reconstruction
      central_term = 0.5 * ((fxL*nx + fyL*ny) + (fxR*nx + fyR*ny));

      // upwind = {
      //     0.5 * abs(l1) * alp1 * (T[0][0] + T[0][1] + T[0][2] + T[0][3]),
      //     0.5 * abs(l2) * alp2 * (T[1][0] + T[1][1] + T[1][2] + T[1][3]),
      //     0.5 * abs(l3) * alp3 * (T[2][0] + T[2][1] + T[2][2] + T[2][3]),
      //     0.5 * abs(l4) * alp4 * (T[3][0] + T[3][1] + T[3][2] + T[3][3]),
      // };
       upwind = {
         0.5*std::abs(l1)*alp1*T[0][0] + 0.5*std::abs(l2)*alp2*T[0][1] + 0.5*std::abs(l3)*alp3*T[0][2] + 0.5*std::abs(l4)*alp4*T[0][3],
         0.5*std::abs(l1)*alp1*T[1][0] + 0.5*std::abs(l2)*alp2*T[1][1] + 0.5*std::abs(l3)*alp3*T[1][2] + 0.5*std::abs(l4)*alp4*T[1][3],
         0.5*std::abs(l1)*alp1*T[2][0] + 0.5*std::abs(l2)*alp2*T[2][1] + 0.5*std::abs(l3)*alp3*T[2][2] + 0.5*std::abs(l4)*alp4*T[2][3],
         0.5*std::abs(l1)*alp1*T[3][0] + 0.5*std::abs(l2)*alp2*T[3][1] + 0.5*std::abs(l3)*alp3*T[3][2] + 0.5*std::abs(l4)*alp4*T[3][3],
      };
      //e->computational->Fcfp[dirL][fn.local] = (central_term - upwind);
      
      //ed.computational->Fcfp[dirL][fn.id] = J_L*(central_term - upwind);
      auto fnx = (central_term - upwind)*nx;
      auto fny = (central_term - upwind)*ny;
      ed.computational->Fcfp[dirL][fn.id] = ds_dxL*fnx + ds_dyL*fny;
      
      // std::cout << "J_L = " << J_L << "\n";
      // std::cout << "rhoL = " << rhoL << "\n";
      // std::cout << "uL = " << uL << "\n";
      // std::cout << "vL = " << vL << "\n";
      // std::cout << "EL = " << EL << "\n";
      // std::cout << "ed.ghost = " << ed.ghost << "\n";
      // std::cout << "\n";
      // std::cout << "J_R = " << J_R << "\n";
      // std::cout << "rhoR = " << rhoR << "\n";
      // std::cout << "uR = " << uR << "\n";
      // std::cout << "vR = " << vR << "\n";
      // std::cout << "ER = " << ER << "\n";
      // std::cout << "\n";
      // std::cout << "fnR[0] = " << J_R*fxR[0]*nx + J_R*fyR[0]*ny << "\n";
      // std::cout << "fnR[1] = " << J_R*fxR[1]*nx + J_R*fyR[1]*ny << "\n";
      // std::cout << "fnR[2] = " << J_R*fxR[2]*nx + J_R*fyR[2]*ny << "\n";
      // std::cout << "fnR[3] = " << J_R*fxR[3]*nx + J_R*fyR[3]*ny << "\n";
      // std::cout << "fnL[0] = " << J_L*fxL[0]*nx + J_L*fyL[0]*ny << "\n";
      // std::cout << "fnL[1] = " << J_L*fxL[1]*nx + J_L*fyL[1]*ny << "\n";
      // std::cout << "fnL[2] = " << J_L*fxL[2]*nx + J_L*fyL[2]*ny << "\n";
      // std::cout << "fnL[3] = " << J_L*fxL[3]*nx + J_L*fyL[3]*ny << "\n";
      // std::cout << "\n";
      // std::cout << "central_term[0] = " << central_term[0] << "\n";
      // std::cout << "central_term[1] = " << central_term[1] << "\n";
      // std::cout << "central_term[2] = " << central_term[2] << "\n";
      // std::cout << "central_term[3] = " << central_term[3] << "\n";
      // std::cout << "upwind[0] = " << upwind[0] << "\n";
      // std::cout << "upwind[1] = " << upwind[1] << "\n";
      // std::cout << "upwind[2] = " << upwind[2] << "\n";
      // std::cout << "upwind[3] = " << upwind[3] << "\n";
      // std::cout << "ed.computational->Fcfp[0] = " << fnx[0] + fny[0] << "\n";
      // std::cout << "ed.computational->Fcfp[1] = " << fnx[1] + fny[1] << "\n";
      // std::cout << "ed.computational->Fcfp[2] = " << fnx[2] + fny[2] << "\n";
      // std::cout << "ed.computational->Fcfp[3] = " << fnx[3] + fny[3] << "\n";
      // std::cin.get();
      // std::cout << "rho_ROE = " << rho_ROE << "\n";
      // std::cout << "ux_ROE = " << ux_ROE << "\n";
      // std::cout << "uy_ROE = " << uy_ROE << "\n";
      // std::cout << "h_ROE = " << h_ROE << "\n";
      // std::cout << "ek_ROE = " << ek_ROE << "\n";
      // std::cout << "un_ROE = " << un_ROE << "\n";
      // std::cout << "c_ROE = " << c_ROE << "\n";
      // std::cout << "nx = " << nx << "\n";
      // std::cout << "ny = " << ny << "\n";

      // std::cout << "abs(l1) = " << abs(l1) << "\n";
      // std::cout << "abs(l2) = " << abs(l2) << "\n";
      // std::cout << "abs(l3) = " << abs(l3) << "\n";
      // std::cout << "abs(l4) = " << abs(l4) << "\n";
      // std::cout << "alp1 = " << alp1 << "\n";
      // std::cout << "alp2 = " << alp2 << "\n";
      // std::cout << "alp3 = " << alp3 << "\n";
      // std::cout << "alp4 = " << alp4 << "\n";
      // std::cout << "T[0][0] = " << T[0][0] << "\n";
      // std::cout << "T[0][1] = " << T[0][1] << "\n";
      // std::cout << "T[0][2] = " << T[0][2] << "\n";
      // std::cout << "T[0][2] = " << T[0][3] << "\n";

      // std::cout << "T[1][0] = " << T[1][0] << "\n";
      // std::cout << "T[1][1] = " << T[1][1] << "\n";
      // std::cout << "T[1][2] = " << T[1][2] << "\n";
      // std::cout << "T[1][2] = " << T[1][3] << "\n";
      
      // std::cout << "T[2][0] = " << T[2][0] << "\n";
      // std::cout << "T[2][1] = " << T[2][1] << "\n";
      // std::cout << "T[2][2] = " << T[2][2] << "\n";
      // std::cout << "T[2][2] = " << T[2][3] << "\n";

      // std::cout << "T[3][0] = " << T[3][0] << "\n";
      // std::cout << "T[3][1] = " << T[3][1] << "\n";
      // std::cout << "T[3][2] = " << T[3][2] << "\n";
      // std::cout << "T[3][2] = " << T[3][3] << "\n";
      // std::cout << "ed.computational->Fcfp[" << dirL << "][" << fn.id << "][0] = " << ed.computational->Fcfp[dirL][fn.id][0] << "\n";
      // std::cout << "ed.computational->Fcfp[" << dirL << "][" << fn.id << "][1] = " << ed.computational->Fcfp[dirL][fn.id][1] << "\n";
      // std::cout << "ed.computational->Fcfp[" << dirL << "][" << fn.id << "][2] = " << ed.computational->Fcfp[dirL][fn.id][2] << "\n";
      // std::cout << "ed.computational->Fcfp[" << dirL << "][" << fn.id << "][3] = " << ed.computational->Fcfp[dirL][fn.id][3] << "\n";
      // std::cin.get();
    }
    ll_edge++;
  }
}

template <>
void SD<NavierStokes>::riemann_solver(std::shared_ptr<Element> &e, const std::vector<std::shared_ptr<Element>> &elems, const std::vector<Ghost> &ghosts)
{
  //pass
}

// 5) INTERPOLATE FROM FLUX POINTS TO SOLUTION POINTS
// Euler
template <>
void SD<Euler>::interpolate_fp2sp(std::shared_ptr<Element> &e)
{
  double csi = 0.0, eta = 0.0;
  double Lcsi, Leta;
  double dLcsi, dLeta;
  unsigned int i, j;
  unsigned int index, s_index, f_index;

  s_index = 0;
  for (auto &node : this->snodes)
  {
    // for each solution point
    // I will calculate the lagrange polynomial at its position
    // interpolate the flux in (x/y) direction from fps
    s_index++;

    csi = node.coords[0];
    eta = node.coords[1];
    auto J = e->J[0][s_index-1];
    index = 0;
    for (auto &vec_lines : this->fnodes)
    {
      index++;

      e->computational->dFcsp[index - 1][s_index - 1] = 0.0;

      f_index = 0;
      for (auto &node : vec_lines)
      {
        // index-1 is related to the flux direction (x/0 or y/1)
        if (index == 1) // csi
        { 
          // csi
          Helpers<GLL>::init();
          Helpers<GLL>::set_nodes(this->order+1);
          Helpers<Lagrange>::init();
          Helpers<Lagrange>::set_nodes(Helpers<GLL>::get_nodes());

          i = (int)f_index % (this->order + 1);
          dLcsi = Helpers<Lagrange>::dPn(i, csi);

          Helpers<Lagrange>::delete_nodes();
          Helpers<GLL>::delete_nodes();

          // eta
          Helpers<GL>::init();
          Helpers<GL>::set_nodes(this->order);

          Helpers<Lagrange>::init();
          Helpers<Lagrange>::set_nodes(Helpers<GL>::get_nodes());

          j = (int)f_index / (this->order+1);
          Leta = Helpers<Lagrange>::Pn(j, eta);
          
          Helpers<Lagrange>::delete_nodes();
          Helpers<GL>::delete_nodes();
          f_index++;

          // std::cout << "dLcsi = " << dLcsi << "\n";
          // std::cout << "Leta = " << Leta << "\n";
          // std::cout << "e->computational->Fcfp[" << index - 1 << "][" << f_index - 1 << "][0] = " << e->computational->Fcfp[index - 1][f_index - 1][0] << "\n";
          // std::cout << "e->computational->Fcfp[" << index - 1 << "][" << f_index - 1 << "][1] = " << e->computational->Fcfp[index - 1][f_index - 1][1] << "\n";
          // std::cout << "e->computational->Fcfp[" << index - 1 << "][" << f_index - 1 << "][2] = " << e->computational->Fcfp[index - 1][f_index - 1][2] << "\n";
          // std::cout << "e->computational->Fcfp[" << index - 1 << "][" << f_index - 1 << "][3] = " << e->computational->Fcfp[index - 1][f_index - 1][3] << "\n";
          // std::cin.get();
          e->computational->dFcsp[index - 1][s_index - 1] += ((dLcsi * Leta) * J *  e->computational->Fcfp[index - 1][f_index - 1]);
        } 
        else // eta
        {
          // csi
          Helpers<GL>::init();
          Helpers<GL>::set_nodes(this->order);

          Helpers<Lagrange>::init();
          Helpers<Lagrange>::set_nodes(Helpers<GL>::get_nodes());

          i = (int)f_index / (this->order+1);
          Lcsi = Helpers<Lagrange>::Pn(i, csi);

          Helpers<Lagrange>::delete_nodes();
          Helpers<GL>::delete_nodes();

          // eta
          Helpers<GLL>::init();
          Helpers<GLL>::set_nodes(this->order+1);
          Helpers<Lagrange>::init();
          Helpers<Lagrange>::set_nodes(Helpers<GLL>::get_nodes());

          j = f_index % (this->order+1);
          dLeta = Helpers<Lagrange>::dPn(j, eta);
          
          Helpers<Lagrange>::delete_nodes();
          Helpers<GLL>::delete_nodes();

          f_index++;
          // std::cout << "Lcsi = " << Lcsi << "\n";
          // std::cout << "dLeta = " << dLeta << "\n";
          // std::cout << "e->computational->Fcfp[" << index - 1 << "][" << f_index - 1 << "][0] = " << e->computational->Fcfp[index - 1][f_index - 1][0] << "\n";
          // std::cout << "e->computational->Fcfp[" << index - 1 << "][" << f_index - 1 << "][1] = " << e->computational->Fcfp[index - 1][f_index - 1][1] << "\n";
          // std::cout << "e->computational->Fcfp[" << index - 1 << "][" << f_index - 1 << "][2] = " << e->computational->Fcfp[index - 1][f_index - 1][2] << "\n";
          // std::cout << "e->computational->Fcfp[" << index - 1 << "][" << f_index - 1 << "][3] = " << e->computational->Fcfp[index - 1][f_index - 1][3] << "\n";
          // std::cin.get();
          e->computational->dFcsp[index - 1][s_index - 1] += ((Lcsi * dLeta) * J * e->computational->Fcfp[index - 1][f_index - 1]);
        }
        //e->physical->dFcsp[index-1][s_index-1] += ((dLcsi*dLeta)*e->physical->Fcfp[index-1][f_index-1]);
      }
    }
  }
  
}

template <>
void SD<NavierStokes>::interpolate_fp2sp(std::shared_ptr<Element> &e)
{
  Helpers<GLL>::init();
  Helpers<GLL>::set_nodes(this->order + 1);

  Helpers<Lagrange>::init();
  Helpers<Lagrange>::set_nodes(Helpers<GLL>::get_nodes());

  double csi = 0.0, eta = 0.0;
  //double Lcsi, Leta;
  double dLcsi, dLeta;
  unsigned int i, j;
  unsigned int index, s_index, f_index;

  s_index = 0;
  for (auto &node : this->snodes)
  {
    // for each solution point
    // I will calculate the lagrange polynomial at its position
    // interpolate the flux in (x/y) direction from fps
    s_index++;

    csi = node.coords[0];
    eta = node.coords[1];

    index = 0;
    for (auto &vec_lines : this->fnodes)
    {
      index++;

      //e.Fsp[index-1][s_index-1] = 0.0;
      e->computational->dFcsp[index - 1][s_index - 1] = 0.0;
      e->computational->dFdsp[index - 1][s_index - 1] = 0.0;

      f_index = 0;
      for (auto &node : vec_lines)
      {
        i = (int)f_index / this->order;
        j = f_index % this->order;
        f_index++;
        //Lcsi = Helpers<Lagrange>::Pn(i, csi);
        //Leta = Helpers<Lagrange>::Pn(j, eta);
        dLcsi = Helpers<Lagrange>::dPn(i, csi);
        dLeta = Helpers<Lagrange>::dPn(j, eta);

        // index-1 is related to the flux direction (x/0 or y/1)
        //e.Fsp[index-1][s_index-1] += Lcsi*Leta*e.Ffp[index-1][f_index-1];
        e->computational->dFcsp[index - 1][s_index - 1] += ((dLcsi * dLeta) * e->computational->Fcfp[index - 1][f_index - 1]);
        //e->physical->dFdsp[index-1][s_index-1] += ((dLcsi*dLeta)*e->physical->Fdfp[index-1][f_index-1]);
      }
    }
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GLL>::delete_nodes();
}


template <>
void SD<Euler>::update_fluxes(std::shared_ptr<Element> &e, std::vector<std::shared_ptr<Element>> &elems)
{
  /* 
      After all fluxes have been corrected by the Riemann Solver at all cells interfaces,
      now one may update these flux points for each cell.
  */
  int local_edge_id=0;
  int dir=0, dirR=0;

  for (auto& ed : e->edges) 
  {
    dir = (local_edge_id == 0 || local_edge_id == 2) ? 1 : 0; // 0: x, 1: y
    dirR = (ed.lr_edge == 0 || ed.lr_edge == 2) ? 1 : 0; // 0: x, 1: y 

    for (auto& fn : ed.fnodes) 
    {
      
      // if (ed.right != -1)
      // {
      //   std::cout << "fR - fL \n";
      //   std::cout << "ed.computational->Fcfp[" << dir << "][" << fn.id << "][0] = " << ed.computational->Fcfp[dir][fn.id][0] << "\n";
      //   std::cout << "ed.computational->Fcfp[" << dir << "][" << fn.id << "][1] = " << ed.computational->Fcfp[dir][fn.id][1] << "\n";
      //   std::cout << "ed.computational->Fcfp[" << dir << "][" << fn.id << "][2] = " << ed.computational->Fcfp[dir][fn.id][2] << "\n";
      //   std::cout << "ed.computational->Fcfp[" << dir << "][" << fn.id << "][3] = " << ed.computational->Fcfp[dir][fn.id][3] << "\n";
      //   std::cout << "eL->Fcfp[" << dir  << "][" << fn.local << "][0] = " << e->computational->Fcfp[dir][fn.local][0] << "\n";
      //   std::cout << "eL->Fcfp[" << dir  << "][" << fn.local << "][1] = " << e->computational->Fcfp[dir][fn.local][1] << "\n";
      //   std::cout << "eL->Fcfp[" << dir  << "][" << fn.local << "][2] = " << e->computational->Fcfp[dir][fn.local][2] << "\n";
      //   std::cout << "eL->Fcfp[" << dir  << "][" << fn.local << "][3] = " << e->computational->Fcfp[dir][fn.local][3] << "\n";
      //   std::cout << "eR->Fcfp[" << dirR << "][" << fn.right << "][0] = " << elems[ed.right]->computational->Fcfp[dirR][fn.right][0] << "\n";
      //   std::cout << "eR->Fcfp[" << dirR << "][" << fn.right << "][1] = " << elems[ed.right]->computational->Fcfp[dirR][fn.right][1] << "\n";
      //   std::cout << "eR->Fcfp[" << dirR << "][" << fn.right << "][2] = " << elems[ed.right]->computational->Fcfp[dirR][fn.right][2] << "\n";
      //   std::cout << "eR->Fcfp[" << dirR << "][" << fn.right << "][3] = " << elems[ed.right]->computational->Fcfp[dirR][fn.right][3] << "\n";
      //   std::cout << " \n";
      //   std::cin.get();
      // }
      // if (ed.ghost != -1) 
      // {
      //   std::cout << "ghostR - fL \n";
      //   std::cout << "ed.computational->Fcfp[" << dir << "][" << fn.id << "][0] = " << ed.computational->Fcfp[dir][fn.id][0] << "\n";
      //   std::cout << "ed.computational->Fcfp[" << dir << "][" << fn.id << "][1] = " << ed.computational->Fcfp[dir][fn.id][1] << "\n";
      //   std::cout << "ed.computational->Fcfp[" << dir << "][" << fn.id << "][2] = " << ed.computational->Fcfp[dir][fn.id][2] << "\n";
      //   std::cout << "ed.computational->Fcfp[" << dir << "][" << fn.id << "][3] = " << ed.computational->Fcfp[dir][fn.id][3] << "\n";
      //   std::cout << "eL->Fcfp[" << dir  << "][" << fn.local << "][0] = " << e->computational->Fcfp[dir][fn.local][0] << "\n";
      //   std::cout << "eL->Fcfp[" << dir  << "][" << fn.local << "][1] = " << e->computational->Fcfp[dir][fn.local][1] << "\n";
      //   std::cout << "eL->Fcfp[" << dir  << "][" << fn.local << "][2] = " << e->computational->Fcfp[dir][fn.local][2] << "\n";
      //   std::cout << "eL->Fcfp[" << dir  << "][" << fn.local << "][3] = " << e->computational->Fcfp[dir][fn.local][3] << "\n";
      //   // std::cout << "eR->Fcfp[" << dirR << "][" << fn.right << "][0] = " << elems[ed.right]->computational->Fcfp[dirR][fn.right][0] << "\n";
      //   // std::cout << "eR->Fcfp[" << dirR << "][" << fn.right << "][1] = " << elems[ed.right]->computational->Fcfp[dirR][fn.right][1] << "\n";
      //   // std::cout << "eR->Fcfp[" << dirR << "][" << fn.right << "][2] = " << elems[ed.right]->computational->Fcfp[dirR][fn.right][2] << "\n";
      //   // std::cout << "eR->Fcfp[" << dirR << "][" << fn.right << "][3] = " << elems[ed.right]->computational->Fcfp[dirR][fn.right][3] << "\n";
      //   std::cout << " \n";
      //   std::cin.get();
      // }
        
      e->computational->Fcfp[dir][fn.local] = ed.computational->Fcfp[dir][fn.id];
      
      if (ed.right != -1)
        elems[ed.right]->computational->Fcfp[dirR][fn.right] = ed.computational->Fcfp[dir][fn.id];
    }
      
      
    local_edge_id++;
  }
}

template <>
void SD<NavierStokes>::update_fluxes(std::shared_ptr<Element> &e, std::vector<std::shared_ptr<Element>> &elems)
{
  // PASS
}

// 6) RESIDUE
// Euler
template <>
void SD<Euler>::residue(std::shared_ptr<Element> &e)
{
  unsigned int s_index;
  //double ds = 0.0;

  s_index = 0;
  for (auto &node : this->snodes)
  {
    e->computational->res[s_index] = 0.0;
    for (auto dir=0; dir < this->dimension; dir++)
    {
      //ds = e->Ji[0][s_index][3 * dir];
      // std::cout << "e->computational->dFcsp[" << dir << "][" << s_index << "][0] = " << e->computational->dFcsp[dir][s_index][0] << "\n";
      // std::cout << "e->computational->dFcsp[" << dir << "][" << s_index << "][1] = " << e->computational->dFcsp[dir][s_index][1] << "\n";
      // std::cout << "e->computational->dFcsp[" << dir << "][" << s_index << "][2] = " << e->computational->dFcsp[dir][s_index][2] << "\n";
      // std::cout << "e->computational->dFcsp[" << dir << "][" << s_index << "][3] = " << e->computational->dFcsp[dir][s_index][3] << "\n";
      e->computational->res[s_index] += (-1.0*e->computational->dFcsp[dir][s_index]);
    }
    s_index++;
  }
  // std::cin.get();
}

// Navier-Stokes
template <>
void SD<NavierStokes>::residue(std::shared_ptr<Element> &e)
{
  unsigned int index, s_index;

  s_index = 0;
  for (auto &node : this->snodes)
  {
    s_index++;

    e->computational->res[s_index - 1] = 0.0;
    index = 0;
    for (auto &vec_lines : this->fnodes)
    {
      index++;
      e->computational->res[s_index - 1] += (-e->computational->dFcsp[index - 1][s_index - 1] + (this->M / this->Re) * e->computational->dFdsp[index - 1][s_index - 1]);
    }
  }
}

// 7) SOLVE
template <typename Equation>
void SD<Equation>::solve(std::shared_ptr<Mesh> &mesh)
{
  // std::vector<std::thread> threads;

  // std::cout << "Interpolate solution from SPs to FPs\n";
  // for (auto &e : mesh->elems)
  //   threads.push_back(std::thread(&SD<Equation>::interpolate_sp2fp, this, std::ref(e)));
  // for (auto &th : threads)
  //   th.join();

  // std::cout << "Calculate fluxes at the FPs\n";
  // for (auto &e : mesh->elems)
  //   threads.push_back(std::thread(&SD<Equation>::calculate_fluxes_fp, this, std::ref(e)));
  // for (auto &th : threads)
  //   th.join();


  // std::cout << "Apply boundary conditions\n";
  // for (auto &g : mesh->ghosts)
  //   threads.push_back(std::thread(
  //     &SD<Equation>::boundary_condition, 
  //     this, 
  //     std::ref(g), 
  //     std::ref(mesh->elems),
  //     std::ref(mesh->nodes)
  //   ));
  // for (auto &th : threads)
  //   th.join();

  // std::cout << "Calculate fluxes at the ghosts FPs\n";
  // for (auto &g : mesh->ghosts)
  //   threads.push_back(
  //     std::thread(
  //       &SD<Equation>::calculate_fluxes_fp, 
  //       this, 
  //       std::ref(g), 
  //       std::ref(mesh->elems)
  //     )
  //   );
  // for (auto &th : threads)
  //   th.join();

  // // Step 2)
  // std::cout << "Compute Riemann Solver\n";
  // for (auto &e : mesh->elems)
  //   threads.push_back(
  //     std::thread(
  //       &SD<Equation>::riemann_solver, 
  //       this, 
  //       std::ref(e), 
  //       std::ref(mesh->elems), 
  //       std::ref(mesh->ghosts)
  //     )
  //   );
  // for (auto &th : threads)
  //   th.join();

  // std::cout << "Update fluxes\n";
  // for (auto &e : mesh->elems)
  //   threads.push_back(std::thread(&SD<Equation>::update_fluxes, this, std::ref(e)));
  // for (auto &th : threads)
  //   th.join();

  // std::cout << "Interpolate flux derivatives from FPs to SPs\n";
  // for (auto &e : mesh->elems)
  //   threads.push_back(std::thread(&SD<Equation>::interpolate_fp2sp, this, std::ref(e)));
  // for (auto &th : threads)
  //   th.join();

  // std::cout << "Calculate Residue\n";
  // for (auto &e : mesh->elems)
  //   threads.push_back(std::thread(&SD<Equation>::residue, this, std::ref(e)));
  // for (auto &th : threads)
  //   th.join();

  // Step 0)
  for (auto &e : mesh->elems)
  {
    this->interpolate_sp2fp(e);
    //this->calculate_fluxes_sp(e);
    this->calculate_fluxes_fp(e);
  }

  // Step 1)
  for (auto &g : mesh->ghosts)
  {
    this->boundary_condition(g, mesh->elems, mesh->nodes);
    this->calculate_fluxes_fp(g, mesh->elems);
  }

  // Step 2)
  for (auto &e : mesh->elems)
  {
    this->riemann_solver(e, mesh->elems, mesh->ghosts);
  }

  // Step 3)
  for (auto &e : mesh->elems)
  {
    this->update_fluxes(e, mesh->elems);
    this->interpolate_fp2sp(e);
    this->residue(e);
  }
}

// 8) SAVE SOLUTION INTO A VTK FILE
template <typename Equation>
void SD<Equation>::to_vtk(const std::shared_ptr<Mesh> &mesh, const std::string &filename)
{
  std::ofstream output;
  output.open(filename);

  // std::unordered_set<long> visited;
  // std::unordered_map<long, long> fp_dict;
  std::vector<double> coords={};
  std::vector<std::string> header, points, cells, cells_header, cells_types;
  std::vector<std::vector<std::string>> point_data;

  DVector vec={};
  Node n={};

  // long neighbor, global_ind = -1, local_ind = -1;
  long num_elems = 0, num_vertices = 0, sp = 0, num_cpoints = 0;
  long n1, n2, n3, n4, n5, n6, n7, n8;
  // double q1, q2, q3, q4, q5, q6, q7, q8;
  long s_index=0, f_index=0, v_index=0, v_idx=0, which_node=-1;
  double rho=0.0, u=0.0, v=0.0, E=0.0, P=0.0, J=0.0, Jacobian=0.0, Residue=0.0, Mach=0.0;
  double el_num=0.0;

  long pd_index = 0;

 /*
   0.1) All current vertices in the current order
  */
  
  
  long nverts=0;                                       // total number of vertices
  long nsps = mesh->Nel * (this->order * this->order); // total number of solution points
  long nxfps = 0; // total number of x-flux points (global)
  long nyfps = 0; // total number of y-flux points (global)
  long global_id = 0;
  int dir=0, m=0, direction=0;
  long lr_edge;
  long xfp_R, ghost_xfp_R, elem_xfp_R, yfp_R, ghost_yfp_R, elem_yfp_R;
  std::unordered_map<long, long> xfp_map = {};
  std::unordered_map<long, long> yfp_map = {};

  std::unordered_map<long, long> v_map = {};
  for (auto &e : mesh->elems)
  {
    // Mapping vertices
    v_idx=0;
    for (auto n : e->nodes)
    {
      if (v_idx >= e->NUM_VERTICES) break;
      if (v_map.find(n) == v_map.end()) 
      {
        v_map.insert({n, nverts});
        points.push_back(std::to_string(mesh->nodes[n].coords[0]) + std::string(" ") +
                         std::to_string(mesh->nodes[n].coords[1]) + std::string(" ") +
                         std::to_string(mesh->nodes[n].coords[2]));
        nverts++;
      }
      v_idx++;
    }
  }

  /* 
    check if the flux point is over an edge:
    |
    |-- T: check if its neighbor has already been mapped
    |   |
    |   |-- T: get the global index
    |   |
    |   |__ F: global index := nxfps, save this in the map {e.id*(order*(order+1)) + f_index, global index}
    |
    |__ F: global index := nxfps, save this in the map {e.id*(order*(order+1)) + f_index, global index}
    
    nxfps++;

  */  

  // Mapping x-flux points
  for (auto &e : mesh->elems) 
  {
    f_index = 0;
    while (f_index < this->order * (this->order + 1))
    {    
      // check if the flux point is over an edge and get which edge owns it
      dir = -1;
      if (f_index % (this->order+1) == 0)
      {
        dir = 0; // for x-fps dir=0 will mean 4th edge        
        m = (int)(f_index % this->order);
        m = (int)(-m + this->order-1); // reverse order at 3rd edge
        xfp_R = e->edges[3].fnodes[m].right;
        ghost_xfp_R = e->edges[3].ghost;
        elem_xfp_R = e->edges[3].right;
        lr_edge = e->edges[3].lr_edge;
        
      }
      else if (f_index % (this->order+1) == this->order)
      {
        dir = 1;
        m = (int)(f_index % this->order);
        xfp_R = e->edges[1].fnodes[m].right;
        ghost_xfp_R = e->edges[1].ghost;
        elem_xfp_R = e->edges[1].right;
        lr_edge = e->edges[1].lr_edge;
      }
      direction = (lr_edge==0 || lr_edge==2) ? 1 : 0; // 0: x,  1: y

      global_id = nxfps;
      // if dir > 0, then we are at an edge's flux point
      if (dir != -1) 
      {
        // check if its neighbor has already been mapped
        // if the neighbor is an element
        if (ghost_xfp_R == -1) 
        {
          if (direction==0)
          {
            auto search = xfp_map.find(elem_xfp_R*(this->order*(this->order+1))+xfp_R);
            if (search != xfp_map.end())
            {
              global_id = search->second;
              nxfps--;
            }
          }
        }
      }
      xfp_map.insert({e->id*(this->order*(this->order+1))+f_index, global_id});
      f_index++;
      nxfps++;
    }
  }

  // Mapping y-flux points
  for (auto &e : mesh->elems) 
  {
    f_index = 0;
    while (f_index < this->order * (this->order + 1))
    {    
      // check if the flux point is over an edge and get which edge owns it
      dir = -1;
      if (f_index % (this->order+1) == 0)
      {
        dir = 0; // for x-fps dir=0 will mean 4th edge        
        m = (int)(f_index % this->order);
        yfp_R = e->edges[0].fnodes[m].right;
        ghost_yfp_R = e->edges[0].ghost;
        elem_yfp_R = e->edges[0].right;
        lr_edge = e->edges[0].lr_edge;
      }
      else if (f_index % (this->order+1) == this->order)
      {
        dir = 1;
        m = (int)(f_index % this->order);
        m = (int)(-m + this->order-1); // reverse order at 2nd edge
        yfp_R = e->edges[2].fnodes[m].right;
        ghost_yfp_R = e->edges[2].ghost;
        elem_yfp_R = e->edges[2].right;
        lr_edge = e->edges[2].lr_edge;
      }
      direction = (lr_edge==0 || lr_edge==2) ? 1 : 0; // 0: x,  1: y

      global_id = nyfps;
      // if dir > 0, then we are at an edge's flux point
      if (dir != -1) 
      {
        // check if its neighbor has already been mapped
        // if the neighbor is an element
        if (ghost_yfp_R == -1) 
        {
          if (direction==0)
          {
            auto search = xfp_map.find(elem_yfp_R*(this->order*(this->order+1))+yfp_R);
            if (search != xfp_map.end())
            {
              global_id = -(search->second+1); // Add 1 to include the 0 case
              nyfps--;
            }
          }
          else
          {
            auto search = yfp_map.find(elem_yfp_R*(this->order*(this->order+1))+yfp_R);
            if (search != yfp_map.end())
            {
              global_id = search->second;
              nyfps--;
            }
          }
        }
      }
      yfp_map.insert({e->id*(this->order*(this->order+1))+f_index, global_id});
      f_index++;
      nyfps++;
    }
  }

  //auto total_num_nodes = mesh->N + mesh->Nel * ((this->order * this->order) + 2 * (this->order * (this->order + 1)));
  auto total_num_nodes = nverts + nsps + nxfps + nyfps;

  header.push_back("# vtk DataFile Version 3.0");
  header.push_back("PosProc Mesh");
  header.push_back("ASCII");
  header.push_back("DATASET UNSTRUCTURED_GRID");
  header.push_back(std::string{"POINTS "} + std::to_string(total_num_nodes) + std::string{" float"});

  point_data.clear();
  point_data.resize(7); // rho, P, Mach, Jacobian, Residue, Umag, BCtype
  for (auto& pd_vec : point_data)
    pd_vec.resize(total_num_nodes + 3);

  // Point data header
  // rho
  point_data[0][0] = std::string{"POINT_DATA "} + std::to_string(total_num_nodes);
  point_data[0][1] = std::string{"SCALARS rho float 1"};
  point_data[0][2] = std::string{"LOOKUP_TABLE default"};
  // P
  point_data[1][0] = std::string{" "};
  point_data[1][1] = std::string{"SCALARS pressure float 1"};
  point_data[1][2] = std::string{"LOOKUP_TABLE default"};
  // Mach
  point_data[2][0] = std::string{" "};
  point_data[2][1] = std::string{"SCALARS Mach float 1"};
  point_data[2][2] = std::string{"LOOKUP_TABLE default"};
  // Jacobian
  point_data[3][0] = std::string{" "};
  point_data[3][1] = std::string{"SCALARS Jacobian float 1"};
  point_data[3][2] = std::string{"LOOKUP_TABLE default"};
  // Residue
  point_data[4][0] = std::string{" "};
  point_data[4][1] = std::string{"SCALARS Residue float 1"};
  point_data[4][2] = std::string{"LOOKUP_TABLE default"};
  // Umag
  point_data[5][0] = std::string{" "};
  point_data[5][1] = std::string{"SCALARS velocity float 3"};
  point_data[5][2] = std::string{"LOOKUP_TABLE default"};
  // BCtype
  point_data[6][0] = std::string{" "};
  point_data[6][1] = std::string{"SCALARS BCtype int 1"};
  point_data[6][2] = std::string{"LOOKUP_TABLE default"};

  pd_index += nverts+3;

  /*
    0) Writing all vertices:
      0.1) All current vertices in the current order
           [0:mesh->N-1] x y z
    
      0.2) Looping through all elements:
        0.2.1) Print all solution nodes (id based on local enumeration
                                        and the current total number of vertices)
              [mesh->N + e.id*(this->order*this->order)+0 : mesh->N + (mesh->Nel-1)*(this->order*this->order)] x y z
                                        
        0.2.2) Print all x-flux nodes (id based on local enumeration
                                       and the current total number of vertices)
              [mesh->N + (mesh->Nel-1)*(this->order*this->order) + e.id*(this->order*(this->order+1)) : mesh->N + (mesh->Nel-1)*(this->order*this->order) ]

        0.2.3) Print all y-flux nodes (id based on local enumeration
                                       and the current total number of vertices)
              [mesh->N + (mesh->Nel-1)*(this->order*this->order) + e.id*(this->order*(this->order+1)) : mesh->N + (mesh->Nel-1)*(this->order*this->order) + (mesh->Nel-1)*(this->order*(this->order+1)) ]
  */

 
  
  std::unordered_map<long, std::vector<double>> vert_to_comp = {
    {0, {-1.0, -1.0,  0.0}},
    {1, { 1.0, -1.0,  0.0}},
    {2, { 1.0,  1.0,  0.0}},
    {3, {-1.0,  1.0,  0.0}}
  };

  /* 0.2) Looping through all elements */
  for (auto &e : mesh->elems)
  {
    // element vertices
    //auto evertices = std::vector<long>(e->nodes.begin(), e->nodes.begin() + 4);
    // Loop though element vertices
    v_index = 0;
    for (auto n_id : e->nodes)
    {
      if (v_index >= e->NUM_VERTICES) break;
      // Interpolating solution from solution points to all element vertices
      auto& vertice = mesh->nodes[n_id];
      auto v_id = v_map.find(vertice.id)->second + 3;
      if (point_data[0][v_id].size() < 1) 
      {
        rho=0.0;
        u=0.0;
        v=0.0;
        E=0.0;
        P=0.0;
        Jacobian=0.0;
        Residue = 0.0;
        Mach=0.0;
        el_num=0.0;
        for (auto e_id : vertice.elems) 
        {
          auto& el = mesh->elems[e_id];
          which_node = -1;
          
          for (auto k=0; k<el->NUM_VERTICES; k++)
            if (el->nodes[k] == vertice.id) which_node=k;
          // std::cout << "k=" << which_node << "\n";
          if (which_node<0) 
          {
            std::cout << "which_node :" << which_node << "\n";
            std::cout << "Vertice ["<<vertice.id<<"]\n";
            std::cin.get();
          }

          // Find vertice computational coordinates
          
          coords = vert_to_comp.find(which_node)->second;
          // std::cout << "Coords: x = " << coords[0] << "; y = " << coords[1] << "\n";
          vec = this->interpolate_solution_to_node(el, Node{coords[0], coords[1], coords[2]});
          J = el->calculate_jacobian_at_node(Node{coords[0], coords[1], coords[2]}, mesh->nodes);
          if (J<=0.0) 
          {
            std::cout << "Jacobian :" << J << "\n";
            std::cout << "Vertice ["<<vertice.id<<"]\n";
            std::cin.get();
          }
          rho      += vec[0]/J;
          u        += vec[1]/vec[0]; // vec[1]/((vec[0]/J)*J)
          v        += vec[2]/vec[0]; // vec[2]/((vec[0]/J)*J)
          E        += vec[3]/J;
          Jacobian += J;
          el_num   += 1.0;
        }
        u = u/el_num;
        v = v/el_num;
        rho = rho/el_num;
        E = E/el_num;
        Jacobian = Jacobian/el_num;
        Residue = 0.0;
        P = (this->GAMMA-1.0)*(E - rho*(u*u + v*v)/2.0);
        Mach = std::sqrt(u*u + v*v)/std::sqrt(this->GAMMA*P/rho);

        // std::cout << "u :" << u << "\n";
        // std::cout << "v :" << v << "\n";
        // std::cout << "E :" << E << "\n";
        // std::cout << "P :" << P << "\n";
        // std::cout << "rho :" << rho << "\n";
        // std::cout << "Mach :" << Mach << "\n";
        // if (Mach<=0.0 || isinf(abs(Mach)) || isnan(abs(Mach)))
        // {
        //   std::cin.get();
        // }

        //std::cout << "el_num :" << el_num << "\n";
        
        point_data[0][v_id] = std::to_string(rho);
        point_data[1][v_id] = std::to_string(P);
        point_data[2][v_id] = std::to_string(Mach);
        point_data[3][v_id] = std::to_string(Jacobian);
        point_data[4][v_id] = std::to_string(Residue);
        point_data[5][v_id] = std::to_string(u) + std::string(" ") +
                              std::to_string(v) + std::string(" ") +
                              std::to_string(0.0);
        point_data[6][v_id] = std::to_string(-1);
      }
      v_index++;
    }
  }

  for (auto &e : mesh->elems)
  {
    /*  0.2.1) Print all solution nodes (id based on local enumeration
                                     and the current total number of vertices)
      [mesh->N + e.id*(this->order*this->order)+0 : mesh->N + (mesh->Nel-1)*(this->order*this->order)] x y z */
    s_index = 0;
    while (s_index < this->order * this->order)
    {
      n = e->transform(this->snodes[s_index], mesh->nodes);
      points.push_back(std::to_string(n.coords[0]) + std::string(" ") +
                       std::to_string(n.coords[1]) + std::string(" ") +
                       std::to_string(n.coords[2]));
      
      vec = this->interpolate_solution_to_node(e, this->snodes[s_index]);
      J = e->calculate_jacobian_at_node(this->snodes[s_index], mesh->nodes);
      if (J<=0.0) 
      {
        std::cout << "Jacobian :" << J << "\n";
        std::cout << "this->snodes["<<s_index<<"]\n";
        std::cin.get();
      }
      rho      = vec[0]/J;
      u        = vec[1]/vec[0]; // vec[1]/((vec[0]/J)*J)
      v        = vec[2]/vec[0]; // vec[2]/((vec[0]/J)*J)
      E        = vec[3]/J;
      P        = (this->GAMMA-1.0)*(E - rho*(u*u + v*v)/2.0);
      Jacobian = J;
      Residue  = std::sqrt(
          e->computational->res[s_index][0]*e->computational->res[s_index][0]+
          e->computational->res[s_index][1]*e->computational->res[s_index][1]+
          e->computational->res[s_index][2]*e->computational->res[s_index][2]+
          e->computational->res[s_index][3]*e->computational->res[s_index][3]
      );
      Mach     = std::sqrt(u*u + v*v)/std::sqrt(this->GAMMA*P/rho);

      // std::cout << "u :" << u << "\n";
      // std::cout << "v :" << v << "\n";
      // std::cout << "E :" << E << "\n";
      // std::cout << "P :" << P << "\n";
      // std::cout << "rho :" << rho << "\n";
      // std::cout << "Mach :" << Mach << "\n";
      // if (Mach<=0.0 || isinf(abs(Mach)) || isnan(abs(Mach)))
      // {
      //   std::cout << "Solution Point: " << s_index << "\n";
      //   std::cout << "Element: " << e->id << "\n";
      //   std::cin.get();
      // }
      point_data[0][pd_index] = std::to_string(rho);
      point_data[1][pd_index] = std::to_string(P);
      point_data[2][pd_index] = std::to_string(Mach);
      point_data[3][pd_index] = std::to_string(Jacobian);
      point_data[4][pd_index] = std::to_string(Residue);
      point_data[5][pd_index] = std::to_string(u) + std::string(" ") +
                                std::to_string(v) + std::string(" ") +
                                std::to_string(0.0);
      point_data[6][pd_index] = std::to_string(-1);
      s_index++;
      pd_index++;
    }
  }
  /* 0.2.2) Print all x-flux nodes (id based on local enumeration
                                      and the current total number of vertices)
              [mesh->N + (mesh->Nel-1)*(this->order*this->order) + e.id*(this->order*(this->order+1)) : mesh->N + (mesh->Nel-1)*(this->order*this->order) ] */

  // Add the cell into the set

  // for (auto &ed : e->edges)
  // {
  //   /*
  //     Given a flux point:
  //     - Check if the flux point is an interface FP
  //     - if it is an interface:
  //       - Get the neighbor cell id
  //       - Check if the neighbor cell is inside the set of visited
  //       - if it is inside, it has been visited
  //         - So look up for its sibling FP id at the unordered_map of FPs
  //       - else
  //         - fp_counter ++;
  //         - Add the new fp to the unordered_map of FPs
  //         - add
  //     - else
  //       - fp_counter ++;
  //         - flux point id = fp_counter;
  //   */

  //   dir
  //       global_ind++;
  //   neighbor = ed.right;
  //   if (visited.find(neighbor) != visited.end())
  //   {

  //     local_ind = fp_dict.find(global_ind).second();
  //     auto n = e->transform(this->fnodes[dir][f_index]);
  //     points.push_back(std::to_string(n.coords[0]) + std::string(" ") +
  //                      std::to_string(n.coords[1]) + std::string(" ") +
  //                      std::to_string(n.coords[2]));
  //   }
  //   else
  //   {
  //     fp_counter++;
  //   }
  // }
  for (auto &e : mesh->elems)
  {
    f_index = 0;
    while (f_index < this->order * (this->order + 1))
    {
      global_id = nverts + 3 + nsps + xfp_map.find(e->id*(this->order*(this->order+1))+f_index)->second;
      if (point_data[0][global_id].size() < 1) 
      {
        n = e->transform(this->fnodes[0][f_index], mesh->nodes);
        points.push_back(std::to_string(n.coords[0]) + std::string(" ") +
                         std::to_string(n.coords[1]) + std::string(" ") +
                         std::to_string(n.coords[2]));
        
        vec = this->interpolate_solution_to_fp(e, this->fnodes[0][f_index], f_index, 0);
        J = e->calculate_jacobian_at_node(this->fnodes[0][f_index], mesh->nodes);
        if (J<=0.0) 
        {
          std::cout << "Jacobian :" << J << "\n";
          std::cout << "this->fnodes[0]["<<f_index<<"]\n";
          std::cin.get();
        }

        rho      = vec[0]/J;
        u        = vec[1]/vec[0]; // vec[1]/((vec[0]/J)*J)
        v        = vec[2]/vec[0]; // vec[2]/((vec[0]/J)*J)
        E        = vec[3]/J;
        P        = (this->GAMMA-1.0)*(E - rho*(u*u + v*v)/2.0);
        Jacobian = J;
        Residue = 0.0;
        Mach     = std::sqrt(u*u + v*v)/std::sqrt(this->GAMMA*P/rho);
        
        // std::cout << "u :" << u << "\n";
        // std::cout << "v :" << v << "\n";
        // std::cout << "E :" << E << "\n";
        // std::cout << "P :" << P << "\n";
        // std::cout << "rho :" << rho << "\n";
        // std::cout << "Mach :" << Mach << "\n";
        // if (Mach<=0.0 || isinf(abs(Mach)) || isnan(abs(Mach)))
        // {
        //   std::cin.get();
        // }
        auto BC = -1;
        if (f_index % (this->order + 1) == 0) 
        {
          if (e->edges[3].ghost != -1) 
          {
            BC = mesh->ghosts[e->edges[3].ghost].type;
          }
        }
        else if (f_index % (this->order + 1) == 2) 
        {
          if (e->edges[1].ghost != -1) 
          {
            BC = mesh->ghosts[e->edges[1].ghost].type;
          }
        }
        point_data[0][global_id] = std::to_string(rho);
        point_data[1][global_id] = std::to_string(P);
        point_data[2][global_id] = std::to_string(Mach);
        point_data[3][global_id] = std::to_string(Jacobian);
        point_data[4][global_id] = std::to_string(Residue);
        point_data[5][global_id] = std::to_string(u) + std::string(" ") +
                                   std::to_string(v) + std::string(" ") +
                                   std::to_string(0.0);
        point_data[6][global_id] = std::to_string(BC);
      }
      f_index++;
    }
  }

  for (auto &e : mesh->elems)
  {
    /* 0.2.3) Print all y-flux nodes (id based on local enumeration
                                       and the current total number of vertices)
              [mesh->N + (mesh->Nel-1)*(this->order*this->order) + e.id*(this->order*(this->order+1)) : mesh->N + (mesh->Nel-1)*(this->order*this->order) + (mesh->Nel-1)*(this->order*(this->order+1)) ] */
    f_index = 0;
    while (f_index < this->order * (this->order + 1))
    {
      global_id = yfp_map.find(e->id*(this->order*(this->order+1))+f_index)->second;
      if (global_id<0)
      {
        global_id = nverts + 3 + nsps - (global_id+1);
      }
      else
      {
        global_id = nverts + 3 + nsps + nxfps + global_id;
      }
      
      if (point_data[0][global_id].size() < 1) 
      {
        n = e->transform(this->fnodes[1][f_index], mesh->nodes);
        points.push_back(std::to_string(n.coords[0]) + std::string(" ") +
                         std::to_string(n.coords[1]) + std::string(" ") +
                         std::to_string(n.coords[2]));
        
        
        vec = this->interpolate_solution_to_fp(e, this->fnodes[1][f_index], f_index, 1);
        J = e->calculate_jacobian_at_node(this->fnodes[1][f_index], mesh->nodes);
        if (J<=0.0) 
        {
          std::cout << "Jacobian :" << J << "\n";
          std::cin.get();
        }
        rho      = vec[0]/J;
        u        = vec[1]/vec[0]; // vec[1]/((vec[0]/J)*J)
        v        = vec[2]/vec[0]; // vec[2]/((vec[0]/J)*J)
        E        = vec[3]/J;
        P        = (this->GAMMA-1.0)*(E - rho*(u*u + v*v)/2.0);
        Jacobian = J;
        Residue = 0.0;
        Mach     = std::sqrt(u*u + v*v)/std::sqrt(this->GAMMA*P/rho);
        
        // std::cout << "u :" << u << "\n";
        // std::cout << "v :" << v << "\n";
        // std::cout << "E :" << E << "\n";
        // std::cout << "P :" << P << "\n";
        // std::cout << "rho :" << rho << "\n";
        // std::cout << "Mach :" << Mach << "\n";
        // if (Mach<=0.0 || isinf(abs(Mach)) || isnan(abs(Mach)))
        // {
        //   std::cin.get();
        // }
        auto BC = -1;
        if (f_index % (this->order + 1) == 0) 
        {
          if (e->edges[0].ghost != -1) 
          {
            BC = mesh->ghosts[e->edges[0].ghost].type;
          }
        }
        else if (f_index % (this->order + 1) == 2) 
        {
          if (e->edges[2].ghost != -1) 
          {
            BC = mesh->ghosts[e->edges[2].ghost].type;
          }
        }

        point_data[0][global_id] = std::to_string(rho);
        point_data[1][global_id] = std::to_string(P);
        point_data[2][global_id] = std::to_string(Mach);
        point_data[3][global_id] = std::to_string(Jacobian);
        point_data[4][global_id] = std::to_string(Residue);
        point_data[5][global_id] = std::to_string(u) + std::string(" ") +
                                   std::to_string(v) + std::string(" ") +
                                   std::to_string(0.0);
        point_data[6][global_id] = std::to_string(BC);
      }
      f_index++;
    }
  }

  cells_types.push_back(std::string{"CELL_TYPES "} + std::to_string((this->order + 1) * (this->order + 1) * mesh->Nel));

  for (auto &e : mesh->elems)
  {
    // 1) define new elements
    // 1.1) First row:
    //      - Quad 1:
    auto evertices = std::vector<long>(e->nodes.begin(), e->nodes.begin() + 4);
    n1 = v_map.find(mesh->get_closest(e->transform(this->fnodes[1][0], mesh->nodes), evertices))->second; // nearest vertice from 0 y-FP
    global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + 0)->second;                        
    if (global_id<0) 
      n2 = nverts + nsps - (global_id+1);                                                                    // x-FP
    else 
      n2 = nverts + nsps + nxfps + global_id;                                                             // y-FP
    n3 = nverts + e->id * (this->order * this->order) + 0;                                                // SP
    n4 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + 0)->second;             // x-FP

    cells.push_back(std::string{"4 "} +
                    std::to_string(n1) + std::string{" "} +
                    std::to_string(n2) + std::string{" "} +
                    std::to_string(n3) + std::string{" "} +
                    std::to_string(n4));
    cells_types.push_back("9");
    num_cpoints += 5; // 4 points + 1 indice

    //      - Pentagon k:
    for (auto k = 0; k < this->order - 1; k++)
    {
      global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + k * (this->order + 1))->second;
      if (global_id<0)
        n1 = nverts + nsps - (global_id+1);                                                                                           // x-FP
      else
        n1 = nverts + nsps + nxfps + global_id;                                                                                   // y-FP
      global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + (k + 1) * (this->order + 1))->second;
      if (global_id<0)
        n2 = nverts + nsps - (global_id+1);                                                                                           // x-FP
      else
        n2 = nverts + nsps + nxfps + global_id;                                                                                   // y-FP
      n3 = nverts + e->id * (this->order * this->order) + (k + 1) * this->order;                                                  // SP
      n4 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + k + 1)->second;                               // x-FP
      n5 = nverts + e->id * (this->order * this->order) + k * this->order;                                                        // SP

      cells.push_back(std::string{"5 "} +
                      std::to_string(n1) + std::string{" "} +
                      std::to_string(n2) + std::string{" "} +
                      std::to_string(n3) + std::string{" "} +
                      std::to_string(n4) + std::string{" "} +
                      std::to_string(n5));
      cells_types.push_back("7");
      num_cpoints += 6; // 5 points + 1 indice
    }

    //      - Quad 2:
    evertices = std::vector<long>(e->nodes.begin(), e->nodes.begin() + 4);
    global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + ((this->order + 1) * (this->order - 1)))->second;
    if (global_id<0)
      n1 = nverts + nsps - (global_id+1);                                                                                                           // x-FP
    else
      n1 = nverts + nsps + nxfps + global_id;                                                                                                   // y-FP
    n2 = v_map.find(mesh->get_closest(e->transform(this->fnodes[1][(this->order + 1) * (this->order - 1)], mesh->nodes), evertices))->second;   // nearest vertice from 0 y-FP
    n3 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + this->order)->second;                                         // x-FP
    n4 = nverts + e->id * (this->order * this->order) + this->order * (this->order - 1);                                                        // SP

    cells.push_back(std::string{"4 "} +
                    std::to_string(n1) + std::string{" "} +
                    std::to_string(n2) + std::string{" "} +
                    std::to_string(n3) + std::string{" "} +
                    std::to_string(n4));
    cells_types.push_back("9");
    num_cpoints += 5; // 4 points + 1 indice

    // 1.2) Mid rows:
    for (auto m = 0; m < this->order - 1; m++)
    {
      //      - Pentagon m (1):
      n1 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + m * (this->order + 1))->second;       // x-FP
      n2 = nverts + e->id * (this->order * this->order) + m;                                                              // SP
      global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + m + 1)->second;
      if (global_id<0)
        n3 = nverts + nsps - (global_id+1);                                                                                   // x-FP
      else
        n3 = nverts + nsps + nxfps + global_id;                                                                           // y-FP
      n4 = nverts + e->id * (this->order * this->order) + m + 1;                                                          // SP
      n5 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + (m + 1) * (this->order + 1))->second; // x-FP

      cells.push_back(std::string{"5 "} +
                      std::to_string(n1) + std::string{" "} +
                      std::to_string(n2) + std::string{" "} +
                      std::to_string(n3) + std::string{" "} +
                      std::to_string(n4) + std::string{" "} +
                      std::to_string(n5));
      cells_types.push_back("7");
      num_cpoints += 6; // 5 points + 1 indice

      for (auto k = 0; k < this->order - 1; k++)
      {
        //      - Octogon m-k:
        n1 = nverts + e->id * (this->order * this->order) + (m + k * this->order);                                                          // SP
        n2 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + m * (this->order + 1) + k + 1)->second;               // x-FP
        n3 = nverts + e->id * (this->order * this->order) + (m + (k + 1) * this->order);                                                    // SP

        global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + m + 1 + (k + 1) * (this->order + 1))->second;
        if (global_id<0) 
          n4 = nverts + nsps - (global_id+1);                                                                                                   // x-FP
        else
          n4 = nverts + nsps + nxfps + global_id;                                                                                           // y-FP
        n5 = nverts + e->id * (this->order * this->order) + m + 1 + (k + 1) * this->order;                                                  // SP
        n6 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + (m + 1) * (this->order + 1) + k + 1)->second;         // x-FP
        n7 = nverts + e->id * (this->order * this->order) + m + 1 + k * this->order;                                                        // SP
        global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + m + 1 + k * (this->order + 1))->second;
        if (global_id<0)
          n8 = nverts + nsps - (global_id+1);                                                                                                   // x-FP
        else
          n8 = nverts + nsps + nxfps + global_id;                                                                                           // y-FP

        cells.push_back(std::string{"8 "} +
                        std::to_string(n1) + std::string{" "} +
                        std::to_string(n2) + std::string{" "} +
                        std::to_string(n3) + std::string{" "} +
                        std::to_string(n4) + std::string{" "} +
                        std::to_string(n5) + std::string{" "} +
                        std::to_string(n6) + std::string{" "} +
                        std::to_string(n7) + std::string{" "} +
                        std::to_string(n8));
        cells_types.push_back("7");
        num_cpoints += 9; // 8 points + 1 indice
      }

      //      - Pentagon m (2):
      n1 = nverts + e->id * (this->order * this->order) + (this->order * this->order - this->order + m);                                                // SP
      n2 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + m * (this->order + 1) + this->order)->second;                       // x-FP
      n3 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + (m + 1) * (this->order + 1) + this->order)->second;                 // x-FP
      n4 = nverts + e->id * (this->order * this->order) + (this->order * this->order - this->order + m + 1);                                            // SP
      global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + (this->order + 1) * this->order - this->order + m)->second;
      if (global_id<0)
        n5 = nverts + nsps - (global_id+1);                                                                                                                 // x-FP
      else
        n5 = nverts + nsps + nxfps + global_id;                                                                                                         // y-FP

      cells.push_back(std::string{"5 "} +
                      std::to_string(n1) + std::string{" "} +
                      std::to_string(n2) + std::string{" "} +
                      std::to_string(n3) + std::string{" "} +
                      std::to_string(n4) + std::string{" "} +
                      std::to_string(n5));
      cells_types.push_back("7");
      num_cpoints += 6; // 5 points + 1 indice
    }

    // 1.3) High row:
    //      - Quad 1:
    evertices = std::vector<long>(e->nodes.begin(), e->nodes.begin() + 4);
    n1 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + (this->order + 1) * (this->order - 1))->second; // x-FP
    n2 = nverts + e->id * (this->order * this->order) + (this->order - 1);                                                        // SP

    global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + this->order)->second;
    if (global_id<0)
      n3 = nverts + nsps - (global_id+1);                                                                                             // x-FP
    else
      n3 = nverts + nsps + nxfps + global_id;                                                                                     // y-FP
    n4 = v_map.find(mesh->get_closest(e->transform(this->fnodes[1][this->order], mesh->nodes), evertices))->second;               // nearest vertice from n y-FP

    cells.push_back(std::string{"4 "} +
                    std::to_string(n1) + std::string{" "} +
                    std::to_string(n2) + std::string{" "} +
                    std::to_string(n3) + std::string{" "} +
                    std::to_string(n4));
    cells_types.push_back("9");
    num_cpoints += 5; // 4 points + 1 indice
    //      - Pentagon k:
    for (auto k = 0; k < this->order - 1; k++)
    {
      n1 = nverts + e->id * (this->order * this->order) + (k + 1) * this->order - 1;                                                          // SP
      n2 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + (k + 1) + (this->order - 1) * (this->order + 1))->second; // x-FP
      n3 = nverts + e->id * (this->order * this->order) + (k + 2) * this->order - 1;                                                          // SP
      global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + (k + 2) * (this->order + 1) - 1)->second;
      if (global_id<0)
        n4 = nverts + nsps - (global_id+1);                                                                                                       // x-FP
      else
        n4 = nverts + nsps + nxfps + global_id;                                                                                               // y-FP
      global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + (k + 1) * (this->order + 1) - 1)->second;
      if (global_id<0)
        n5 = nverts + nsps - (global_id+1);                                                                                                       // x-FP
      else 
        n5 = nverts + nsps + nxfps + global_id;                                                                                               // y-FP

      cells.push_back(std::string{"5 "} +
                      std::to_string(n1) + std::string{" "} +
                      std::to_string(n2) + std::string{" "} +
                      std::to_string(n3) + std::string{" "} +
                      std::to_string(n4) + std::string{" "} +
                      std::to_string(n5));
      cells_types.push_back("7");
      num_cpoints += 6; // 5 points + 1 indice
    }

    //      - Quad 2:
    evertices = std::vector<long>(e->nodes.begin(), e->nodes.begin() + 4);
    n1 = nverts + e->id * (this->order * this->order) + this->order * this->order - 1;                                                       // SP
    n2 = nverts + nsps + xfp_map.find(e->id * (this->order * (this->order + 1)) + this->order * (this->order + 1) - 1)->second;              // x-FP
    n3 = v_map.find(mesh->get_closest(e->transform(this->fnodes[0][this->order * (this->order + 1) - 1], mesh->nodes), evertices))->second;  // nearest vertice from n*(n+1)-1 x-FP
    global_id = yfp_map.find(e->id * (this->order * (this->order + 1)) + this->order * (this->order + 1) - 1)->second;
    if (global_id<0)
      n4 = nverts + nsps - (global_id+1);                                                                                                         // x-FP
    else
      n4 = nverts + nsps + nxfps + global_id;                                                                                                 // y-FP

    cells.push_back(std::string{"4 "} +
                    std::to_string(n1) + std::string{" "} +
                    std::to_string(n2) + std::string{" "} +
                    std::to_string(n3) + std::string{" "} +
                    std::to_string(n4));
    cells_types.push_back("9");
    num_cpoints += 5; // 4 points + 1 indice
  }

  cells_header.push_back(std::string{"CELLS "} +
                         std::to_string((this->order + 1) * (this->order + 1) * mesh->Nel) +
                         std::string{" "} +
                         std::to_string(num_cpoints));

  // header, points, point_data, cells, cells_header
  // HEADER
  for (auto& h : header)
    output << h << std::endl;

  // POINTS
  for (auto& p : points)
    output << p << std::endl;

  // CELLS HEADER
  for (auto& ch : cells_header)
    output << ch << std::endl;

  // CELLS
  for (auto& c : cells)
    output << c << std::endl;

  // CELL_TYPES
  for (auto& ct : cells_types)
    output << ct << std::endl;

  // POINT DATA
  for (auto& pd_vec : point_data)
  {
    for (auto& pd : pd_vec)
      output << pd << std::endl;
  }
    
  // CELL_DATA 841
  // SCALARS boundary int 1
  // LOOKUP_TABLE default

  output.close();
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
