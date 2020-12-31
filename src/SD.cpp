#pragma once

#include <algorithm>
#include <fstream>
#include <memory>
#include <string>
#include <sstream>
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

// Auxiliar function
template <typename Equation>
void SD<Equation>::_init_dvec(std::vector<DVector> &vec, size_t num_nodes)
{
  std::vector<double> init_vec(this->dimension + 2, 1.0);

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
  std::vector<double> init_vec(this->dimension + 2, 1.0);

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
void SD<Euler>::initialize_properties(Ghost &g)
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
void SD<NavierStokes>::initialize_properties(Ghost &g)
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
void SD<Euler>::initialize_properties(std::shared_ptr<Element> &e, const std::vector<Vertice> &enodes)
{
  // PHYSICAL
  // Conservative Properties
  this->_init_dvec(e->physical->Qsp, this->snodes.size());
  this->_init_dvec(e->physical->Qfp, this->fnodes[0].size());

  // Convective Fluxes
  this->_init_dvec(e->physical->Fcsp, this->snodes.size());
  this->_init_dvec(e->physical->Fcfp, this->fnodes[0].size());

  // Gradients
  // Convective Fluxes
  this->_init_dvec(e->physical->dFcsp, this->snodes.size());
  this->_init_dvec(e->physical->dFcfp, this->fnodes[0].size());

  // Metrics
  e->allocate_jacobian(this->order);
  e->calculate_jacobian(this->snodes, this->fnodes, enodes);

  std::size_t global_fn = 0;
  std::size_t local_face = 0;
  std::size_t dir = 0, indice = 0;
  std::vector<fNode> fn;

  Node n;
  for (auto &ed : e->edges)
  {
    this->_init_dvec(ed.physical->Qfp, this->order);
    this->_init_dvec(ed.physical->Fcfp, this->order);
    this->_init_dvec(ed.physical->dFcfp, this->order);

    fn.reserve(this->order);
    if (local_face == 0 || local_face == 2)
      dir = 1;
    else
      dir = 0;

    if (local_face == 0 || local_face == 3)
      indice = 0;
    else
      indice = 1;

    for (size_t idx = 0; idx < this->order; idx++)
    {
      global_fn = this->order * indice + (this->order + 1) * idx;
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

  // Gradients
  // Convective Fluxes
  this->_init_dvec(e->computational->dFcsp, this->snodes.size());
  this->_init_dvec(e->computational->dFcfp, this->fnodes[0].size());

  // Residue
  this->_init_dvec(e->computational->res, this->snodes.size());
}

template <>
void SD<NavierStokes>::initialize_properties(std::shared_ptr<Element> &e, const std::vector<Vertice> &enodes)
{
  // PHYSICAL
  // Conservative Properties
  this->_init_dvec(e->physical->Qsp, this->snodes.size());
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

  // Metrics
  e->allocate_jacobian(this->order);
  e->calculate_jacobian(this->snodes, this->fnodes, enodes);

  std::size_t global_fn = 0;
  std::size_t local_face = 0;
  std::size_t dir = 0, indice = 0;
  std::vector<fNode> fn;

  Node n;
  for (auto &ed : e->edges)
  {
    this->_init_dvec(ed.physical->Qfp, this->order);
    this->_init_dvec(ed.physical->dQfp, this->order);
    this->_init_dvec(ed.physical->Fcfp, this->order);
    this->_init_dvec(ed.physical->dFcfp, this->order);
    this->_init_dvec(ed.physical->Fdfp, this->order);
    this->_init_dvec(ed.physical->dFdfp, this->order);

    fn.reserve(this->order);
    if (local_face == 0 || local_face == 2)
      dir = 1;
    else
      dir = 0;

    if (local_face == 0 || local_face == 3)
      indice = 0;
    else
      indice = 1;

    for (size_t idx = 0; idx < this->order; idx++)
    {
      global_fn = this->order * indice + (this->order + 1) * idx;
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
void SD<Equation>::update_edges(std::shared_ptr<Element> &e, std::vector<std::shared_ptr<Element>> &elems, std::vector<Ghost> &ghosts)
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

            if (abs(fn1.coords[0] - fn2.coords[0]) <= 1E-7 &&
                abs(fn1.coords[1] - fn2.coords[1]) <= 1E-7 &&
                abs(fn1.coords[2] - fn2.coords[2]) <= 1E-7)
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
void SD<Equation>::setup(std::shared_ptr<Mesh> &mesh)
{
  /* 
    On setup method, all solution and flux nodes
    are allocated and calculated in create_nodes.

    Then initialize_properties will allocate and
    populate the mesh with an initial value
  */
  this->create_nodes();

  for (auto &e : mesh->elems)
    this->initialize_properties(e, mesh->nodes);

  for (auto &g : mesh->ghosts)
    this->initialize_properties(g);

  // Update edge communication
  for (auto &e : mesh->elems)
    this->update_edges(e, mesh->elems, mesh->ghosts);
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
void SD<Euler>::boundary_condition(Ghost &g, const std::vector<std::shared_ptr<Element>> &elems)
{
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
      g.computational->Qfp[dir][fn.local][0] = elems[g.elm_id]->computational->Qfp[dir][fn.right][0]; // rho
      //  Uwall,t = some value (SLIP - there's no boundary layer here)
      if (dir == 0) // csi
      {
        g.computational->Qfp[dir][fn.local][1] = -elems[g.elm_id]->computational->Qfp[dir][fn.right][1]; // rho*u
        g.computational->Qfp[dir][fn.local][2] = elems[g.elm_id]->computational->Qfp[dir][fn.right][2];  // rho*v
      }
      else // eta
      {
        g.computational->Qfp[dir][fn.local][1] = elems[g.elm_id]->computational->Qfp[dir][fn.right][1];  // rho*u
        g.computational->Qfp[dir][fn.local][2] = -elems[g.elm_id]->computational->Qfp[dir][fn.right][2]; // rho*v
      }
      g.computational->Qfp[dir][fn.local][3] = elems[g.elm_id]->computational->Qfp[dir][fn.right][3]; // E
    }
  // Supersonic Inlet BC or Far Field
  case 1:
    for (auto &fn : g.fnodes)
    {

      auto Qbnd = Ghost::Qbnds[g.group];
      auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
      auto J = elems[g.elm_id]->J[dir + 1][fn.right]; // J[0] for snodes, J[1] for csi-fnodes and J[2] for eta-fnodes
      // Inplace
      std::transform(Qbnd.begin(), Qbnd.end(), Qbnd.begin(), [&J](auto &c) { return c * J; });
      // Ghost's Computational Properties
      g.computational->Qfp[dir][fn.local] = 2.0 * Qbnd - Qint;
    }
  // Supersonic Outlet BC (Neumann BC)
  case 2:
    for (auto &fn : g.fnodes)
    {
      // Ghost's Computational Properties
      auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
      g.computational->Qfp[dir][fn.local] = Qint;
    }
  }
}

template <>
void SD<NavierStokes>::boundary_condition(Ghost &g, const std::vector<std::shared_ptr<Element>> &elems)
{
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
      auto Qbnd = Ghost::Qbnds[g.group];
      auto Qint = elems[g.elm_id]->computational->Qfp[dir][fn.right];
      auto dQint = elems[g.elm_id]->computational->dQfp[dir][fn.right];
      auto J = elems[g.elm_id]->J[dir + 1][fn.right]; // J[0] for snodes, J[1] for csi-fnodes and J[2] for eta-fnodes
      // Inplace
      std::transform(Qbnd.begin(), Qbnd.end(), Qbnd.begin(), [&J](auto &c) { return c * J; });
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

// Interpolation from:
template <typename Equation>
void SD<Equation>::interpolate_interface(Mesh &mesh, std::shared_ptr<Element> &e)
{

  /* 
      PASS
  */
  // Helpers<GL>::init();
  // Helpers<GL>::set_nodes(this->order);

  // Helpers<Lagrange>::init();
  // Helpers<Lagrange>::set_nodes(Helpers<GL>::get_nodes());

  // double x=0.0, y=0.0;
  // double csi=0.0, eta=0.0;
  // double Lcsi, Leta;
  // unsigned int i, j;
  // unsigned int index, s_index, f_index;

  // // node coordinates

  // for (auto& ed: e->edges)
  // {
  //   for ()
  //   {

  //   }
  //   x = node.coords[0];
  //   y = node.coords[1];
  // }

  // // initialize flux nodes solution
  // e->computational->Qfp[f_index-1] = 0.0;
  // e->physical->Qfp[f_index-1] = 0.0;

  // s_index = 0;
  // for (auto& n : this->snodes)
  // {
  //   s_index++;
  //   // s_index = (j+1) + this->order*i;

  //   i = (int) s_index / this->order;
  //   j = s_index % this->order;

  //   Lcsi = Helpers<Lagrange>::Pn(i, csi);
  //   Leta = Helpers<Lagrange>::Pn(j, eta);
  //   e->computational->Qfp[f_index-1] += ((Lcsi*Leta)*e->computational->Qsp[s_index-1]);
  //   e->physical->Qfp[f_index-1] += ((Lcsi*Leta)*e->physical->Qsp[s_index-1]);
  // }

  // index = 0;
  // for (auto& vec_lines : this->fnodes)
  // {
  //   index++;

  //   f_index = 0;
  //   for (auto& node : vec_lines) // line nodes for a specific direction (x, y)
  //   {
  //     f_index++;
  //     // flux node coordinates
  //     csi = node.coords[0];
  //     eta = node.coords[1];

  //     // initialize flux nodes solution
  //     e->computational->Qfp[f_index-1] = 0.0;
  //     e->physical->Qfp[f_index-1] = 0.0;

  //     s_index = 0;
  //     for (auto& n : this->snodes)
  //     {
  //       s_index++;
  //       // s_index = (j+1) + this->order*i;

  //       i = (int) s_index / this->order;
  //       j = s_index % this->order;

  //       Lcsi = Helpers<Lagrange>::Pn(i, csi);
  //       Leta = Helpers<Lagrange>::Pn(j, eta);
  //       e->computational->Qfp[f_index-1] += ((Lcsi*Leta)*e->computational->Qsp[s_index-1]);
  //       e->physical->Qfp[f_index-1] += ((Lcsi*Leta)*e->physical->Qsp[s_index-1]);
  //     }
  //   }
  // }
  // Helpers<Lagrange>::delete_nodes();
  // Helpers<GL>::delete_nodes();
}

// ---------------------- //

// 2) INTERPOLATE FROM SOLUTION POINTS TO FLUX POINTS
template <typename Equation>
void SD<Equation>::interpolate_sp2fp(std::shared_ptr<Element> &e)
{
  Helpers<GL>::init();
  Helpers<GL>::set_nodes(this->order);

  Helpers<Lagrange>::init();
  Helpers<Lagrange>::set_nodes(Helpers<GL>::get_nodes());

  double csi = 0.0, eta = 0.0;
  double Lcsi, Leta;
  unsigned int i, j;
  unsigned int index, s_index, f_index;

  index = 0;
  for (auto &vec_lines : this->fnodes)
  {
    index++;

    f_index = 0;
    for (auto &node : vec_lines) // line nodes for a specific direction (x, y)
    {
      f_index++;
      // flux node coordinates
      csi = node.coords[0];
      eta = node.coords[1];

      // reinitialize flux nodes solution
      e->computational->Qfp[index - 1][f_index - 1] = 0.0;
      e->physical->Qfp[index - 1][f_index - 1] = 0.0;

      s_index = 0;
      for (auto &n : this->snodes)
      {
        s_index++;
        // s_index = (j+1) + this->order*i;

        i = (int)s_index / this->order;
        j = s_index % this->order;

        Lcsi = Helpers<Lagrange>::Pn(i, csi);
        Leta = Helpers<Lagrange>::Pn(j, eta);
        e->computational->Qfp[index - 1][f_index - 1] += ((Lcsi * Leta) * e->computational->Qsp[s_index - 1]);
        e->physical->Qfp[index - 1][f_index - 1] += ((Lcsi * Leta) * e->physical->Qsp[s_index - 1]);
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
void SD<Euler>::calculate_fluxes_sp(std::shared_ptr<Element> &e)
{
  double q1, q2, q3, q4, q5;
  double dcsi_dx, dcsi_dy, deta_dx, deta_dy;
  unsigned int s_index;

  s_index = 0;
  // for each solution point, calculate the flux vector
  for (auto &node : this->snodes)
  {
    s_index++;

    e->physical->Qsp[s_index - 1][0] = e->computational->Qsp[s_index - 1][0] * (1.0 / e->J[0][s_index - 1]);
    e->physical->Qsp[s_index - 1][1] = e->computational->Qsp[s_index - 1][1] * (1.0 / e->J[0][s_index - 1]);
    e->physical->Qsp[s_index - 1][2] = e->computational->Qsp[s_index - 1][2] * (1.0 / e->J[0][s_index - 1]);
    e->physical->Qsp[s_index - 1][3] = e->computational->Qsp[s_index - 1][3] * (1.0 / e->J[0][s_index - 1]);

    q1 = e->physical->Qsp[s_index - 1][0];
    q2 = e->physical->Qsp[s_index - 1][1];
    q3 = e->physical->Qsp[s_index - 1][2];
    q4 = e->physical->Qsp[s_index - 1][3];

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
      dcsi_dx = e->Ji[0][s_index - 1][0];
      dcsi_dy = e->Ji[0][s_index - 1][1];
      deta_dx = e->Ji[0][s_index - 1][2];
      deta_dy = e->Ji[0][s_index - 1][3];
      auto gamma = this->GAMMA;
      e->computational->Fcsp[0][s_index - 1] = {dcsi_dx * q2 + dcsi_dy * q3,
                                                dcsi_dx * (q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + dcsi_dy * q3 * q2 / q1,
                                                dcsi_dx * (q2 * q3 / q1) + dcsi_dy * (q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
                                                dcsi_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + dcsi_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};

      e->computational->Fcsp[1][s_index - 1] = {deta_dx * q2 + deta_dy * q3,
                                                deta_dx * (q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + deta_dy * q3 * q2 / q1,
                                                deta_dx * (q2 * q3 / q1) + deta_dy * (q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
                                                deta_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + deta_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};

      e->physical->Fcsp[0][s_index - 1] = {q2,
                                           q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
                                           q2 * q3 / q1,
                                           (q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)};

      e->physical->Fcsp[1][s_index - 1] = {q3,
                                           q3 * q2 / q1,
                                           q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
                                           (q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)};
    }
    else if (this->dimension + 2 == 5)
    {
      q5 = e->physical->Qsp[s_index - 1][4];

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

      // e->computational->Fcsp[0][s_index-1] = {e->Ji[1][1][s_index]*(q2) + e->Ji[1][2][s_index]*(q3) + e->Ji[1][3][s_index]*(q4),
      //                                         e->Ji[1][1][s_index]*(q2*q2/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[1][2][s_index]*(q3*q2/q1) + e->Ji[1][3][s_index]*(q4*q2/q1),
      //                                         e->Ji[1][1][s_index]*(q2*q3/q1) + e->Ji[1][2][s_index]*(q3*q3/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[1][3][s_index]*(q4*q3/q1),
      //                                         e->Ji[1][1][s_index]*(q2*q4/q1) + e->Ji[1][2][s_index]*(q3*q4/q1) + e->Ji[1][3][s_index]*(q4*q4/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)),
      //                                         e->Ji[1][1][s_index]*((q2/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[1][2][s_index]*((q3/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[1][3][s_index]*((q4/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1))};

      // e->computational->Fcsp[1][s_index-1] = {e->Ji[2][1][s_index]*(q2) + e->Ji[2][2][s_index]*(q3) + e->Ji[2][3][s_index]*(q4),
      //                                         e->Ji[2][1][s_index]*(q2*q2/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[2][2][s_index]*(q3*q2/q1) + e->Ji[2][3][s_index]*(q4*q2/q1),
      //                                         e->Ji[2][1][s_index]*(q2*q3/q1) + e->Ji[2][2][s_index]*(q3*q3/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[2][3][s_index]*(q4*q3/q1),
      //                                         e->Ji[2][1][s_index]*(q2*q4/q1) + e->Ji[2][2][s_index]*(q3*q4/q1) + e->Ji[2][3][s_index]*(q4*q4/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)),
      //                                         e->Ji[2][1][s_index]*((q2/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[2][2][s_index]*((q3/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)) + e->Ji[2][3][s_index]*((q4/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1))};

      // e->computational->Fcsp[2][s_index-1] = {e->Ji[3][1][s_index]* (q2)+e->Ji[3][2][s_index]* (q3)+e->Ji[3][3][s_index]* (q4),
      //                                         e->Ji[3][1][s_index]* (q2 * q2 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1)) + e->Ji[3][2][s_index]* (q3 * q2 / q1) + e->Ji[3][3][s_index]* (q4 * q2 / q1),
      //                                         e->Ji[3][1][s_index]* (q2 * q3 / q1) + e->Ji[3][2][s_index]* (q3 * q3 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1)) + e->Ji[3][3][s_index]* (q4 * q3 / q1),
      //                                         e->Ji[3][1][s_index]* (q2 * q4 / q1) + e->Ji[3][2][s_index]* (q3 * q4 / q1) + e->Ji[3][3][s_index]* (q4 * q4 / q1 + (this->GAMMA - 1.0) * (q5 - 0.5 * (q2 * q2 + q3 * q3 + q4 * q4) / q1)),
      //                                         e->Ji[3][1][s_index]* ((q2 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1)) + e->Ji[3][2][s_index]* ((q3 / q1) * (this->GAMMA * q5 - 0.5 * (this->GAMMA - 1.0) * (q2 * q2 + q3 * q3 + q4 * q4) / q1)) + e->Ji[3][3][s_index]*((q4/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1))};

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
  double q1, q2, q3, q4, q5;
  double ds_dx, ds_dy;
  //double dcsi_dx, dcsi_dy, deta_dx, deta_dy;
  unsigned int f_index = 0, index = 0;
  auto gamma = this->GAMMA;

  for (auto &vec_lines : this->fnodes)
  {
    index++;

    f_index = 0;
    for (auto &node : vec_lines) // line nodes for a specific direction (x, y)
    {
      e->physical->Qfp[index - 1][f_index - 1] = e->computational->Qfp[index - 1][f_index - 1] * (1.0 / e->J[index][f_index - 1]);

      q1 = e->physical->Qfp[index - 1][f_index - 1][0];
      q2 = e->physical->Qfp[index - 1][f_index - 1][1];
      q3 = e->physical->Qfp[index - 1][f_index - 1][2];
      q4 = e->physical->Qfp[index - 1][f_index - 1][3];

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
        // dcsi_dx = e->Ji[1][f_index - 1][0];
        // dcsi_dy = e->Ji[1][f_index - 1][1];
        // deta_dx = e->Ji[2][f_index - 1][2];
        // deta_dy = e->Ji[2][f_index - 1][3];

        ds_dx = e->Ji[index][f_index - 1][2 * (index - 1)];
        ds_dy = e->Ji[index][f_index - 1][2 * (index - 1) + 1];

        e->computational->Fcfp[index - 1][f_index - 1] = {ds_dx * q2 + ds_dy * q3,
                                                          ds_dx * (q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + ds_dy * q3 * q2 / q1,
                                                          ds_dx * (q2 * q3 / q1) + ds_dy * (q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
                                                          ds_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + ds_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};

        e->physical->Fcfp[index - 1][f_index - 1] = {q2,
                                                     q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
                                                     q2 * q3 / q1,
                                                     (q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)};

        // e->computational->Fcfp[0][f_index - 1] = {dcsi_dx * q2 + dcsi_dy * q3,
        //                                           dcsi_dx * (q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + dcsi_dy * q3 * q2 / q1,
        //                                           dcsi_dx * (q2 * q3 / q1) + dcsi_dy * (q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
        //                                           dcsi_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + dcsi_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};

        // e->computational->Fcfp[1][f_index - 1] = {deta_dx * q2 + deta_dy * q3,
        //                                           deta_dx * (q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)) + deta_dy * q3 * q2 / q1,
        //                                           deta_dx * (q2 * q3 / q1) + deta_dy * (q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1)),
        //                                           deta_dx * ((q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)) + deta_dy * ((q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1))};

        // e->physical->Fcfp[0][f_index - 1] = {q2,
        //                                      q2 * q2 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
        //                                      q2 * q3 / q1,
        //                                      (q2 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)};

        // e->physical->Fcfp[1][f_index - 1] = {q3,
        //                                      q3 * q2 / q1,
        //                                      q3 * q3 / q1 + (gamma - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1),
        //                                      (q3 / q1) * (gamma * q4 - 0.5 * (gamma - 1.0) * (q2 * q2 + q3 * q3) / q1)};
      }
      else if (this->dimension + 2 == 5)
      {
        // q5 = e->physical->Qsp[f_index-1][4];

        // e->computational->Fcsp[0][f_index-1] = {q2,
        //                          q2*q2/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
        //                          q2*q3/q1,
        //                          q2*q4/q1,
        //                          (q2/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};

        // e->computational->Fcsp[1][s_index-1] = {q3,
        //                          q3*q2/q1,
        //                          q3*q3/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
        //                          q3*q4/q1,
        //                          (q3/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};

        // e->computational->Fcsp[2][s_index-1] = {q4,
        //                          q4*q2/q1,
        //                          q4*q3/q1,
        //                          q4*q4/q1 + (this->GAMMA-1.0)*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1),
        //                          (q4/q1)*(this->GAMMA*q5 - 0.5*(this->GAMMA-1.0)*(q2*q2 + q3*q3 + q4*q4)/q1)};
      }
    }
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
  double q1, q2, q3, q4, q5;
  double dcsi_dx, dcsi_dy, deta_dx, deta_dy;
  unsigned int f_index = 0;

  // ghost face direction (0: csi/x and 1: eta/y)
  int dir = (g.lr_edge == 0 || g.lr_edge == 2) ? 1 : 0;

  // for each ghost flux point, calculate the flux vector
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
  double rho_ROE, ux_ROE, uy_ROE, h_ROE, ek_ROE, un_ROE, c_ROE;
  double nx, ny;
  double l1, l2, l3, l4;
  double alp1, alp2, alp3, alp4;
  DVector central_term;
  DVector upwind;
  int dir;
  auto gamma = this->GAMMA;
  std::vector<std::vector<double>> T;

  for (auto &ed : e->edges)
  {
    dir = (ed.lr_edge == 0 || ed.lr_edge == 2) ? 1 : 0; // 0: x, 1: y
    for (auto &fn : ed.fnodes)
    {
      rhoL = e->computational->Qfp[dir][fn.local][0];
      uL = e->computational->Qfp[dir][fn.local][1] / rhoL;
      vL = e->computational->Qfp[dir][fn.local][2] / rhoL;
      EL = e->computational->Qfp[dir][fn.local][3];
      pL = (gamma - 1.0) * (EL - 0.5 * rhoL * (uL * uL + vL * vL));
      hL = (EL + pL) / rhoL;
      if (fn.right != -1) // if not a boundary
      {
        rhoR = elems[ed.right]->computational->Qfp[dir][fn.right][0];
        uR = elems[ed.right]->computational->Qfp[dir][fn.right][1] / rhoR;
        vR = elems[ed.right]->computational->Qfp[dir][fn.right][2] / rhoR;
        ER = elems[ed.right]->computational->Qfp[dir][fn.right][3];
        pR = (gamma - 1.0) * (ER - 0.5 * rhoR * (uR * uR + vR * vR));
        hR = (elems[ed.right]->computational->Qfp[dir][fn.right][3] + pR) / rhoR;
        central_term = 0.5 * (e->computational->Fcfp[dir][fn.local] + elems[ed.right]->computational->Fcfp[dir][fn.right]);
      }
      else
      {
        rhoR = ghosts[ed.right].computational->Qfp[dir][fn.local][0];
        uR = ghosts[ed.right].computational->Qfp[dir][fn.local][1] / rhoR;
        vR = ghosts[ed.right].computational->Qfp[dir][fn.local][2] / rhoR;
        ER = ghosts[ed.right].computational->Qfp[dir][fn.local][3];
        pR = (gamma - 1.0) * (ER - 0.5 * rhoR * (uR * uR + vR * vR));
        hR = (ghosts[ed.right].computational->Qfp[dir][fn.local][3] + pR) / rhoR;
        central_term = 0.5 * (e->computational->Fcfp[dir][fn.local] + ghosts[ed.right].computational->Fcfp[dir][fn.right]);
      }

      // Roe-averages
      rho_ROE = std::sqrt(rhoL) * std::sqrt(rhoR);
      ux_ROE = (uL * std::sqrt(rhoL) + uR * std::sqrt(rhoR)) / rho_ROE;
      uy_ROE = (vL * std::sqrt(rhoL) + vR * std::sqrt(rhoR)) / rho_ROE;
      h_ROE = (hL * std::sqrt(rhoL) + hR * std::sqrt(rhoR)) / rho_ROE;
      ek_ROE = (ux_ROE * ux_ROE + uy_ROE * uy_ROE) / 2.0;
      un_ROE = (dir == 1) ? uy_ROE : ux_ROE;
      c_ROE = std::sqrt((gamma - 1.0) * (h_ROE - ek_ROE));

      // Roe lambdas
      l1 = un_ROE;
      l2 = un_ROE;
      l3 = un_ROE + c_ROE;
      l4 = un_ROE - c_ROE;

      /*
      T = {{  1.0,             0.0,              1.0,                1.0         },
           { ux_ROE,           ny,         ux_ROE+c_ROE*nx,     ux_ROE-c_ROE*nx  },
           { uy_ROE,          -nx,         uy_ROE+c_ROE*ny,     uy_ROE-c_ROE*ny  },
           { ek_ROE, ux_ROE*ny-uy_ROE*nx, h_ROE+c_ROE*un_ROE, h_ROE-c_ROE*un_ROE }};
      */
      alp1 = (rhoR - rhoL) - (pR - pL) / (c_ROE * c_ROE);

      if (dir == 1) // y
      {
        nx = 0.0;
        ny = 1.0;
        T = {{1.0, 0.0, 1.0, 1.0},
             {ux_ROE, 1.0, ux_ROE, ux_ROE},
             {uy_ROE, 0.0, uy_ROE + c_ROE, uy_ROE - c_ROE},
             {ek_ROE, ux_ROE, h_ROE + c_ROE * un_ROE, h_ROE - c_ROE * un_ROE}};
        alp2 = rho_ROE * (uR - uL);
        alp3 = (pR - pL + rho_ROE * c_ROE * (vR - vL)) / (2.0 * c_ROE * c_ROE);
        alp4 = (pR - pL - rho_ROE * c_ROE * (vR - vL)) / (2.0 * c_ROE * c_ROE);
      }
      else // x
      {
        nx = 1.0;
        ny = 0.0;
        T = {{1.0, 0.0, 1.0, 1.0},
             {ux_ROE, 0.0, ux_ROE + c_ROE, ux_ROE - c_ROE},
             {uy_ROE, -1.0, uy_ROE, uy_ROE},
             {ek_ROE, -uy_ROE, h_ROE + c_ROE * un_ROE, h_ROE - c_ROE * un_ROE}};
        alp2 = rho_ROE * (vR - vL);
        alp3 = (pR - pL + rho_ROE * c_ROE * (uR - uL)) / (2.0 * c_ROE * c_ROE);
        alp4 = (pR - pL - rho_ROE * c_ROE * (uR - uL)) / (2.0 * c_ROE * c_ROE);
      }

      // Flux reconstruction
      upwind = {
          0.5 * abs(l1) * alp1 * (T[0][0] + T[0][1] + T[0][2] + T[0][3]),
          0.5 * abs(l2) * alp2 * (T[1][0] + T[1][1] + T[1][2] + T[1][3]),
          0.5 * abs(l3) * alp3 * (T[2][0] + T[2][1] + T[2][2] + T[2][3]),
          0.5 * abs(l4) * alp4 * (T[3][0] + T[3][1] + T[3][2] + T[3][3]),
      };
      e->computational->Fcfp[dir][fn.local] = (central_term - upwind);
    }
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

      e->computational->dFcsp[index - 1][s_index - 1] = 0.0;

      f_index = 0;
      for (auto &node : vec_lines)
      {
        f_index++;
        i = (int)f_index / this->order;
        j = f_index % this->order;

        //Lcsi = Helpers<Lagrange>::Pn(i, csi);
        //Leta = Helpers<Lagrange>::Pn(j, eta);
        dLcsi = Helpers<Lagrange>::dPn(i, csi);
        dLeta = Helpers<Lagrange>::dPn(j, eta);

        // index-1 is related to the flux direction (x/0 or y/1)
        e->computational->dFcsp[index - 1][s_index - 1] += ((dLcsi * dLeta) * e->computational->Fcfp[index - 1][f_index - 1]);
        //e->physical->dFcsp[index-1][s_index-1] += ((dLcsi*dLeta)*e->physical->Fcfp[index-1][f_index-1]);
      }
    }
  }
  Helpers<Lagrange>::delete_nodes();
  Helpers<GLL>::delete_nodes();
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
        f_index++;
        i = (int)f_index / this->order;
        j = f_index % this->order;

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

// 6) RESIDUE
// Euler
template <>
void SD<Euler>::residue(std::shared_ptr<Element> &e)
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
      e->computational->res[s_index - 1] += (-e->computational->dFcsp[index - 1][s_index - 1]);
    }
  }
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

  // Step 0)
  for (auto &g : mesh->ghosts)
  {
    this->boundary_condition(g, mesh->elems);
    this->calculate_fluxes_fp(g, mesh->elems);
  }

  // Step 1)
  for (auto &e : mesh->elems)
  {
    this->interpolate_sp2fp(e);
    //this->calculate_fluxes_sp(e);
    this->calculate_fluxes_fp(e);
  }

  // Step 2)
  for (auto &e : mesh->elems)
  {
    //this->calculate_interface_fluxes(e, mesh->elems, mesh->ghosts);
    this->riemann_solver(e, mesh->elems, mesh->ghosts);
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
  std::vector<std::string> header, points, point_data, cells, cells_header, cells_types;
  // long neighbor, global_ind = -1, local_ind = -1;
  long num_elems = 0, num_vertices = 0, sp = 0, num_cpoints = 0;
  long n1, n2, n3, n4, n5, n6, n7, n8;
  // double q1, q2, q3, q4, q5, q6, q7, q8;

  auto nverts = mesh->N;                                      // total number of vertices
  auto nsps = mesh->Nel * (this->order * this->order);        // total number of solution points
  auto nxfps = mesh->Nel * (this->order * (this->order + 1)); // total number of x-flux points (global)
  auto nyfps = nxfps;                                         // total number of y-flux points (global)

  header.push_back("# vtk DataFile Version 3.0");
  header.push_back("PosProc Mesh");
  header.push_back("ASCII");
  header.push_back("DATASET UNSTRUCTURED_GRID");
  header.push_back(std::string{"POINTS "} + std::to_string(mesh->N + mesh->Nel * ((this->order * this->order) + 2 * (this->order * (this->order + 1)))) + std::string{" float"});

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

  /*
   0.1) All current vertices in the current order
*/
  for (auto &n : mesh->nodes)
  {
    points.push_back(std::to_string(n.coords[0]) + std::string(" ") +
                     std::to_string(n.coords[1]) + std::string(" ") +
                     std::to_string(n.coords[2]));
  }

  /* 0.2) Looping through all elements */
  for (auto &e : mesh->elems)
  {
    // Mark element as visited
    // visited.insert(e->id);

    /*  0.2.1) Print all solution nodes (id based on local enumeration
                                     and the current total number of vertices)
      [mesh->N + e.id*(this->order*this->order)+0 : mesh->N + (mesh->Nel-1)*(this->order*this->order)] x y z */
    auto s_index = 0;
    while (s_index < this->order * this->order)
    {
      auto n = e->transform(this->snodes[s_index], mesh->nodes);
      points.push_back(std::to_string(n.coords[0]) + std::string(" ") +
                       std::to_string(n.coords[1]) + std::string(" ") +
                       std::to_string(n.coords[2]));
      s_index++;
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
    auto f_index = 0;
    while (f_index < this->order * (this->order + 1))
    {
      auto n = e->transform(this->fnodes[0][f_index], mesh->nodes);
      points.push_back(std::to_string(n.coords[0]) + std::string(" ") +
                       std::to_string(n.coords[1]) + std::string(" ") +
                       std::to_string(n.coords[2]));
      f_index++;
    }
  }

  for (auto &e : mesh->elems)
  {
    /* 0.2.3) Print all y-flux nodes (id based on local enumeration
                                       and the current total number of vertices)
              [mesh->N + (mesh->Nel-1)*(this->order*this->order) + e.id*(this->order*(this->order+1)) : mesh->N + (mesh->Nel-1)*(this->order*this->order) + (mesh->Nel-1)*(this->order*(this->order+1)) ] */
    auto f_index = 0;
    while (f_index < this->order * (this->order + 1))
    {
      auto n = e->transform(this->fnodes[1][f_index], mesh->nodes);
      points.push_back(std::to_string(n.coords[0]) + std::string(" ") +
                       std::to_string(n.coords[1]) + std::string(" ") +
                       std::to_string(n.coords[2]));
      f_index++;
    }
  }

  cells_types.push_back(std::string{"CELL_TYPES "} + std::to_string((this->order + 1) * (this->order + 1) * mesh->Nel));

  for (auto &e : mesh->elems)
  {
    // 1) define new elements
    // 1.1) First row:
    //      - Quad 1:
    n1 = mesh->get_closest(e->transform(this->fnodes[1][0], mesh->nodes), e->nodes); // nearest vertice from 0 y-FP
    n2 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + 0;      // y-FP
    n3 = nverts + e->id * (this->order * this->order) + 0;                           // SP
    n4 = nverts + nsps + e->id * (this->order * (this->order + 1)) + 0;              // x-FP

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
      n1 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + k * (this->order + 1);       // y-FP
      n2 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + (k + 1) * (this->order + 1); // y-FP
      n3 = nverts + e->id * (this->order * this->order) + (k + 1) * this->order;                            // SP
      n4 = nverts + nsps + e->id * (this->order * (this->order + 1)) + k + 1;                               // x-FP
      n5 = nverts + e->id * (this->order * this->order) + k * this->order;                                  // SP

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
    n1 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + ((this->order + 1) * (this->order - 1));    // y-FP
    n2 = mesh->get_closest(e->transform(this->fnodes[1][(this->order + 1) * (this->order - 1)], mesh->nodes), e->nodes); // nearest vertice from 0 y-FP
    n3 = nverts + nsps + e->id * (this->order * (this->order + 1)) + this->order;                                        // x-FP
    n4 = nverts + e->id * (this->order * this->order) + this->order * (this->order - 1);                                 // SP

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
      n1 = nverts + nsps + e->id * (this->order * (this->order + 1)) + m * (this->order + 1);       // x-FP
      n2 = nverts + e->id * (this->order * this->order) + m;                                        // SP
      n3 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + m + 1;               // y-FP
      n4 = nverts + e->id * (this->order * this->order) + m + 1;                                    // SP
      n5 = nverts + nsps + e->id * (this->order * (this->order + 1)) + (m + 1) * (this->order + 1); // x-FP

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
        n1 = nverts + e->id * (this->order * this->order) + (m + k * this->order);                                    // SP
        n2 = nverts + nsps + e->id * (this->order * (this->order + 1)) + m * (this->order + 1) + k + 1;               // x-FP
        n3 = nverts + e->id * (this->order * this->order) + (m + (k + 1) * this->order);                              // SP
        n4 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + m + 1 + (k + 1) * (this->order + 1); // y-FP
        n5 = nverts + e->id * (this->order * this->order) + m + 1 + (k + 1) * this->order;                            // SP
        n6 = nverts + nsps + e->id * (this->order * (this->order + 1)) + (m + 1) * (this->order + 1) + k + 1;         // x-FP
        n7 = nverts + e->id * (this->order * this->order) + m + 1 + k * this->order;                                  // SP
        n8 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + m + 1 + k * (this->order + 1);       // y-FP

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
      n1 = nverts + e->id * (this->order * this->order) + (this->order * this->order - this->order + m);                          // SP
      n2 = nverts + nsps + e->id * (this->order * (this->order + 1)) + m * (this->order + 1) + this->order;                       // x-FP
      n3 = nverts + nsps + e->id * (this->order * (this->order + 1)) + (m + 1) * (this->order + 1) + this->order;                 // x-FP
      n4 = nverts + e->id * (this->order * this->order) + (this->order * this->order - this->order + m + 1);                      // SP
      n5 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + (this->order + 1) * this->order - this->order + m; // y-FP

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
    n1 = nverts + nsps + e->id * (this->order * (this->order + 1)) + (this->order + 1) * (this->order - 1); // x-FP
    n2 = nverts + e->id * (this->order * this->order) + (this->order - 1);                                  // SP
    n3 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + this->order;                   // y-FP
    n4 = mesh->get_closest(e->transform(this->fnodes[1][this->order], mesh->nodes), e->nodes);              // nearest vertice from n y-FP

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
      n1 = nverts + e->id * (this->order * this->order) + (k + 1) * this->order - 1;                                    // SP
      n2 = nverts + nsps + e->id * (this->order * (this->order + 1)) + (k + 1) + (this->order - 1) * (this->order + 1); // x-FP
      n3 = nverts + e->id * (this->order * this->order) + (k + 2) * this->order - 1;                                    // SP
      n4 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + (k + 2) * (this->order + 1) - 1;         // y-FP
      n5 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + (k + 1) * (this->order + 1) - 1;         // y-FP

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
    n1 = nverts + e->id * (this->order * this->order) + this->order * this->order - 1;                                 // SP
    n2 = nverts + nsps + e->id * (this->order * (this->order + 1)) + this->order * (this->order + 1) - 1;              // x-FP
    n3 = mesh->get_closest(e->transform(this->fnodes[0][this->order * (this->order + 1) - 1], mesh->nodes), e->nodes); // nearest vertice from n*(n+1)-1 x-FP
    n4 = nverts + nsps + nxfps + e->id * (this->order * (this->order + 1)) + this->order * (this->order + 1) - 1;      // y-FP

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
  for (auto &h : header)
  {
    output << h << std::endl;
  }

  // POINTS
  for (auto &p : points)
  {
    output << p << std::endl;
  }

  // POINT DATA
  // for (auto & pd : point_data)
  // {
  //   output << pd << std::endl;
  // }

  // CELLS HEADER
  for (auto &ch : cells_header)
  {
    output << ch << std::endl;
  }

  // CELLS
  for (auto &c : cells)
  {
    output << c << std::endl;
  }

  // CELL_TYPES
  for (auto &ct : cells_types)
  {
    output << ct << std::endl;
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
