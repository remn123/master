#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <iterator>
#include <algorithm>

#include <Element.h>
#include <Helpers.h>
#include <Mesh.h>
#include <Node.h>
#include <Poly.h>

// Type Definitions
typedef std::vector<std::vector<std::vector<double>>> DoubleArr3D;
typedef std::vector<std::vector<double>> DoubleArr2D;

// ------------------------------------------------------------------------ //

// Edge
Edge::Edge(const std::vector<long> &e_nodes,
           const long &edge_id,
           const long &left,
           const long &right)
{
  //num_edges++;
  this->id = edge_id;
  this->boundary = 0;

  // copying edge[edge_id] vector of nodes to the new edge object
  this->nodes = e_nodes;

  // neighbors
  this->left = left;
  this->right = right;

  this->physical = std::make_shared<Property>();
  this->computational = std::make_shared<Property>();

  this->group = -1;
  this->type = -1;
  this->ghost = -1;
  this->lr_edge = -1;

  //if (right < 0)
  //{
  //	this->boundary = 1;
  //}

  //std::cout << "Edge(" << this->id << ") has been created with nodes: ";
  //for (const auto& n : (this->nodes))
  //	std::cout << n << " ";
  //std::cout << std::endl;
}

Edge::~Edge(void)
{
  //std::cout << "Edge (" << this->id <<") has been deleted!" << std::endl;
}

void Edge::print_nodes(void)
{
  std::cout << "Edge(" << this->id << "): ";
  for (const auto &v : (this->nodes))
    std::cout << v << " ";
  std::cout << std::endl;
}
// ------------------------------------------------------------------------ //

// Face
Face::Face(std::vector<std::string> node_list)
{
  //num_faces++;
  //this->id = num_faces;

  // transforming vector<string> to vector<int>
  std::transform(node_list.begin(), node_list.end(), std::back_inserter((*nodes)),
                 [](const std::string &str) { return std::stol(str); });

  //std::cout << "Face(" << this->id << ") has been created with nodes: ";
  //for (const auto& v : (*nodes))
  //	std::cout << v << " ";
  //std::cout << std::endl;
}

Face::~Face(void)
{
  //std::cout << "Face has been deleted!" << std::endl;
}

void Face::print_nodes(void)
{
  std::cout << "Face(" << this->id << "): ";
  for (const auto &v : (*this->nodes))
    std::cout << v << " ";
  std::cout << std::endl;
}
// ------------------------------------------------------------------------ //

/* CONSTRUCTORS */
Element::Element(const std::vector<std::string> &node_list)
{
  num_elems++;
  this->id = num_elems;
  this->boundary = 0;
  this->fringe = 0;
  //this->J = 0.0; // Jacobian
  this->physical = std::make_shared<Property>();
  this->computational = std::make_shared<Property>();
  //std::cout << "I am the (" << this->id << ") Element" << std::endl;

  // transforming vector<string> to vector<int>
  // subtract 1 from the node_id because the first id for GMSH is 1, but as I am using the node id to consult it on a vector
  // I need it to be 0-started.
  std::transform(node_list.begin(), node_list.end(), std::back_inserter(this->nodes),
                 [](const std::string &str) { return std::stol(str) - 1; });

  // this->print_nodes();
}

/* DESTRUCTORS */
Element::~Element(void)
{
  if (Element::num_elems > -1) Element::num_elems--;
}

/* ELEMENT METHODS */
void Element::print_nodes(void)
{
  std::cout << "Element(" << this->id << "): ";
  for (const auto &v : this->nodes)
    std::cout << v << " ";
  std::cout << std::endl;
}

void Element::get_nodes(void)
{
}

bool Element::was_enumerated(const std::vector<long> &vertices)
{
  auto edges_map_iterator = this->edges_map.find(vertices);
  if (edges_map_iterator == this->edges_map.end())
    return false;
  return true;
}

std::vector<long> Triangle::get_ordered_nodes_by_local_edge_id(long edge_id)
{
  //auto num_nodes = nodes.size();
  auto first = this->nodes.begin() + edge_id;
  auto last = first + 2;

  if (first == this->nodes.end())
  {
    last = this->nodes.begin() + 0;
  }

  std::vector<long> edge_vertices(first, last);
  std::sort(edge_vertices.begin(), edge_vertices.end());
  return edge_vertices;
}

void Triangle::enumerate_edges(void)
{
  // 2D Triangle with 3 nodes
  // vertices: 0 1 2
  // edges:
  //		  (0) 0-1
  //		  (1) 1-2
  //		  (2) 2-0
  std::cout << "Enumerating edges of a Triangle" << std::endl;

  //long local_id = 0;
  //while (local_id < NUM_FACES)
  //{
  //	std::vector<long> ordered_vertices(this->get_ordered_nodes_by_local_edge_id(local_id++));
  //
  //	if (!was_enumerated(ordered_vertices))
  //	{
  //		// if edge does not exist yet
  //		edges_map[ordered_vertices]++;
  //	}
  //	ordered_vertices.clear();
  //}
}

double Triangle::calculate_jacobian_at_node(const Node& node)
{
  return 0.0;
}

/* ------------------------------ */
/*            QUADRANGLE          */
/* ------------------------------ */
std::vector<long> Quadrangle::get_nodes_by_local_edge_id(long edge_id, bool sorted)
{
  //auto num_nodes = nodes.size();
  auto first = this->nodes.begin() + edge_id;

  if (first == this->nodes.end() - 1)
  {
    std::vector<long> edge_vertices;
    edge_vertices.push_back(this->nodes.back());
    edge_vertices.push_back(this->nodes.front());
    //for (const auto& v : edge_vertices)
    //	std::cout << v << " ";
    //std::cout << std::endl;
    if (sorted)
    {
      std::sort(edge_vertices.begin(), edge_vertices.end());
    }
    return edge_vertices;
  }
  else
  {
    auto last = first + 2;
    std::vector<long> edge_vertices(first, last);
    //for (const auto& v : edge_vertices)
    //	std::cout << v << " ";
    //std::cout << std::endl;
    if (sorted)
    {
      std::sort(edge_vertices.begin(), edge_vertices.end());
    }
    return edge_vertices;
  }
}

void Quadrangle::enumerate_edges(void)
{
  // 2D Quad with 4 nodes
  // vertices: 0 1 2 3
  // edges:
  //		  (0) 0-1
  //		  (1) 1-2
  //		  (2) 2-3
  //		  (3) 3-0
  //std::cout << "Enumerating edges of a Quadrangle" << std::endl;

  long local_id = 0;
  while (local_id < NUM_FACES)
  {
    std::vector<long> ordered_vertices(this->get_nodes_by_local_edge_id(local_id, true));    // true returns sorted vector
    std::vector<long> unordered_vertices(this->get_nodes_by_local_edge_id(local_id, false)); // false returns unsorted vector

    if (!was_enumerated(ordered_vertices))
    {
      this->num_edges++;
      long edge_id = this->num_edges;
      // if edge does not exist yet
      edges_map.emplace(
          std::make_pair(ordered_vertices, std::vector<long>{edge_id, this->id, -1}));

      //std::cout << "Checking errors in edges_map: odered (" ;
      //for (auto& v : ordered_vertices)
      //{
      //	std::cout << v << " ";
      //}
      //std::cout << ") unordered (";
      //for (auto& v : unordered_vertices)
      //{
      //	std::cout << v << " ";
      //}
      //std::cout << std::endl;
      //for (auto& e : edges_map.find(ordered_vertices)->second)
      //{
      //	std::cout << e << " ";
      //}
      //std::cout << std::endl;

      this->edges.emplace_back(Edge{unordered_vertices, edge_id, this->id, -1});

      //std::cout << "Element(" << this->id << ")->Edge(" << local_id << "/" << edge_id << "/left=" << edges_map.find(ordered_vertices)->second[1] << "/right=" << edges_map.find(ordered_vertices)->second[2] << "): ";
    }
    else // if the edge exist, find its id
    {
      long edge_id = edges_map.find(ordered_vertices)->second[0];
      long neighbor_id = edges_map.find(ordered_vertices)->second[1];
      //std::cout << "Checking errors in edges_map: ordered (";
      //for (auto& v : ordered_vertices)
      //{
      //	std::cout << v << " ";
      //}
      //std::cout << ") unordered (";
      //for (auto& v : unordered_vertices)
      //{
      //	std::cout << v << " ";
      //}
      //std::cout << std::endl;
      this->edges.emplace_back(Edge{unordered_vertices, edge_id, this->id, neighbor_id});

      //std::cout << "Element(" << this->id << ")->Edge(" << local_id << "/" << edge_id << "/left=" << this->id << "/right=" << neighbor_id << "): ";
    }
    ordered_vertices.clear();
    unordered_vertices.clear();
    local_id++;
    //std::cin.get();
  }
  //edges_map.clear();
}

void Quadrangle::allocate_jacobian(int order)
{

  // SPs = order*order
  // FPs = (order+1)*order
  // Jm[0][0:SPs][0:3]
  // Jm[1][0:FPs][0:3] x
  // Jm[2][0:FPs][0:3] y

  // Ji[0][0:SPs][0:3]
  // Ji[1][0:FPs][0:3] x
  // Ji[2][0:FPs][0:3] y

  // this->Jm = DoubleArr3D{3, DoubleArr2D{(order + 1) * order, std::vector<double>{4}}};
  // this->Ji = DoubleArr3D{3, DoubleArr2D{(order + 1) * order, std::vector<double>{4}}};
  this->J = DoubleArr2D(3, std::vector<double>(order * (order + 1), 0.0));
  this->Jm = DoubleArr3D(3, DoubleArr2D((order + 1) * order, std::vector<double>(4, 0.0)));
  this->Ji = DoubleArr3D(3, DoubleArr2D((order + 1) * order, std::vector<double>(4, 0.0)));
}

void Quadrangle::calculate_jacobian(const std::vector<Node> &snodes,
                                    const std::vector<std::vector<Node>> &fnodes,
                                    const std::vector<Vertice> &enodes)
{
  double x1 = enodes[this->nodes[0]].coords[0];
  double x2 = enodes[this->nodes[1]].coords[0];
  double x3 = enodes[this->nodes[2]].coords[0];
  double x4 = enodes[this->nodes[3]].coords[0];

  double y1 = enodes[this->nodes[0]].coords[1];
  double y2 = enodes[this->nodes[1]].coords[1];
  double y3 = enodes[this->nodes[2]].coords[1];
  double y4 = enodes[this->nodes[3]].coords[1];

  double a1, a2, b1, b2, c1, c2, d1, d2, csi, eta;

  //this->J = 0.5 * (x2 - x1) + 0.5 * (y4 - y3);

  a1 = x2 - x1 + x3 - x4;
  b1 = -x2 - x1 + x3 + x4;
  c1 = -x2 + x1 + x3 - x4;
  d1 = x2 + x1 + x3 + x4;
  a2 = y2 - y1 + y3 - y4;
  b2 = -y2 - y1 + y3 + y4;
  c2 = -y2 + y1 + y3 - y4;
  d2 = y2 + y1 + y3 + y4;

  this->metrics = {a1, b1, c1, d1, a2, b2, c2, d2};

  std::size_t index = 0, s_index = 0, f_index = 0;

  for (auto &node : snodes)
  {
    s_index++;
    // flux node coordinates
    csi = node.coords[0];
    eta = node.coords[1];

    /* 
      Jm = [dx_dcsi  dx_deta;
            dy_dcsi  dy_deta]
      
      Ji = [dcsi_dx  dcsi_dy;
            deta_dx  deta_dy]
    */
    double dx_dcsi, dx_deta, dy_dcsi, dy_deta;
    double dcsi_dx, dcsi_dy, deta_dx, deta_dy;

    dx_dcsi = (a1 + c1 * eta);
    dx_deta = (b1 + c1 * csi);
    dy_dcsi = (a2 + c2 * eta);
    dy_deta = (b2 + c2 * csi);

    this->J[index][s_index - 1] = dx_dcsi * dy_deta - dx_deta * dy_dcsi;
    this->Jm[index][s_index - 1] = {dx_dcsi, dx_deta,
                                    dy_dcsi, dy_deta};

    dcsi_dx = +dy_deta / this->J[index][s_index - 1];
    dcsi_dy = -dx_deta / this->J[index][s_index - 1];
    deta_dx = -dy_dcsi / this->J[index][s_index - 1];
    deta_dy = +dx_dcsi / this->J[index][s_index - 1];

    this->Ji[index][s_index - 1] = {dcsi_dx, dcsi_dy,
                                    deta_dx, deta_dy};
  }

  index = 0;
  for (auto &vec_lines : fnodes)
  {
    index++;

    f_index = 0;
    for (auto &node : vec_lines) // line nodes for a specific direction (x, y)
    {
      f_index++;
      // flux node coordinates
      csi = node.coords[0];
      eta = node.coords[1];

      /* 
        Jm = [dx_dcsi  dx_deta;
              dy_dcsi  dy_deta]
        
        Ji = [dcsi_dx  dcsi_dy;
              deta_dx  deta_dy]
      */
      double dx_dcsi, dx_deta, dy_dcsi, dy_deta;
      double dcsi_dx, dcsi_dy, deta_dx, deta_dy;

      dx_dcsi = (a1 + c1 * eta);
      dx_deta = (b1 + c1 * csi);
      dy_dcsi = (a2 + c2 * eta);
      dy_deta = (b2 + c2 * csi);

      this->J[index][f_index - 1] = dx_dcsi * dy_deta - dx_deta * dy_dcsi;
      this->Jm[index][f_index - 1] = {dx_dcsi, dx_deta,
                                      dy_dcsi, dy_deta};

      dcsi_dx = +dy_deta / this->J[index][f_index - 1];
      dcsi_dy = -dx_deta / this->J[index][f_index - 1];
      deta_dx = -dy_dcsi / this->J[index][f_index - 1];
      deta_dy = +dx_dcsi / this->J[index][f_index - 1];

      this->Ji[index][f_index - 1] = {dcsi_dx, dcsi_dy,
                                      deta_dx, deta_dy};
    }
  }
}

double Quadrangle::calculate_jacobian_at_node(const Node& node)
{
  auto a1 = this->metrics[0];
  auto b1 = this->metrics[1];
  auto c1 = this->metrics[2];
  auto d1 = this->metrics[3];
  auto a2 = this->metrics[4];
  auto b2 = this->metrics[5];
  auto c2 = this->metrics[6];
  auto d2 = this->metrics[7];

  auto csi = node.coords[0];
  auto eta = node.coords[1];

  auto dx_dcsi = (a1 + c1 * eta);
  auto dx_deta = (b1 + c1 * csi);
  auto dy_dcsi = (a2 + c2 * eta);
  auto dy_deta = (b2 + c2 * csi);

  return dx_dcsi * dy_deta - dx_deta * dy_dcsi;
}

Node Quadrangle::transform(const Node &n, const std::vector<Vertice> &enodes)
{
  /* 
    Apply mapping from csi, eta (computational) -> x, y (physical)
   */
  double x = 0.0, y = 0.0, z = 0.0;
  double csi = 0.0, eta = 0.0;

  // Physical to Computational
  csi = n.coords[0];
  eta = n.coords[1];

  x = 0.25 * (this->metrics[0] * csi + this->metrics[1] * eta + this->metrics[2] * csi * eta + this->metrics[3]);
  y = 0.25 * (this->metrics[4] * csi + this->metrics[5] * eta + this->metrics[6] * csi * eta + this->metrics[7]);
  return Node{x, y, z};
}

/* ------------------------------ */
/*     QUADRANGLE 25-Nodes        */
/* ------------------------------ */
std::vector<long> QuadrangleHO::get_nodes_by_local_edge_id(long edge_id, bool sorted)
{
  // This functions returns a vector of vertices ids
  // of a required local edge.
  // First it gets the vertice at the edge_id position
  // Then validates if it's
  // Example:
  //   cell.nodes = {5, 12, 4, 0, ....} // only the first NUM_VERTICES are vertices

  std::vector<long> edge_vertices;
  edge_vertices.push_back(this->nodes[edge_id]);
  edge_vertices.push_back(this->nodes[(edge_id == NUM_VERTICES - 1) ? 0 : edge_id + 1]);
  if (sorted)
  {
    std::sort(edge_vertices.begin(), edge_vertices.end());
  }
  return edge_vertices;
}

void QuadrangleHO::enumerate_edges(void)
{
  // 2D Quad with N nodes
  // vertices: 0 1 2 3 4 ... N
  // edges:
  //		  (0) 0-1
  //		  (1) 1-2
  //		  (2) 2-3
  //		  (3) 3-0
  //      (4) ...
  //std::cout << "Enumerating edges of a QuadrangleHO" << std::endl;
  this->NUM_NODES = this->nodes.size();
  this->ORDER = static_cast<int>(sqrt(this->NUM_NODES)) - 1;

  // Calculates computational map that creates a unordered_map (computational_map)
  // whose keys are the cell's node and its value gives a vector
  // for the i, j coordinates in the isoparametric space
  this->calculate_computational_map();

  long local_id = 0;
  while (local_id < this->NUM_FACES)
  {
    std::vector<long> ordered_vertices(
        this->get_nodes_by_local_edge_id(local_id, true)); // true returns sorted vector
    std::vector<long> unordered_vertices(
        this->get_nodes_by_local_edge_id(local_id, false)); // false returns unsorted vector

    if (!was_enumerated(ordered_vertices))
    {
      this->num_edges++;
      long edge_id = this->num_edges;
      // if edge does not exist yet
      edges_map.emplace(
          std::make_pair(ordered_vertices, std::vector<long>{edge_id, this->id, -1}));

      this->edges.emplace_back(Edge{unordered_vertices, edge_id, this->id, -1});
    }
    else // if the edge exist, find its id
    {
      long edge_id = edges_map.find(ordered_vertices)->second[0];
      long neighbor_id = edges_map.find(ordered_vertices)->second[1];

      this->edges.emplace_back(Edge{unordered_vertices, edge_id, this->id, neighbor_id});
    }
    ordered_vertices.clear();
    unordered_vertices.clear();
    local_id++;
  }
  //edges_map.clear();
}

void QuadrangleHO::allocate_jacobian(int order)
{

  // SPs = order*order
  // FPs = (order+1)*order
  // Jm[0][0:SPs][0:3]
  // Jm[1][0:FPs][0:3] x
  // Jm[2][0:FPs][0:3] y

  // Ji[0][0:SPs][0:3]
  // Ji[1][0:FPs][0:3] x
  // Ji[2][0:FPs][0:3] y

  this->J = DoubleArr2D(3, std::vector<double>(order * (order + 1), 0.0));
  this->Jm = DoubleArr3D(3, DoubleArr2D((order + 1) * order, std::vector<double>(4, 0.0)));
  this->Ji = DoubleArr3D(3, DoubleArr2D((order + 1) * order, std::vector<double>(4, 0.0)));
  //this->Jm = DoubleArr3D{3, DoubleArr2D{(order + 1) * order, std::vector<double>{4}}};
  //this->Ji = DoubleArr3D{3, DoubleArr2D{(order + 1) * order, std::vector<double>{4}}};
}

void QuadrangleHO::calculate_jacobian(const std::vector<Node> &snodes,
                                      const std::vector<std::vector<Node>> &fnodes,
                                      const std::vector<Vertice> &enodes)
{
  std::vector<double> x, y;
  x.reserve(this->NUM_NODES);
  y.reserve(this->NUM_NODES);

  this->ce_space.reserve(sqrt(this->NUM_NODES));
  double ds = 2.0 / double(this->ORDER);
  double Lcsi = 0.0, Leta = 0.0, dLcsi = 0.0, dLeta = 0.0;

  for (size_t i = 0; i <= this->ORDER; i++)
    this->ce_space.push_back(-1.0 + ds * i);

  // Shape functions in 2D will be constructed based on
  // the Lagrange Polynomials applied to the isoparametric space
  // due to all nodes are allocated in a specific coordinate
  // Once these nodes' locations are known, it can be used to interpolate
  // each dimension separetely and multiply them in order to obtain
  // the set of shape functions that map a node from isoparametric space
  // to the physical space. For the inverse mapping, it will be used the
  // inverse of the Jacobian matrix.
  Helpers<Lagrange>::init();
  Helpers<Lagrange>::set_nodes(this->ce_space);

  double csi, eta;
  std::size_t index = 0, s_index = 0, f_index = 0;
  double dx_dcsi = 0.0, dx_deta = 0.0, dy_dcsi = 0.0, dy_deta = 0.0;
  double dcsi_dx = 0.0, dcsi_dy = 0.0, deta_dx = 0.0, deta_dy = 0.0;

  for (auto &node : snodes)
  {
    s_index++;
    // flux node coordinates
    csi = node.coords[0];
    eta = node.coords[1];

    /* 
      Jm = [dx_dcsi  dx_deta;
            dy_dcsi  dy_deta]
      
      Ji = [dcsi_dx  dcsi_dy;
            deta_dx  deta_dy]
    */

    for (int k = 0; k < this->NUM_NODES; k++)
    {
      auto coords = this->computational_map.find(k)->second;
      int i = coords[0];
      int j = coords[1];

      Lcsi = Helpers<Lagrange>::Pn(i, csi);
      Leta = Helpers<Lagrange>::Pn(j, eta);
      dLcsi = Helpers<Lagrange>::dPn(i, csi);
      dLeta = Helpers<Lagrange>::dPn(j, eta);

      dx_dcsi += dLcsi * Leta * enodes[this->nodes[k]].coords[0];
      dx_deta += Lcsi * dLeta * enodes[this->nodes[k]].coords[0];
      dy_dcsi += dLcsi * Leta * enodes[this->nodes[k]].coords[1];
      dy_deta += Lcsi * dLeta * enodes[this->nodes[k]].coords[1];
    }
    this->J[index][s_index - 1] = dx_dcsi * dy_deta - dx_deta * dy_dcsi;
    this->Jm[index][s_index - 1] = {dx_dcsi, dx_deta,
                                    dy_dcsi, dy_deta};

    dcsi_dx = +dy_deta / this->J[index][s_index - 1];
    dcsi_dy = -dx_deta / this->J[index][s_index - 1];
    deta_dx = -dy_dcsi / this->J[index][s_index - 1];
    deta_dy = +dx_dcsi / this->J[index][s_index - 1];

    this->Ji[index][s_index - 1] = {dcsi_dx, dcsi_dy,
                                    deta_dx, deta_dy};
  }

  index = 0;
  for (auto &vec_lines : fnodes)
  {
    index++;

    f_index = 0;
    for (auto &node : vec_lines) // line nodes for a specific direction (x, y)
    {
      f_index++;
      // flux node coordinates
      csi = node.coords[0];
      eta = node.coords[1];

      for (int k = 0; k < this->NUM_NODES; k++)
      {
        auto coords = this->computational_map.find(k)->second;
        int i = coords[0];
        int j = coords[1];

        Lcsi = Helpers<Lagrange>::Pn(i, csi);
        Leta = Helpers<Lagrange>::Pn(j, eta);
        dLcsi = Helpers<Lagrange>::dPn(i, csi);
        dLeta = Helpers<Lagrange>::dPn(j, eta);

        dx_dcsi += dLcsi * Leta * enodes[this->nodes[k]].coords[0];
        dx_deta += Lcsi * dLeta * enodes[this->nodes[k]].coords[0];
        dy_dcsi += dLcsi * Leta * enodes[this->nodes[k]].coords[1];
        dy_deta += Lcsi * dLeta * enodes[this->nodes[k]].coords[1];
      }
      /* 
        Jm = [dx_dcsi  dx_deta;
              dy_dcsi  dy_deta]
        
        Ji = [dcsi_dx  dcsi_dy;
              deta_dx  deta_dy]
      */
      this->J[index][f_index - 1] = dx_dcsi * dy_deta - dx_deta * dy_dcsi;
      this->Jm[index][f_index - 1] = {dx_dcsi, dx_deta,
                                      dy_dcsi, dy_deta};

      dcsi_dx = +dy_deta / this->J[index][f_index - 1];
      dcsi_dy = -dx_deta / this->J[index][f_index - 1];
      deta_dx = -dy_dcsi / this->J[index][f_index - 1];
      deta_dy = +dx_dcsi / this->J[index][f_index - 1];

      this->Ji[index][f_index - 1] = {dcsi_dx, dcsi_dy,
                                      deta_dx, deta_dy};
    }
  }
}

void QuadrangleHO::recursive_computational_map(std::vector<int> indices, int i0, int j0)
{
  double N = indices.size();
  double v;
  int order = static_cast<int>(sqrt(N)) - 1;

  if (N == 1)
  {
    v = indices.back();
    indices.pop_back();
    this->computational_map.insert({v, {i0, j0}});
  }

  if (N != 0)
  {
    std::vector<long> vertices;
    int i, j;
    // Nodes at the vertices
    while (vertices.size() < 4)
    {
      vertices.push_back(indices.back());
      indices.pop_back();
      int idx = 0;
      for (auto v : vertices)
      {
        if (idx == 0)
        {
          i = i0;
          j = j0;
        }
        if (idx == 1)
          i += order;
        if (idx == 2)
          j += order;
        if (idx == 3)
          j -= order;

        this->computational_map.insert({v, {i, j}});
      }
    }

    if (indices.size() > 0)
    {
      // Nodes at edges
      std::vector<int> edges;
      while (edges.size() < 4 * (order - 1))
      {
        edges.push_back(indices.back());
        indices.pop_back();
      }

      i = i0;
      j = j0;
      double L = order - 1;
      int idx = 0;
      for (auto e : edges)
      {
        if (idx < L)
          i += 1;
        else if (idx < 2 * L)
          j += 1;
        else if (idx < 3 * L)
          i -= 1;
        else if (idx < 4 * L)
          j -= 1;

        this->computational_map.insert({e, {i, j}});

        if (idx == L - 1)
          i += 1;
        if (idx == 2 * L - 1)
          j += 1;
        if (idx == 3 * L - 1)
          i -= 1;
        if (idx == 4 * L - 1)
          i += 1;
      }
    }
    this->recursive_computational_map(indices, i, j);
  }
}

void QuadrangleHO::calculate_computational_map(void)
{
  double ds = 2.0 / double(this->ORDER);
  std::vector<int> indices;
  indices.reserve((this->ORDER + 1) * (this->ORDER + 1));
  for (auto i = (this->ORDER + 1) * (this->ORDER + 1) - 1; i < 0; i--)
    indices.push_back(i);

  this->recursive_computational_map(indices, 0, 0);
}

Node QuadrangleHO::transform(const Node &n, const std::vector<Vertice> &enodes)
{
  /* 
    Apply mapping from csi, eta (computational) -> x, y (physical)
   */
  double x = 0.0, y = 0.0, z = 0.0;
  double csi = 0.0, eta = 0.0;
  double Lcsi = 0.0, Leta = 0.0;

  Helpers<Lagrange>::init();
  Helpers<Lagrange>::set_nodes(this->ce_space);

  csi = n.coords[0];
  eta = n.coords[1];

  // Computational to Physical
  for (int k = 0; k < this->NUM_NODES; k++)
  {
    auto coords = this->computational_map.find(k)->second;
    int i = coords[0];
    int j = coords[1];

    Lcsi = Helpers<Lagrange>::Pn(i, csi);
    Leta = Helpers<Lagrange>::Pn(j, eta);

    x += Lcsi * Leta * enodes[this->nodes[k]].coords[0];
    y += Lcsi * Leta * enodes[this->nodes[k]].coords[1];
  }

  return Node{x, y, z};
  //return Node();
}


void Triangle::allocate_jacobian(int order)
{
}

void Tetrahedron::allocate_jacobian(int order)
{
}

void Hexahedron::allocate_jacobian(int order)
{
}

void Prism::allocate_jacobian(int order)
{
}

void Pyramid::allocate_jacobian(int order)
{
}

void Triangle::calculate_jacobian(const std::vector<Node> &snodes,
                                  const std::vector<std::vector<Node>> &fnodes,
                                  const std::vector<Vertice> &enodes)
{
}

void Tetrahedron::calculate_jacobian(const std::vector<Node> &snodes,
                                     const std::vector<std::vector<Node>> &fnodes,
                                     const std::vector<Vertice> &enodes)
{
}

void Hexahedron::calculate_jacobian(const std::vector<Node> &snodes,
                                    const std::vector<std::vector<Node>> &fnodes,
                                    const std::vector<Vertice> &enodes)
{
}

void Prism::calculate_jacobian(const std::vector<Node> &snodes,
                               const std::vector<std::vector<Node>> &fnodes,
                               const std::vector<Vertice> &enodes)
{
}

void Pyramid::calculate_jacobian(const std::vector<Node> &snodes,
                                 const std::vector<std::vector<Node>> &fnodes,
                                 const std::vector<Vertice> &enodes)
{
}

Node Triangle::transform(const Node &n, const std::vector<Vertice> &enodes)
{
}

Node Tetrahedron::transform(const Node &n, const std::vector<Vertice> &enodes)
{
}

Node Hexahedron::transform(const Node &n, const std::vector<Vertice> &enodes)
{
}

Node Prism::transform(const Node &n, const std::vector<Vertice> &enodes)
{
}

Node Pyramid::transform(const Node &n, const std::vector<Vertice> &enodes)
{
}

double QuadrangleHO::calculate_jacobian_at_node(const Node& node)
{
  return 0.0;
}

double Tetrahedron::calculate_jacobian_at_node(const Node& node)
{
  return 0.0;
}

double Hexahedron::calculate_jacobian_at_node(const Node& node)
{
  return 0.0;
}

double Prism::calculate_jacobian_at_node(const Node& node)
{
  return 0.0;
}

double Pyramid::calculate_jacobian_at_node(const Node& node)
{
  return 0.0;
}

void Tetrahedron::enumerate_faces(void)
{
  // 3D Tetrahedron with 4 nodes
  // vertices: 0 1 2 3
  // faces:
  //		  (0) 0-1-2
  //		  (1) 0-3-1
  //		  (2) 0-3-2
  std::cout << "I am a Tetrahedron" << std::endl;
}

void Hexahedron::enumerate_faces(void)
{
  // 3D Hexahedron with 8 nodes
  // vertices: 0 1 2 3 4 5 6 7
  // faces:
  //		  (0) 0-3-2-1
  //		  (1) 0-1-5-4
  //		  (2) 0-4-7-3
  //		  (3) 1-2-6-5
  //		  (4) 1-2-7-6
  //		  (5) 4-5-6-7
  std::cout << "I am a Hexahedron" << std::endl;
}

void Prism::enumerate_faces(void)
{
  // 3D Prism with 6 nodes
  // vertices: 0 1 2 3 4 5
  // faces:
  //		  (0) 0-2-1
  //		  (1) 0-1-4-3
  //		  (2) 0-3-5-2
  //		  (3) 1-2-5-4
  //		  (4) 3-4-5
  std::cout << "I am a Prism" << std::endl;
}

void Pyramid::enumerate_faces(void)
{
  // 3D Pyramid with 5 nodes
  // vertices: 0 1 2 3 4
  // faces:
  //		  (0) 0-3-2-1
  //		  (1) 0-1-4
  //		  (2) 0-4-3
  //		  (3) 1-2-4
  //		  (4) 2-3-4
  std::cout << "I am a Pyramid" << std::endl;
}

// Print_vertices
void Triangle::print_vertices(void)
{
  std::cout << "I am a Triangle" << std::endl;
}

void Quadrangle::print_vertices(void)
{
  std::cout << "I am a Quadrangle" << std::endl;
}

void QuadrangleHO::print_vertices(void)
{
  std::cout << "I am a High-Order Quadrangle" << std::endl;
}

void Tetrahedron::print_vertices(void)
{
  std::cout << "I am a Tetrahedron" << std::endl;
}

void Hexahedron::print_vertices(void)
{
  std::cout << "I am a Hexahedron" << std::endl;
}

void Prism::print_vertices(void)
{
  std::cout << "I am a Prism" << std::endl;
}

void Pyramid::print_vertices(void)
{
  std::cout << "I am a Pyramid" << std::endl;
}

// get_vertices
void Triangle::get_vertices(void)
{
  std::cout << "I am a Triangle" << std::endl;
}

void Quadrangle::get_vertices(void)
{
  std::cout << "I am a Quadrangle" << std::endl;
}

void QuadrangleHO::get_vertices(void)
{
  std::cout << "I am a  High-Order Quadrangle" << std::endl;
}

void Tetrahedron::get_vertices(void)
{
  std::cout << "I am a Tetrahedron" << std::endl;
}

void Hexahedron::get_vertices(void)
{
  std::cout << "I am a Hexahedron" << std::endl;
}

void Prism::get_vertices(void)
{
  std::cout << "I am a Prism" << std::endl;
}

void Pyramid::get_vertices(void)
{
  std::cout << "I am a Pyramid" << std::endl;
}
// ------------------------------------------------------------------------ //
