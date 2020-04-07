#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <iterator> 
#include <algorithm>  

#include <Element.h>
#include <Mesh.h>
#include <Node.h>

// Nodes
Node::Node(const std::vector<std::string>& coordinates)
{
  num_nodes++;
  this->id = num_nodes;

  // transforming vector<string> to vector<double> using lambda
  std::transform(coordinates.begin(), coordinates.end(), std::back_inserter((this->coords)),
    [](const std::string& str) { return std::stod(str); });

  this->left = -1;
  this->right = -1;

  //std::cout << "Node(" << this->id << ") has been created with coords: ";
  //for (const auto& v : (this->coords))
  //	std::cout << v << " ";
  //std::cout << std::endl;

  //std::cout << "Node has been created!" << std::endl;
}

Node::Node(double x, double y, double z)
{
  this->id = -1;
  this->coords.resize(3);

  this->coords[0] = x;
  this->coords[1] = y;
  this->coords[2] = z;

  this->left = -1;
  this->right = -1;

  //std::cout << "Test Node(" << this->id << ") has been created with coords: ";
  //for (const auto& v : (this->coords))
  //	std::cout << v << " ";
  //std::cout << std::endl;

  //std::cout << "Test Node has been created!" << std::endl;
}

Node::Node()
{
  this->id = -1;
  this->coords.resize(3);

  this->coords[0] = 0.0;
  this->coords[1] = 0.0;
  this->coords[2] = 0.0;

  this->left = -1;
  this->right = -1;

  //std::cout << "Test Node(" << this->id << ") has been created with coords: ";
  //for (const auto& v : (this->coords))
  //	std::cout << v << " ";
  //std::cout << std::endl;

  //std::cout << "Test Node has been created!" << std::endl;
}

Node::~Node(void)
{
  //std::cin.get();
  //std::cout << "Node(" << this->id << ") has been deleted!" << std::endl;
  //std::cin.get();
  
}

void Node::print_coords(void)
{
  std::cout << "Node(" << this->id << "): ";
  for (const auto& v : (this->coords))
    std::cout << v << " ";
  std::cout << std::endl;
}

void Node::print_elements(void)
{
  std::cout << "Node(" << this->id << ").elements: ";
  for (const auto& e : (this->elems))
    std::cout << e << " ";
  std::cout << std::endl;
}

// ------------------------------------------------------------------------ //

// Edge
Edge::Edge(const std::vector<long>& e_nodes, const long& edge_id, const long& left, const long& right)
{
  //num_edges++;
  this->id = edge_id;
  this->boundary = 0;

  // copying edge[edge_id] vector of nodes to the new edge class
  this->nodes = e_nodes;

  // neighbors
  this->left = left;
  this->right = right;
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
  for (const auto& v : (this->nodes))
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
    [](const std::string& str) { return std::stol(str); });


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
  for (const auto& v : (*this->nodes))
    std::cout << v << " ";
  std::cout << std::endl;
}
// ------------------------------------------------------------------------ //



/* CONSTRUCTORS */
Element::Element(const std::vector<std::string>& node_list)
{
  num_elems++;
  this->id = num_elems;
  this->boundary = 0;
  this->fringe = 0;
  this->J = 0.0; // Jacobian

  //std::cout << "I am the (" << this->id << ") Element" << std::endl;

  // transforming vector<string> to vector<int>
  // subtract 1 from the node_id because the first id for GMSH is 1, but as I am using the node id to consult it on a vector
  // I need it to be 0-started.
  std::transform(node_list.begin(), node_list.end(), std::back_inserter(this->nodes),
    [](const std::string& str) { return std::stol(str)-1; }); 

  // this->print_nodes();
}

/* DESTRUCTORS */
Element::~Element(void)
{
  //std::cout << "Element (" << this->id << ") has been deleted!" << std::endl;
  //std::cin.get();
}


/* ELEMENT METHODS */
void Element::print_nodes(void)
{
  std::cout << "Element(" << this->id << "): ";
  for (const auto& v : this->nodes)
    std::cout << v << " ";
  std::cout << std::endl;
}

void Element::get_nodes(void)
{

}



bool Element::was_enumerated(const std::vector<long>& vertices)
{
  auto edges_map_iterator = this->edges_map.find(vertices);
  if (edges_map_iterator == this->edges_map.end()) return false;
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
    std::vector<long> ordered_vertices(this->get_nodes_by_local_edge_id(local_id, true)); // true returns sorted vector
    std::vector<long> unordered_vertices(this->get_nodes_by_local_edge_id(local_id, false)); // false returns unsorted vector
    
    if (!was_enumerated(ordered_vertices))
    {
      this->num_edges++;
      long edge_id = this->num_edges;
      // if edge does not exist yet
      edges_map.emplace(std::make_pair(ordered_vertices, std::vector<long>{edge_id, this->id, -1}));
      
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

      this->edges.emplace_back(Edge{ unordered_vertices, edge_id, this->id, -1});
      
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
      this->edges.emplace_back(Edge{ unordered_vertices, edge_id, this->id, neighbor_id});

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
  
  this->Jm = std::vector<std::vector<std::vector<double>>>{3, std::vector<std::vector<double>>{(order+1)*order, std::vector<double>{4}}};
  this->Ji = std::vector<std::vector<std::vector<double>>>{3, std::vector<std::vector<double>>{(order+1)*order, std::vector<double>{4}}};
  
}

void Quadrangle::calculate_jacobian(const std::vector<Node>& snodes, 
                                    const std::vector<std::vector<Node>>& fnodes, 
                                    const std::vector<Node>& enodes)
{
  double x1 = enodes[this->nodes[0]].coords[0];
  double x2 = enodes[this->nodes[1]].coords[0];
  double x3 = enodes[this->nodes[2]].coords[0];
  double x4 = enodes[this->nodes[3]].coords[0];

  double y1 = enodes[this->nodes[0]].coords[1];
  double y2 = enodes[this->nodes[1]].coords[1];
  double y3 = enodes[this->nodes[2]].coords[1];
  double y4 = enodes[this->nodes[3]].coords[1];


  double a1, a2, b1, b2, c1, c2, csi, eta;

  this->J = 0.5*(x2-x1)*0.5*(y4-y3);
  
  a1 =  x2-x1+x3-x4;
  b1 = -x2-x1+x3+x4;
  c1 = -x2+x1+x3-x4;
  a2 =  y2-y1+y3-y4;
  b2 = -y2-y1+y3+y4;
  c2 = -y2+y1+y3-y4;

  std::size_t index=0, s_index=0, f_index=0;
  
  for (auto& node: snodes)
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
   
    this->Jm[index][s_index-1] = {0.25*(a1 + c1*eta), 0.25*(b1 + c1*csi), 
                                  0.25*(a2 + c2*eta), 0.25*(b2 + c2*csi)};
    this->Ji[index][s_index-1] = {+this->Jm[index][s_index-1][3], -this->Jm[index][s_index-1][1], 
                                  -this->Jm[index][s_index-1][2], +this->Jm[index][s_index-1][0]};
  }

  index=0;
  for (auto& vec_lines : fnodes)
  {
    index++;

    f_index = 0;
    for (auto& node : vec_lines) // line nodes for a specific direction (x, y)
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
      this->Jm[index][f_index-1] = {0.25*(a1 + c1*eta), 0.25*(b1 + c1*csi), 
                                    0.25*(a2 + c2*eta), 0.25*(b2 + c2*csi)};
      this->Ji[index][f_index-1] = {+this->Jm[index][f_index-1][3], -this->Jm[index][f_index-1][1], 
                                    -this->Jm[index][f_index-1][2], +this->Jm[index][f_index-1][0]};
    }
  }
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



void Triangle::calculate_jacobian(const std::vector<Node>& snodes, 
                  const std::vector<std::vector<Node>>& fnodes, 
                  const std::vector<Node>& enodes)
{
}

void Tetrahedron::calculate_jacobian(const std::vector<Node>& snodes, 
                   const std::vector<std::vector<Node>>& fnodes, 
                   const std::vector<Node>& enodes)
{
}

void Hexahedron::calculate_jacobian(const std::vector<Node>& snodes, 
                  const std::vector<std::vector<Node>>& fnodes, 
                  const std::vector<Node>& enodes)
{
}

void Prism::calculate_jacobian(const std::vector<Node>& snodes, 
                 const std::vector<std::vector<Node>>& fnodes, 
                 const std::vector<Node>& enodes)
{
}

void Pyramid::calculate_jacobian(const std::vector<Node>& snodes, 
                 const std::vector<std::vector<Node>>& fnodes, 
                 const std::vector<Node>& enodes)
{
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

