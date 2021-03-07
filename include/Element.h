#pragma once

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

#include <DVector.h>
#include <Edge.h>
#include <Node.h>
#include <Property.h>

/* ELEMENT INTERFACE CLASS */
class Element
{
public:
  int NUM_FACES = 0;
  size_t NUM_VERTICES = 0;
  size_t NUM_NODES = 0;

  //std::vector<int> nodes;
  std::vector<long> nodes;
  std::vector<Edge> edges;
  std::vector<Edge> faces;
  std::vector<double> metrics;

  std::shared_ptr<Property> physical;
  std::shared_ptr<Property> computational;

  static std::unordered_map<std::vector<long>, std::vector<long>> faces_map;
  static std::unordered_map<std::vector<long>, std::vector<long>> edges_map;
  long id;
  int fringe;
  int boundary;
  //double J; // Jacobian

  // Jm[0][0:SPs][0:3/8]
  // Jm[1][0:FPs][0:3/8] x
  // Jm[2][0:FPs][0:3/8] y
  // Jm[3][0:FPs][0:3/8] z

  // Ji[0][0:SPs][0:3/8]
  // Ji[1][0:FPs][0:3/8] x
  // Ji[2][0:FPs][0:3/8] y
  // Ji[3][0:FPs][0:3/8] z
  std::vector<std::vector<double>> J;               // Jacobian
  std::vector<std::vector<std::vector<double>>> Jm; // Jacobian Matrix
  std::vector<std::vector<std::vector<double>>> Ji; // Jacobian Inverse Matrix

  static long num_elems;
  static long num_faces;
  static long num_edges;

public:
  Element(const std::vector<std::string> &);
  virtual ~Element();

  void print_nodes(void);
  void get_nodes(void);
  bool was_enumerated(const std::vector<long> &);
  virtual Node transform(const Node &, const std::vector<Vertice> &) = 0;
  virtual void print_vertices(void) = 0;
  virtual void get_vertices(void) = 0;
  virtual void allocate_jacobian(int) = 0;
  virtual void calculate_jacobian(const std::vector<Node> &,
                                  const std::vector<std::vector<Node>> &,
                                  const std::vector<Vertice> &) = 0;
  virtual double calculate_jacobian_at_node(const Node&, const std::vector<Vertice>&) = 0;
  virtual std::vector<double> get_normal_vector(int, int, int) = 0;
  /*
  get_normal_vector: 
    calculate normal vector components in physical space for a specific flux point
        Parameters:
          index:    1 - csi, 2 - eta
          point:    0-(order*(order+1))
          local_ed: 0-4
   */


private:
};

/* DERIVED ELEMENT CLASSES */
// Triangle
class Triangle : public Element
{
public:
  Triangle(const std::vector<std::string> &node_list) : Element(node_list)
  {
    this->NUM_FACES = 3;
    this->NUM_VERTICES = 3;
    this->NUM_NODES = 3;
    this->enumerate_edges();
  }
  ~Triangle(void)
  { /*std::cout << "Triangle has been deleted!" << std::endl; */
  }

  Node transform(const Node &, const std::vector<Vertice> &);
  void print_vertices(void);
  void get_vertices(void);
  void enumerate_edges(void);
  void allocate_jacobian(int);
  void calculate_jacobian(const std::vector<Node> &,
                          const std::vector<std::vector<Node>> &,
                          const std::vector<Vertice> &);
  double calculate_jacobian_at_node(const Node&, const std::vector<Vertice>&);
  std::vector<long> get_ordered_nodes_by_local_edge_id(long);
  std::vector<double> get_normal_vector(int, int, int);
};

// Quadrangle
class Quadrangle : public Element
{
public:
  Quadrangle(const std::vector<std::string> &node_list) : Element(node_list)
  {
    this->NUM_FACES = 4;
    this->NUM_VERTICES = 4;
    this->NUM_NODES = 4;
    this->enumerate_edges();
  }
  ~Quadrangle(void)
  { /*std::cout << "Quadrangle has been deleted!" << std::endl; */
  }

  Node transform(const Node &, const std::vector<Vertice> &);
  void print_vertices(void);
  void get_vertices(void);
  void enumerate_edges(void);
  void allocate_jacobian(int);
  void calculate_jacobian(const std::vector<Node> &,
                          const std::vector<std::vector<Node>> &,
                          const std::vector<Vertice> &);
  double calculate_jacobian_at_node(const Node&, const std::vector<Vertice>&);
  std::vector<long> get_nodes_by_local_edge_id(long, bool);
  std::vector<double> get_normal_vector(int, int, int);
};

// Quadrangle High-Order
class QuadrangleHO : public Element
{
  
  int ORDER;
  std::vector<double> ce_space;
  std::unordered_map<long, std::vector<long>> computational_map;

public:
  QuadrangleHO(const std::vector<std::string> &node_list) : Element(node_list)
  {
    this->NUM_FACES = 4;
    this->NUM_VERTICES = 4;
    this->NUM_NODES = 0;
    this->enumerate_edges();
  }
  ~QuadrangleHO(void)
  { /*std::cout << "QuadrangleHO has been deleted!" << std::endl; */
  }

  Node transform(const Node &, const std::vector<Vertice> &);
  void print_vertices(void);
  void get_vertices(void);
  void enumerate_edges(void);
  void allocate_jacobian(int);
  void calculate_jacobian(const std::vector<Node> &,
                          const std::vector<std::vector<Node>> &,
                          const std::vector<Vertice> &);
  double calculate_jacobian_at_node(const Node&, const std::vector<Vertice>&);
  std::vector<long> get_nodes_by_local_edge_id(long, bool);
  std::vector<double> get_normal_vector(int, int, int);

private:
  void calculate_computational_map(void);
  void recursive_computational_map(std::vector<long>, long, long);
};

// Tetrahedron
class Tetrahedron : public Element
{
  
public:
  Tetrahedron(const std::vector<std::string> &node_list) : Element(node_list) {
    this->NUM_FACES = 4;
    this->NUM_VERTICES = 4;
    this->NUM_NODES = 4;
  }
  ~Tetrahedron(void) { std::cout << "Tetrahedron has been deleted!" << std::endl; }

  Node transform(const Node &, const std::vector<Vertice> &);
  void print_vertices(void);
  void get_vertices(void);
  void enumerate_faces(void);
  void allocate_jacobian(int);
  void calculate_jacobian(const std::vector<Node> &,
                          const std::vector<std::vector<Node>> &,
                          const std::vector<Vertice> &);
  double calculate_jacobian_at_node(const Node&, const std::vector<Vertice>&);
  std::vector<double> get_normal_vector(int, int, int);
  // std::vector<long> get_nodes_by_local_edge_id(long, bool);
};

// Hexahedron
class Hexahedron : public Element
{
  
public:
  Hexahedron(const std::vector<std::string> &node_list) : Element(node_list) {
    this->NUM_FACES = 6;
    this->NUM_VERTICES = 8;
    this->NUM_NODES = 8;
  }
  ~Hexahedron(void) { std::cout << "Hexahedron has been deleted!" << std::endl; }

  Node transform(const Node &, const std::vector<Vertice> &);
  void print_vertices(void);
  void get_vertices(void);
  void enumerate_faces(void);
  void allocate_jacobian(int);
  void calculate_jacobian(const std::vector<Node> &,
                          const std::vector<std::vector<Node>> &,
                          const std::vector<Vertice> &);
  double calculate_jacobian_at_node(const Node&, const std::vector<Vertice>&);
  std::vector<double> get_normal_vector(int, int, int);
  // std::vector<long> get_nodes_by_local_edge_id(long, bool);
};

// Prism
class Prism : public Element
{
  
public:
  Prism(const std::vector<std::string> &node_list) : Element(node_list) {
    this->NUM_FACES = 5;
    this->NUM_VERTICES = 6;
    this->NUM_NODES = 6;
  }
  ~Prism(void) { std::cout << "Prism has been deleted!" << std::endl; }

  Node transform(const Node &, const std::vector<Vertice> &);
  void print_vertices(void);
  void get_vertices(void);
  void enumerate_faces(void);
  void allocate_jacobian(int);
  void calculate_jacobian(const std::vector<Node> &,
                          const std::vector<std::vector<Node>> &,
                          const std::vector<Vertice> &);
  double calculate_jacobian_at_node(const Node&, const std::vector<Vertice>&);
  std::vector<double> get_normal_vector(int, int, int);
  // std::vector<long> get_nodes_by_local_edge_id(long, bool);
};

// Pyramid
class Pyramid : public Element
{
public:
  Pyramid(const std::vector<std::string> &node_list) : Element(node_list) {
    this->NUM_FACES = 5;
    this->NUM_VERTICES = 5;
    this->NUM_NODES = 5;
  }
  ~Pyramid(void) { std::cout << "Pyramid has been deleted!" << std::endl; }

  Node transform(const Node &, const std::vector<Vertice> &);
  void print_vertices(void);
  void get_vertices(void);
  void enumerate_faces(void);
  void allocate_jacobian(int);
  void calculate_jacobian(const std::vector<Node> &,
                          const std::vector<std::vector<Node>> &,
                          const std::vector<Vertice> &);
  double calculate_jacobian_at_node(const Node&, const std::vector<Vertice>&);
  std::vector<double> get_normal_vector(int, int, int);
  // std::vector<long> get_nodes_by_local_edge_id(long, bool);
};

// https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_5_6/Common/GmshDefines.h
enum class elm_type
{
  NODE2_LINE = 1,
  NODE3_TRI,
  NODE4_QUAD,
  NODE4_TETRA,
  NODE8_HEXA,
  NODE6_PRIS,
  NODE5_PYR,
  NODE3_O2_LINE,
  NODE6_O2_TRI,
  NODE9_O2_QUAD,
  NODE10_O2_TETRA,
  NODE27_O2_HEXA,
  NODE18_O2_PRIS,
  NODE14_O2_PYR,
  NODE1_POINT,
  NODE8_O2_QUAD,
  NODE20_O2_HEXA,
  NODE15_O2_PRIS,
  NODE13_O2_PYR,
  NODE9_O3_ITRI,
  NODE10_O3_TRI,
  NODE12_O4_ITRI,
  NODE15_O4_TRI,
  NODE15_O5_ITRI,
  NODE21_O5_TRI,
  NODE4_O3_EDG,
  NODE5_O4_EDG,
  NODE6_O5_EDG,
  NODE20_O3_TETRA,
  NODE35_O4_TETRA,
  NODE56_O5_TETRA,
  NODE22_TETRA,
  NODE28_TETRA,
  NODE_POLYG,
  NODE_POLYH,
  NODE16_O3_QUAD,
  NODE25_O4_QUAD,
  NODE36_O5_QUAD,
  NODE12_O3_QUAD,
  NODE16I_O3_QUAD,
  NODE20_O4_QUAD,
  NODE_TRI_28,
  NODE_TRI_36,
  NODE_TRI_45,
  NODE_TRI_55,
  NODE_TRI_66,
  NODE49_O6_QUAD,
  NODE64_O7_QUAD,
  NODE81_O8_QUAD,
  NODE100_O9_QUAD,
  NODE121_O10_QUAD,
};