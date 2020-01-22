#pragma once

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <unordered_map> 
#include <functional>

#include <DVector.h>

// NODE CLASS
class Node
{
public:

	long id;
	static long num_nodes;
	std::vector<double> coords;
	long left;
	long right;
	int boundary;
	std::vector<long> elems;

public:
	Node(const std::vector<std::string>&);
	Node(double, double, double);
	~Node();

	void print_coords(void);
	void print_elements(void);
	
private:

};

// Edge CLASS
class Edge
{
public:
	std::vector<long> nodes;
	long left;
	long right;
	long id;
	int boundary;
public:
	Edge(const std::vector<long>&, const long&, const long&, const long&);
	~Edge();

	void print_nodes(void);

private:

};

// Face CLASS
class Face
{
public:

	std::shared_ptr<std::vector<long>> nodes = std::make_shared<std::vector<long>>();
	long id;

	
public:
	Face(std::vector<std::string>);
	~Face();

	void print_nodes(void);

private:

};


template<typename T>
struct std::hash<std::vector<T>>
{
	typedef std::vector<T> argument_type;
	typedef size_t result_type;

	result_type operator()(const argument_type& a) const
	{
		std::hash<T> hasher;
		result_type h = 0;
		for (result_type i = 0; i < a.size(); ++i)
		{
			h = h * 31 + hasher(a[i]);
		}
		return h;
	}
};

/* ELEMENT INTERFACE CLASS */
class Element
{
public:
	//std::vector<int> nodes;
	std::vector<long> nodes;
	std::vector<Edge> edges;
	std::vector<Edge> faces;
	
  /* 
    DVector represents a property vector containing all fluid 
    properties or flux properties;
    
    The need of a vector of DVector is due to the high-order 
    polynomial interpolation. So for each (solution/flux) node
    a DVector is allocated as an element inside this vector.

    Although, for flux properties as it has 2 (or 3) directions, 
    the first vector stands for DIRECTION (x, y, z),
    the second vector stands for the number of nodes 
    (as in solution properties).

    At last, for the residue, it has one DVector for each
    node similar to the solution properties.
	 */
   
  // Conservative Properties
  std::vector<DVector> Qsp;                 // Solution at Solution nodes
  std::vector<DVector> Qfp;                 // Solution at Flux nodes

  // Convective Fluxes
  std::vector<std::vector<DVector>> Fcsp;   // Convective Flux at Solution nodes
  std::vector<std::vector<DVector>> Fcfp;   // Convective Flux at Flux nodes 
  
  // Diffusive Fluxes
  std::vector<std::vector<DVector>> Fdsp;   // Diffusive Flux at Solution nodes
  std::vector<std::vector<DVector>> Fdfp;   // Diffusive Flux at Flux nodes
  
  // Gradients
  // Conservative Properties
  std::vector<std::vector<DVector>> dQsp;   // Solution Gradient at Solution nodes
  std::vector<std::vector<DVector>> dQfp;   // Solution Gradient at Flux nodes

  // Convective Fluxes
  std::vector<std::vector<DVector>> dFcsp;   // Convective Flux Gradient at Solution nodes
  std::vector<std::vector<DVector>> dFcfp;   // Convective Flux Gradient at Flux nodes 
  
  // Diffusive Fluxes
  std::vector<std::vector<DVector>> dFdsp;   // Diffusive Flux Gradient at Solution nodes
  std::vector<std::vector<DVector>> dFdfp;   // Diffusive Flux Gradient at Flux nodes
  
  // Residue Vector
  std::vector<DVector> res;

	static std::unordered_map<std::vector<long>, std::vector<long>> faces_map;
	static std::unordered_map<std::vector<long>, std::vector<long>> edges_map;
	long id;
	int fringe;
	int boundary;
  double J; // Jacobian
  std::vector<std::vector<double>> Jm; // Jacobian Matrix
  std::vector<std::vector<double>> Ji; // Jacobian Inverse Matrix


	static long num_elems;
	static long num_faces;
	static long num_edges;
public:
	Element(const std::vector<std::string>&);
	virtual ~Element();

	void print_nodes(void);
	void get_nodes(void);
	bool was_enumerated(const std::vector<long>&);
	virtual void print_vertices(void) = 0;
	virtual void get_vertices(void) = 0;
  virtual void allocate_jacobian(void) = 0;
  virtual void calculate_jacobian(const std::vector<Node>&) = 0;
private:

};

/* DERIVED ELEMENT CLASSES */
// Triangle
class Triangle : public Element
{
	const int NUM_FACES = 3;
public:
	Triangle(const std::vector<std::string>& node_list) : Element(node_list) { this->enumerate_edges(); this->allocate_jacobian();}
	~Triangle(void) { /*std::cout << "Triangle has been deleted!" << std::endl; */}

	void print_vertices(void);
	void get_vertices(void);
	void enumerate_edges(void);
	void allocate_jacobian(void);
  void calculate_jacobian(const std::vector<Node>&);
  std::vector<long> get_ordered_nodes_by_local_edge_id(long);
	
};

// Quadrangle
class Quadrangle : public Element
{
	const int NUM_FACES = 4;
public:
	Quadrangle(const std::vector<std::string>& node_list) : Element(node_list) { this->enumerate_edges(); this->allocate_jacobian();}
	~Quadrangle(void) { /*std::cout << "Quadrangle has been deleted!" << std::endl; */}

	void print_vertices(void);
	void get_vertices(void);
	void enumerate_edges(void);
  void allocate_jacobian(void);
  void calculate_jacobian(const std::vector<Node>&);
	std::vector<long> get_nodes_by_local_edge_id(long, bool);
};

// Tetrahedron
class Tetrahedron : public Element
{

public:
	Tetrahedron(const std::vector<std::string>& node_list) : Element(node_list) {}
	~Tetrahedron(void) { std::cout << "Tetrahedron has been deleted!" << std::endl; }

	void print_vertices(void);
	void get_vertices(void);
	void enumerate_faces(void);
};


// Hexahedron
class Hexahedron : public Element
{
public:
	Hexahedron(const std::vector<std::string>& node_list) : Element(node_list) {}
	~Hexahedron(void) { std::cout << "Hexahedron has been deleted!" << std::endl; }

	void print_vertices(void);
	void get_vertices(void);
	void enumerate_faces(void);

};

// Prism
class Prism : public Element
{
public:
	Prism(const std::vector<std::string>& node_list) : Element(node_list) {}
	~Prism(void) { std::cout << "Prism has been deleted!" << std::endl; }

	void print_vertices(void);
	void get_vertices(void);
	void enumerate_faces(void);

};

// Pyramid
class Pyramid : public Element
{

public:
	Pyramid(const std::vector<std::string>& node_list) : Element(node_list) {}
	~Pyramid(void) { std::cout << "Pyramid has been deleted!" << std::endl; }

	void print_vertices(void);
	void get_vertices(void);
	void enumerate_faces(void);
};


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
	NODE56_O5_TETRA

};


// MESH CLASS
class Mesh
{
protected:
	int id, dimension;
	long N, Nel;

public:
	static int number_meshes;
	std::vector<std::shared_ptr<Element>> elems;
	std::vector<Node> nodes;
public:
	Mesh(int);
	virtual ~Mesh();

	int getID();
	int getDimesion();
	long get_number_nodes();
	long get_number_elements();
	void print_element_by_id(long);
	void print_node_by_id(long);
	void print_node_id(long);
	void read_gmsh(const std::string&);
	double get_area(const std::vector<Node>&);
	double get_area(const std::vector<long>&, const Node&);
	double get_volume(const std::vector<Node>&);
	void update_element_neighbors(void);
	void mark_boundaries(void);
	void to_vtk(const std::string&);
  void calculate_jacobians(std::shared_ptr<Element>&);

private:
	void append_elem_to_nodes(const std::shared_ptr<Element>&);
};

// STATIC_MESH CLASS
class Static_Mesh : public Mesh
{
public:
	long root;
	static int dim;

public:
	Static_Mesh(int);
	~Static_Mesh();

	void createKDtree(void);	
	void build_kdtree(long&, std::vector<long>&, long);
	void get_pivot(long&, const std::vector<long>&);
	long mark_fringes(const Node&);
	void to_graphviz(void);
	void print_tree(long&, long&, std::string, int&);
	

private:
	//long _search(const Node&, long&, int&, long&);
	long _search(const Node&, long&, int, double, long);
	long _find_element(const Node&, std::vector<long>&);
	void _get_kneighbors(long&, std::vector<long>&);
	long _get_closest(const Node&, std::vector<long>&);
	double _get_distance(const Node&, const Node&);
};

// STATUS CLASS
class Status
{
public:
	Status() = delete;
	~Status();

	virtual void printAll() = 0;
private:

};



// DYNAMIC_MESH CLASS
class Dynamic_Mesh : public Mesh
{
public:

public:
	Dynamic_Mesh(int);
	~Dynamic_Mesh();

	void rotate(double);
	void translate(double, int);
	void move(int, double, double, int);
	void start(void);
private:
	
};
