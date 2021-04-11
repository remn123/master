#pragma once

#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include <functional>

#include <DVector.h>
#include <Element.h>
#include <Ghost.h>
#include <Node.h> // Vertice

template <typename T>
struct std::hash<std::vector<T>>
{
	typedef std::vector<T> argument_type;
	typedef size_t result_type;

	result_type operator()(const argument_type &a) const
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

// MESH CLASS
class Mesh
{
protected:
	int id, dimension;

public:
	long N, Nel, Ngh, Ned;
	static int number_meshes;
	std::vector<std::shared_ptr<Element>> elems;
	std::vector<Ghost> ghosts;
	std::vector<Vertice> nodes;
	std::unordered_map<std::vector<long>, std::vector<int>> bc_map;

public:
	Mesh(int);
	virtual ~Mesh();

	int get_id(void);
	int get_dimension(void);
	long get_number_nodes(void);
	long get_number_elements(void);
	void print_element_by_id(long);
	void print_node_by_id(long);
	void print_node_id(long);
	void read_gmsh(const std::string &);
	double get_area(const std::vector<Vertice> &);
	double get_area(const std::vector<long> &, const Vertice &);
	double get_volume(const std::vector<Vertice> &);
	void update_element_neighbors(void);
	void update_physical_tags(void);
	void mark_boundaries(void);
	void to_vtk(const std::string &);
	long get_closest(const Node &, std::vector<long> &);
	double get_residue_norm(int);
	//void calculate_jacobians(std::shared_ptr<Element>&);

private:
	void append_elem_to_nodes(const std::shared_ptr<Element> &);
	void append_boundary_face(const int, const int, const std::vector<std::string> &);

	std::istream &get_line(std::istream &, std::string &);
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

	void create_kdtree(void);
	void build_kdtree(long &, std::vector<long> &, long);
	void get_pivot(long &, const std::vector<long> &);
	long mark_fringes(const Vertice &);
	void to_graphviz(void);
	void print_tree(long &, long &, std::string, int &);

private:
	//long _search(const Node&, long&, int&, long&);
	long _search(const Vertice &, long &, int, double, long);
	long _find_element(const Vertice &, std::vector<long> &);
	void _get_kneighbors(long &, std::vector<long> &);
	double _get_distance(const Vertice &, const Vertice &);
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
