#pragma once

#include <vector>
#include <string>

#include <DVector.h>

// Node Interface
class Node
{
public:
  long id;
  std::vector<double> coords;

public:
  Node(const std::vector<std::string> &);
  Node(double, double, double);
  Node(const Node &);
  Node();
  virtual ~Node();

  void print_coords(void);
};

class Vertice : public Node
{
public:
  inline static long num_nodes = -1;
  long left;
  long right;
  int boundary;
  std::vector<long> elems;

public:
  Vertice(const std::vector<std::string> &);
  Vertice(double, double, double);
  Vertice();
  ~Vertice();

  void print_elements(void);

private:
};

class fNode : public Node
{
public:
  inline static long num_nodes = -1;
  long right;
  long local;
  long donor;
  DVector analytical_solution;
  bool has_analytical_solution;

public:
  fNode(const std::vector<std::string> &);
  fNode(double, double, double);
  fNode(long, long, const Node &);
  fNode(long, long, long, const std::vector<double> &);
  fNode();
  ~fNode();

private:
};
