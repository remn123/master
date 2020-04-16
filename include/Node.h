#pragma once

#include <vector>
#include <string>


// Node Interface
class Node
{
public:

  long id;
  std::vector<double> coords;

public:
  Node(const std::vector<std::string>&);
  Node(double, double, double);
  Node();
  virtual ~Node();

  void print_coords(void);
};




class Vertice : public Node
{
public:
  static long num_nodes;
  long left;
  long right;
  int boundary;
  std::vector<long> elems;

public:
  Vertice(const std::vector<std::string>&);
  Vertice(double, double, double);
  Vertice();
  ~Vertice();

  void print_elements(void);
  
private:

};


