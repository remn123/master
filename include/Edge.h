#pragma once

#include <memory>
#include <string>
#include <vector>


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