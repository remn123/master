#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <Node.h>
#include <DVector.h>

// Nodes
Node::Node(const std::vector<std::string> &coordinates)
{
  this->id = -1;

  // transforming vector<string> to vector<double> using lambda
  std::transform(coordinates.begin(), coordinates.end(), std::back_inserter((this->coords)),
                 [](const std::string &str) { return std::stod(str); });
}

Node::Node(double x, double y, double z)
{
  this->id = -1;
  this->coords.resize(3);

  this->coords[0] = x;
  this->coords[1] = y;
  this->coords[2] = z;
}

Node::Node(const Node &n)
{
  this->id = n.id;
  this->coords.resize(3);

  this->coords[0] = n.coords[0];
  this->coords[1] = n.coords[1];
  this->coords[2] = n.coords[2];
}

Node::Node()
{
  this->id = -1;
  this->coords.resize(3);

  this->coords[0] = 0.0;
  this->coords[1] = 0.0;
  this->coords[2] = 0.0;
}

Node::~Node(void)
{
}

void Node::print_coords(void)
{
  std::cout << "Node(" << this->id << "): ";
  for (const auto &v : (this->coords))
    std::cout << v << " ";
  std::cout << std::endl;
}

// Vertices
Vertice::Vertice(const std::vector<std::string> &coordinates) : Node(coordinates)
{
  num_nodes++;
  this->id = num_nodes;

  this->left = -1;
  this->right = -1;
}

Vertice::Vertice(double x, double y, double z) : Node(x, y, z)
{
  this->id = -1;
  this->left = -1;
  this->right = -1;
}

Vertice::Vertice() : Node()
{
  this->id = -1;
  this->left = -1;
  this->right = -1;
}

Vertice::~Vertice(void)
{
}

void Vertice::print_elements(void)
{
  std::cout << "Vertice(" << this->id << ").elements: ";
  for (const auto &e : (this->elems))
    std::cout << e << " ";
  std::cout << std::endl;
}

// fNodes
fNode::fNode(const std::vector<std::string> &coordinates) : Node(coordinates)
{
  num_nodes++;
  this->id = num_nodes;

  this->right = -1;
  this->analytical_solution = DVector{};
  this->has_analytical_solution = false;
  this->solved = false;
}

fNode::fNode(double x, double y, double z) : Node(x, y, z)
{
  this->id = -1;
  this->right = -1;
  this->analytical_solution = DVector{};
  this->has_analytical_solution = false;
  this->solved = false;
}

fNode::fNode(long id, long local, const Node &n)
{
  this->id = id;
  this->local = local;
  this->right = -1;

  this->coords[0] = n.coords[0];
  this->coords[1] = n.coords[1];
  this->coords[2] = n.coords[2];
  this->analytical_solution = DVector{};
  this->has_analytical_solution = false;
  this->solved = false;
}

fNode::fNode(long id, long local, long right, const std::vector<double> &coords)
{
  this->id = id;
  this->local = local;
  this->right = right;

  this->coords[0] = coords[0];
  this->coords[1] = coords[1];
  this->coords[2] = coords[2];
  this->analytical_solution = DVector{};
  this->has_analytical_solution = false;
  this->solved = false;
}

fNode::fNode() : Node()
{
  this->id = -1;
  this->right = -1;
  this->analytical_solution = DVector{};
  this->has_analytical_solution = false;
  this->solved = false;
}

fNode::~fNode(void)
{
}