#pragma once

#include <any>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include <Property.h>
#include <Node.h>

enum PhysicalEnum
{
  WALL = 0,
  SUPERSONIC_INLET,
  SUPERSONIC_OUTLET,
};

class Ghost
{
public:
  std::shared_ptr<Property> physical;
  std::shared_ptr<Property> computational;
  std::vector<fNode> fnodes;
  
  inline static std::any analytical_solution = {};
  inline static std::unordered_map<int, std::any> Qbnds = {};
  inline static std::unordered_map<int, std::vector<double>> dQbnds = {};
  inline static std::unordered_map<int, PhysicalEnum> tag_name_map = {};
  inline static long num_ghosts = -1;
  long id;
  long elm_id;
  long edg_id;
  int lr_edge;
  int group;
  int tag;
  PhysicalEnum type;

public:
  Ghost(const long, const long, const int, const int, const int);
  ~Ghost();

private:
};


