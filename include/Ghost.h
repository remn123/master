#pragma once

#include <unordered_map>
#include <memory>
#include <string>
#include <vector>

#include <Property.h>
#include <Node.h>
class Ghost
{
public:
  std::shared_ptr<Property> physical;
  std::shared_ptr<Property> computational;
  std::vector<fNode> fnodes;

  static std::unordered_map<int, std::vector<double>> Qbnds;
  static std::unordered_map<int, std::vector<double>> dQbnds;
  inline static long num_ghosts = -1;
  long id;
  long elm_id;
  long edg_id;
  int lr_edge;
  int group;
  int type;
public:
  Ghost(const long, const long, const int, const int, const int);
  ~Ghost();

private:

};