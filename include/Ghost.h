#pragma once

#include <memory>
#include <string>
#include <vector>

#include <Property.h>

class Ghost
{
public:
  std::shared_ptr<Property> physical;
  std::shared_ptr<Property> computational;

  inline static long num_ghosts = -1;
  long id;
  long elm_id;
  long edg_id;
  int local_edg_id;
  int group;
  int type;
public:
  Ghost(const long, const long, const int, const int, const int);
  ~Ghost();

private:

};