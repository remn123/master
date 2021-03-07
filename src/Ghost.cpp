#pragma once

#include <memory>
#include <Ghost.h>
#include <Property.h>

Ghost::Ghost(const long element_id, const long edge_id, const int local_id, const int type, const int group)
{
  Ghost::num_ghosts++;
  this->id = Ghost::num_ghosts;

  this->physical = std::make_shared<Property>();
  this->computational = std::make_shared<Property>();

  this->elm_id = element_id;
  this->edg_id = edge_id;
  this->lr_edge = local_id;
  this->type = type;
  this->group = group;
}

Ghost::~Ghost()
{
  //Ghost::num_ghosts--;
}