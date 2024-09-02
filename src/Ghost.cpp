#pragma once

#include <memory>
#include <Ghost.h>
#include <Property.h>

Ghost::Ghost(const long element_id, const long edge_id, const int local_id, const int tag, const int group)
{
  Ghost::num_ghosts++;
  this->id = Ghost::num_ghosts;

  this->physical = std::make_shared<Property>();
  this->computational = std::make_shared<Property>();

  this->elm_id = element_id;
  this->edg_id = edge_id;
  this->lr_edge = local_id;
  this->tag = tag;
  this->group = group;
}

Ghost::~Ghost()
{
  
}