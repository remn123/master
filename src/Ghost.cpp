#pragma once

#include <Ghost.h>

Ghost::Ghost(const long element_id, const long edge_id, const int local_id, const int type, const int group)
{
  Ghost::num_ghosts++;
  this->id = Ghost::num_ghosts;
  
  this->elm_id = element_id;
  this->edg_id = edge_id;
  this->local_edg_id = local_id;
  this->type = type;
  this->group = group;
}

Ghost::~Ghost()
{
}