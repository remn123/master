#pragma once

#include <algorithm>
#include <fstream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <Dummies.h>
#include <Mesh.h>
// #include <Helpers.h>
// #include <Node.h>
// #include <Poly.h>
#include <SD.h>

// explicit instances
template class TimeInterface<Explicit::RungeKutta4th>;

template <>
void TimeInterface::update(std::shared_ptr<Mesh> mesh,
                           std::shared_ptr<SD> sd)
{
  TimeInterface::update
}
