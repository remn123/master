#pragma once

#include <vector>
#include <memory>

#include <Dummies.h>
#include <Element.h>
#include <Ghost.h>
#include <Mesh.h>
#include <Node.h>

template <typename Method>
class TimeInterface
{
public:
  double CFL;
  long MAX_ITER;

  TimeInterface(double, long);
  virtual ~TimeInterface();

  void update();
};

