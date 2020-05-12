#pragma once

#include <Dict.h>

/* Preproc */
template <typename MeshClass, typename... Args>
class Preproc
{
private:
  MeshClass mesh;
  dict params;

public:
  Preproc(dict&);
  ~Preproc();

  void setup(void);
  void run(void);
};