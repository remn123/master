/* 
  ******************************** SPECTRAL DIFFERENCE **********************************
  ***************************************************************************************

   * 1) Apply initial condition (Q)
   * 2) Find out flux (FP) and solution (SP) points on numerical space
   * 3) Interpolate solution (Q) from SP to FP
   * 4) Solve Riemann-Problem at each edge/face (Rusanov, Roe, etc) ... flux calculation
   * 5) Interpolate corrected flux solution from FP to SP
   * 6) Calculate Residue
   * 7) Iterate in time (explicit/implicit)

  ***************************************************************************************
*/
#pragma once

#include <Dict.h>


/* Solver */
template <template <typename, typename...> class SpaceClass,
          template <typename> class TimeClass,
          typename SpaceType, typename TimeType, typename... Args>
class Solver
{
private:
  dict params;
  dict space_params;
  dict time_params;

  SpaceClass<SpaceType, Args..> space;
  TimeClass<TimeType> time;
  
public:
  Solver(dict&);
  ~Solver();

  void setup(void);
  void run(void);
};


/* TimeDiscretization */
class TimeDiscretization
{
public:
  double delta_t;
  double CFL;

public:
  TimeDiscretization();
  virtual ~TimeDiscretization();

  virtual void setup(void);
  virtual void step(void);
}

class ExplicitEuler : public TimeDiscretization
{
public:

public:
  ExplicitEuler();
  ~ExplicitEuler();

  void setup(void);
  void step(void);
}

