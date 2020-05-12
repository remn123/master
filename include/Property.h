#pragma once

#include <memory>
#include <string>
#include <vector>

#include <DVector.h>


class Property 
{
public:
  /* 
    DVector represents a property vector containing all fluid 
    properties or flux properties;
    
    The need of a vector of DVector is due to the high-order 
    polynomial interpolation. So for each (solution/flux) node
    a DVector is allocated as an element inside this vector.

    Although, for flux properties as it has 2 (or 3) directions, 
    the first vector stands for DIRECTION (x, y, z),
    the second vector stands for the number of nodes 
    (as in solution properties).

    At last, for the residue, it has one DVector for each
    node similar to the solution properties.
   */

  // Conservative Properties
  std::vector<DVector> Qsp;                 // Solution at Solution nodes
  std::vector<DVector> Qfp;                 // Solution at Flux nodes

  // Convective Fluxes
  std::vector<std::vector<DVector>> Fcsp;   // Convective Flux at Solution nodes
  std::vector<std::vector<DVector>> Fcfp;   // Convective Flux at Flux nodes 
  
  // Diffusive Fluxes
  std::vector<std::vector<DVector>> Fdsp;   // Diffusive Flux at Solution nodes
  std::vector<std::vector<DVector>> Fdfp;   // Diffusive Flux at Flux nodes
  
  // Gradients
  // Conservative Properties
  std::vector<std::vector<DVector>> dQsp;   // Solution Gradient at Solution nodes
  std::vector<std::vector<DVector>> dQfp;   // Solution Gradient at Flux nodes

  // Convective Fluxes
  std::vector<std::vector<DVector>> dFcsp;   // Convective Flux Gradient at Solution nodes
  std::vector<std::vector<DVector>> dFcfp;   // Convective Flux Gradient at Flux nodes 
  
  // Diffusive Fluxes
  std::vector<std::vector<DVector>> dFdsp;   // Diffusive Flux Gradient at Solution nodes
  std::vector<std::vector<DVector>> dFdfp;   // Diffusive Flux Gradient at Flux nodes
  
  // Residue Vector
  std::vector<DVector> res;

public:
  Property();
  ~Property();

};


