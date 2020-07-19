#pragma once

#include <memory>
#include <boost/numeric/ublas/vector.hpp>

#include <Dummies.h>
#include <Mesh.h>
// #include <Helpers.h>
// #include <Node.h>
// #include <Poly.h>
#include <SD.h>
#include <Time.h>

// explicit instances
template class Time<Explicit::SSPRungeKutta>;


// Constructor
template<>
Time<Explicit::SSPRungeKutta>::Time(double CFL, long MAX_ITER, 
                                    int stages, int order)
{
  this->CFL = CFL;
  this->MAX_ITER = MAX_ITER;
  this->stages = stages;
  this->order = order;
}

template<>
void Time<Explicit::SSPRungeKutta>::update(std::shared_ptr<Mesh>& mesh, 
                                           void * res)
{
  // Strong-Stability-Preserving Runge Kutta(s,o) s-stage, o-order
  
  // U(0)   = U(n)
  // for i in [1, ..., s]; for k in [0, 1, ..., i-1];
  // U(i)   = U(0) + dt*SUM(cik*Residue(U(k))) 
  // U(n+1) = U(s)

  // Deep Copy U(n) solution into a contiguous vector
  this->read_from_mesh(mesh, this->Q);
  this->get_residue(mesh);

  for (auto i = 1; i <= this->stages; i++)
  {
    for (auto k = 0; k < i ; k++)
    { 
      
    }
    this->Qnew = this->Q + dt*
  }
  U(0) = mesh->U;
  U(1) = U(0) + dt*(c10*res(U(0))) 
  U(2) = U(0) + dt*(c10*res(U(0)) + c11*res(U(1)))
  U(2) = U(0) + dt*(c10*res(U(0)) + c11*res(U(1)) + c12*res(U(2)))
  
  
}


template<typename Method>
void Time<Method>::read_from_mesh(const std::shared_ptr<Mesh>& mesh, ublas::vector<double>& vec)
{
  size_t index = 0;

  for (auto& e : mesh->elems)
  {
    for (auto& q : e->computational->Qsp)
    {
      for (auto& val : q)
      {
        vec(index) = val;
        index++;
      }
    }
  }
}

template<typename Method>
void Time<Method>::write_into_mesh(std::shared_ptr<Mesh>& mesh, ublas::vector<double>& vec)
{
  size_t index = 0;

  for (auto& e : mesh->elems)
  {
    for (auto& q : e->computational->Qsp)
    {
      for (auto& val : q)
      {
        val = vec(index);
        index++;
      }
    }
}


template <typename Method>
void Time<Method>::get_residue(const std::shared_ptr<Mesh>& mesh)
{
  size_t index = 0;

  for (auto& e : mesh->elems)
  {
    for (auto& q : e->computational->res)
    {
      for (auto& val : q)
      {
        this->res_sum(index) = val;
        index++;
      }
    }
  }
}