#pragma once

#include <memory>
#include <math>
// #include <boost/numeric/ublas/vector.hpp>

#include <Dummies.h>
#include <DVector.h>
#include <Mesh.h>
// #include <Helpers.h>
// #include <Node.h>
// #include <Poly.h>
#include <SD.h>
#include <Time.h>

// explicit instances
template class Time<Explicit::SSPRungeKutta>;

// Constructor
template <>
Time<Explicit::SSPRungeKutta>::Time(double CFL, long MAX_ITER,
                                    int stages, int order, size_t size)
{
  this->CFL = CFL;
  this->MAX_ITER = MAX_ITER;
  this->stages = stages;
  this->order = order;
  this->dt = -1.0;
  this->iter = 0;

  // Initialize DVectors
  for (auto i = 0; i < stages; i++)
  {
    this->Q.push_back(DVector(std::vector<double>(size, 0.0)));
    this->res.push_back(DVector(std::vector<double>(size, 0.0)));
  }
  this->res_sum = DVector(std::vector<double>(size, 0.0));

  switch (order)
  {
  case 3:
    if (stages == 3)
    {
      this->alpha = {std::vector<double>{1.0, 0.0, 0.0,
                                         3.0 / 4.0, 1.0 / 4.0, 0.0,
                                         1.0 / 3.0, 0.0, 2.0 / 3.0}};
      this->beta = {std::vector<double>{1.0, 0.0, 0.0,
                                        0.0, 1.0 / 4.0, 0.0,
                                        0.0, 0.0, 2.0 / 3.0}};
    }
    else if (stages == 4)
    {
      this->alpha = {std::vector<double>{1.0, 0.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0, 0.0,
                                         2.0 / 3.0, 0.0, 1.0 / 3.0, 0.0,
                                         0.0, 0.0, 0.0, 1.0}};
      this->beta = {std::vector<double>{1.0 / 2.0, 0.0, 0.0, 0.0,
                                        0.0, 1.0 / 2.0, 0.0, 0.0,
                                        0.0, 0.0, 1.0 / 6.0, 0.0,
                                        0.0, 0.0, 0.0, 1.0 / 2.0}};
    }
    else if (stages == 5)
    {
      this->alpha = {std::vector<double>{1.0, 0.0, 0.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0, 0.0, 0.0,
                                         0.56656131914033, 0.0, 0.43343868085967, 0.0, 0.0,
                                         0.09299483444413, 0.00002090369620, 0.0, 0.90698426185967, 0.0,
                                         0.00736132260920, 0.20127980325145, 0.00182955389682, 0.0, 0.78952932024253}};
      this->beta = {std::vector<double>{0.37726891511710, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.37726891511710, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.16352294089771, 0.0, 0.0,
                                        0.00071997378654, 0.0, 0.0, 0.34217696850008, 0.0,
                                        0.00277719819460, 0.00001567934613, 0.0, 0.0, 0.29786487010104}};
    }

    break;
  case 4:
    if (stages == 5)
    {
      this->alpha = {std::vector<double>{1.0, 0.0, 0.0, 0.0, 0.0,
                                         0.44437049406734, 0.55562950593266, 0.0, 0.0, 0.0,
                                         0.62010185138540, 0.0, 0.37989814861460, 0.0, 0.0,
                                         0.17807995410773, 0.0, 0.0, 0.82192004589227, 0.0,
                                         0.00683325884039, 0.0, 0.51723167208978, 0.12759831133288, 0.34833675773694}};
      this->beta = {std::vector<double>{0.39175222700392, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.36841059262959, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.25189177424738, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.54497475021237, 0.0,
                                        0.0, 0.0, 0.0, 0.08460416338212, 0.22600748319395}};
    }
    break;
  default:
    this->alpha = {std::vector<double>{1.0, 0.0, 0.0, 0.0, 0.0,
                                       0.44437049406734, 0.55562950593266, 0.0, 0.0, 0.0,
                                       0.62010185138540, 0.0, 0.37989814861460, 0.0, 0.0,
                                       0.17807995410773, 0.0, 0.0, 0.82192004589227, 0.0,
                                       0.00683325884039, 0.0, 0.51723167208978, 0.12759831133288, 0.34833675773694}};
    this->beta = {std::vector<double>{0.39175222700392, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.36841059262959, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.25189177424738, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.54497475021237, 0.0,
                                      0.0, 0.0, 0.0, 0.08460416338212, 0.22600748319395}};
    this->order = 4;
    this->stages = 5;
  }
}

template <>
void Time<Explicit::SSPRungeKutta>::update(std::shared_ptr<Mesh> &mesh,
                                           void (*solve)(std::shared_ptr<Mesh> &))
{
  // Strong-Stability-Preserving Runge Kutta(s,o) s-stage, o-order
  double Linf = 0.0;

  // U(0)   = U(n)
  // for i in [1, ..., s]; for k in [0, 1, ..., i-1];
  // U(i)   = U(0) + dt*SUM(cik*Residue(U(k)))
  // U(n+1) = U(s)

  //this->dt = CFL * dx / (U * (2.0 * p + 1.0)); // TO DO
  this->dt = 0.01 * CFL;

  // Diferentes abordagens (preproc):
  // - const dx -> weak
  // - U = max(u) ? evolve during simulation
  // - dx -> sqrt(e->area()); ... etc
  // -

  // Deep Copy U(n) solution into a contiguous vector
  this->read_solution(mesh, 0);
  this->read_residue(mesh, 0);

  for (size_t i = 1; i <= this->stages; i++)
  {
    this->res_sum = 0.0 * this->res_sum; // resets res_sum
    for (size_t k = 0; k < i; k++)
    {
      if (k == i - 1 && k > 0)
      {
        // writes solution into all mesh elements
        this->write_solution(mesh, k);
        // recalculates residue
        solve(mesh);
        // reads residue from mesh
        this->read_residue(mesh, k);
      }
      this->res_sum += this->c(i, k) * this->res[k];
    }
    this->Q[i] = this->Q[0] + dt * this->res_sum;
  }

  // writes solution into all mesh elements
  this->write_solution(mesh, this->stages);
  // re-calculates residue
  solve(mesh);
  this->iter++;

  //L1 = mesh->get_residue_norm(0); // L1-norm
  //L2 = mesh->get_residue_norm(1); // L2-norm
  Linf = mesh->get_residue_norm(2); // Linf-norm
  if (this->iter % 100 == 0)
    std::cout << "Iter[" << iter << "]: Residue " << log10(Linf) << std::endl;
}

// Main time step loop
template <typename Method>
void Time<Method>::loop(std::shared_ptr<Mesh> &mesh,
                        void (*solve)(std::shared_ptr<Mesh> &))
{
  while (this->iter <= this->MAX_ITER)
    this->update(mesh, solve);
}

template <typename Method>
void Time<Method>::read_solution(const std::shared_ptr<Mesh> &mesh, size_t k)
{
  size_t index = 0;

  for (auto &e : mesh->elems)
  {
    for (auto &q : e->computational->Qsp)
    {
      for (auto val : q)
      {
        this->Q[k][index] = val;
        index++;
      }
    }
  }
}

template <typename Method>
void Time<Method>::write_solution(std::shared_ptr<Mesh> &mesh, size_t k)
{
  size_t index = 0;

  for (auto &e : mesh->elems)
  {
    for (auto &q : e->computational->Qsp)
    {
      for (auto &val : q)
      {
        val = this->Q[k][index];
        index++;
      }
    }
  }
}

template <typename Method>
void Time<Method>::read_residue(const std::shared_ptr<Mesh> &mesh, size_t k)
{
  size_t index = 0;

  for (auto &e : mesh->elems)
  {
    for (auto &q : e->computational->res)
    {
      for (auto val : q)
      {
        this->res[k][index] = val;
        index++;
      }
    }
  }
}

template <>
double Time<Explicit::SSPRungeKutta>::c(int i, int k)
{
  size_t n = i + k;
  double c_ = 0.0;

  for (auto j = k + 1; j < i - 1; j++)
  {
    c_ += this->alpha[i * stages + j] * this->c(j, k) + this->beta[i * stages + k];
  }
  return c_;
}
