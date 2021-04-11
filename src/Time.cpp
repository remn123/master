#pragma once

#include <memory>
#include <math.h>
#include <iomanip>
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
  for (auto i = 0; i <= stages; i++)
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
      this->alpha = {std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                         1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                         0.44437049406734, 0.55562950593266, 0.0, 0.0, 0.0, 0.0,
                                         0.62010185138540, 0.0, 0.37989814861460, 0.0, 0.0, 0.0,
                                         0.17807995410773, 0.0, 0.0, 0.82192004589227, 0.0, 0.0,
                                         0.00683325884039, 0.0, 0.51723167208978, 0.12759831133288, 0.34833675773694, 0.0}};
      this->beta = {std::vector<double>{0.0,              0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.39175222700392, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.36841059262959, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.25189177424738, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.54497475021237, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.08460416338212, 0.22600748319395, 0.0}};
    }
    break;
  default:
    this->alpha = {std::vector<double>{1.0,              0.0,              0.0, 0.0, 0.0,
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


template <typename Method>
Time<Method>::~Time()
{
}

template <>
Time<Explicit::SSPRungeKutta>::~Time()
{
}


template <>
void Time<Explicit::SSPRungeKutta>::update(std::shared_ptr<Mesh> &mesh,
                                           std::function<void(std::shared_ptr<Mesh> &)> solve)
{
  // Strong-Stability-Preserving Runge Kutta(s,o) s-stage, o-order
  double L1 = 0.0;
  double L2 = 0.0;
  double Linf = 0.0;

  // U(0)   = U(n)
  // for i in [1, ..., s]; for k in [0, 1, ..., i-1];
  // U(i)   = U(0) + dt*SUM(cik*Residue(U(k)))
  // U(n+1) = U(s)

  //this->dt = CFL * dx / (U * (2.0 * p + 1.0)); // TO DO
  this->dt = 0.001 * CFL;

  // Diferentes abordagens (preproc):
  // - const dx -> weak
  // - U = max(u) ? evolve during simulation
  // - dx -> sqrt(e->area()); ... etc
  // -

  // Deep Copy U(n) solution into a contiguous vector
  this->read_solution(mesh, 0);
  this->read_residue(mesh, 0);
  //L2 = mesh->get_residue_norm(1); // L2-norm
  //Linf = mesh->get_residue_norm(2); // Linf-norm
  //if (this->iter % 100 == 0)
  //std::cout << "Iter[" << iter << "]: Residue " << log10(L2) << std::endl;
  // std::cin.get();
  for (size_t i = 1; i <= this->stages; i++)
  {
    this->res_sum = 0.0 * this->res_sum; // resets res_sum
    //std::cout << "res_sum += ";
    for (size_t k = 0; k < i; k++)
    {
      
      //std::cout << " c(" << i << ", " << k << ") * res[" << k << "] + \n";
      //this->res_sum += this->c(i, k) * this->res[k];
      //std::cout << this->alpha[i*(this->stages+1) + k] << " * U(" << k << ") + " << this->beta[i*(this->stages+1) + k] << " * dt * R(" << k << ")\n";
      this->res_sum += this->alpha[i*(this->stages+1) + k] * this->Q[k] + dt * this->beta[i*(this->stages+1) + k]*this->res[k];
      // if (k == i - 1 && k > 0)
      // {
      
      //}
      // i = 5
      // k = 0
      // i*(s+1)+k => 5*(5+1)+0 => 30
      // i = 1
      // k = 0
      // i*(s+1)+k => 1*(5+1)+0 => 6
    }
    // for (auto &qk : this->Q[0])
    //   std::cout << "qk = " << qk << "\n";
    // for (auto &qk : this->res_sum)
    //   std::cout << "res = " << qk << "\n";
    if (i == this->stages)
    { 
      this->res_sum = 0.517231671970585*this->Q[2] + 0.096059710526147*this->Q[3] + 0.063692468666290*dt*this->res[3] + 0.386708617503268*this->Q[4] + 0.226007483236906*dt*this->res[4];
      this->Q[i] = this->res_sum;
      // writes solution into all mesh elements
      this->write_solution(mesh, i);
      solve(mesh);
      break;
    }
    this->Q[i] = this->res_sum;
    // writes solution into all mesh elements
    this->write_solution(mesh, i);
    // recalculates residue
    solve(mesh);
    // reads residue from mesh
    this->read_residue(mesh, i);
    //std::cout << "U[" << i << "] = res_sum; \n";
    // for (auto &qk : this->Q[i])
    //   std::cout << "this->Q[" << i << "] = " << qk << "\n";
    //std::cin.get();
  } 

  // writes solution into all mesh elements
  //this->write_solution(mesh, this->stages);
  // re-calculates residue
  //solve(mesh);
  this->iter++;

  L1 = mesh->get_residue_norm(0); // L1-norm
  L2 = mesh->get_residue_norm(1); // L2-norm
  Linf = mesh->get_residue_norm(2); // Linf-norm
  if (this->iter % 100 == 0)
    std::cout << "Iter[" << iter << "]: L1 = " << std::log10(L1) << "; L2 = " << log10(L2) << "; Linf = " << log10(Linf) << std::endl;
  //std::cin.get();
}

// Main time step loop
template <typename Method>
void Time<Method>::loop(std::shared_ptr<Mesh> &mesh,
                        std::function<void(std::shared_ptr<Mesh> &)> solve,
                        const std::string &filename,
                        std::function<void(const std::shared_ptr<Mesh> &, const std::string &)> to_vtk)
{
  while (this->iter <= this->MAX_ITER)
  {
    this->update(mesh, solve);
    if (this->iter % 10 == 0)
    {
      std::string tstamp = std::to_string(this->iter);
      tstamp.insert(tstamp.begin(), 5 - tstamp.length(), '0');
      this->save(mesh, filename + tstamp + std::string{".vtk"}, to_vtk);
      //std::cin.get();
    }
  }
}

template <typename Method>
void Time<Method>::save(const std::shared_ptr<Mesh> &mesh,
                        const std::string &filename,
                        std::function<void(const std::shared_ptr<Mesh> &, const std::string &)> to_vtk)
{
  to_vtk(mesh, filename);
}

template <typename Method>
void Time<Method>::read_solution(const std::shared_ptr<Mesh> &mesh, size_t k)
{
  size_t index = 0;

  for (auto &e : mesh->elems)
  {
    for (auto s_index=0; s_index<e->computational->Qsp.size(); s_index++)
    {
      for (auto j=0; j<e->computational->Qsp[s_index].size(); j++)
      {
        this->Q[k][index] = e->computational->Qsp[s_index][j];
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
    for (auto s_index=0; s_index<e->computational->Qsp.size(); s_index++)
    {
      for (auto j=0; j<e->computational->Qsp[s_index].size(); j++)
      {
        e->computational->Qsp[s_index][j] = this->Q[k][index];
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
    for (auto s_index=0; s_index<e->computational->res.size(); s_index++)
    {
      for (auto j=0; j<e->computational->res[s_index].size(); j++)
      {
        //std::cout << "e[" << e->id << "]->computational->res[" << s_index << "][" << j << "] = " << e->computational->res[s_index][j] << "\n";
        this->res[k][index] = e->computational->res[s_index][j];
        //std::cout << "this->res[" << k << "][" << index << "] = " << this->res[k][index] << "\n";
        index++;
      }
    }
  }
  //std::cin.get();
}

template <typename Method>
double Time<Method>::c(int i, int k)
{
  //size_t n = i + k;
  double c_ = 0.0;
  //std::cout << "c[" << i << ", " << k << "] = beta[" << i << ", " << k << "] + ";

  // c(i,k) = SUM(j, a(i, j)*c(j, k) ) + b(i,k)

  // i=1, k=0
  // k+1 = 1
  // i-1 = 1-1 = 0
  // c(1,0) = b(1, 0) 

  // i=2, k=0
  // k+1 = 1
  // i-1 = 2-1 = 1
  // c(2,0) = a(2, 1)*c(1,0) + b(2, 0)

  // i=2, k=1
  // k+1 = 1+1 = 2
  // i-1 = 2-1 = 1
  // c(2,1) =  b(2,1)

  // i=3, k=0
  // k+1 = 1
  // i-1 = 3-1 = 2
  // c(3,0) = a(3, 1)*c(1, 0) + a(3, 2)*c(2, 0) + b(3, 0)
  // c(3,0) = (a(3, 2)*a(2, 1) + a(3, 1) )*b(1,0)  + a(3, 2)*b(2, 0) + b(3, 0)

  c_ = this->beta[i * (stages+1) + k];

  for (auto j = k + 1; j < i; j++)
  {
    //std::cout << " alpha[" << i << ", " << j << "] * c[" << j << ", " << k << ")  +";
    c_ += this->alpha[i * (stages+1) + j] * this->c(j, k);
  }
  //std::cin.get();
  //std::cout << "\n";
  //std::cout << "c_ = " << c_ << "\n";
  return c_;
}
