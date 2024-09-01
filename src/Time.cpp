#pragma once

#include <chrono>
#include <iomanip>
#include <fstream>
#include <future>
#include <math.h>
#include <memory>
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
                                    int stages, int order, size_t size,
                                    int p, double Uinf, double dx, const std::string log_path)
{
  this->CFL = CFL;
  this->MAX_ITER = MAX_ITER;
  this->stages = stages;
  this->order = order;
  //this->dt = -1.0;
  this->iter = 0;
  this->log_path = log_path + "/log/log.txt";
  
  // Ringleb
  //auto deltat = (dx/Uinf);
  auto deltat = 1.0 / (Uinf * (2.0 * p + 1.0));
  //auto deltat = 10.0 / Uinf;
  
  // Vortex
  // auto c_inf = std::sqrt(1.4);
  // auto deltat = (5.0*c_inf);
  this->tinf = deltat;
  //auto deltat = 0.1/Uinf;
  this->dt = this->CFL*deltat;
  std::cout << "dt = " << this->dt << "\n";

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
void Time<Explicit::SSPRungeKutta>::update(
  std::shared_ptr<Mesh> &mesh,
  std::function<void(std::shared_ptr<Mesh> &)> solve,
  std::function<std::vector<double>(std::shared_ptr<Mesh> &, double)> calculate_error
  )
{
  // start timer
  auto start = std::chrono::high_resolution_clock::now();

  // Strong-Stability-Preserving Runge Kutta(s,o) s-stage, o-order
  double L1 = 0.0;
  double L2 = 0.0;
  double Linf = 0.0;

  // U(0)   = U(n)
  // for i in [1, ..., s]; for k in [0, 1, ..., i-1];
  // U(i)   = U(0) + dt*SUM(cik*Residue(U(k)))
  // U(n+1) = U(s)

  //this->dt = CFL * dx / (U * (2.0 * p + 1.0)); // TO DO
  
  //this->dt = 0.01 * CFL;

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

  this->Q[1] = this->Q[0] + dt * this->res[0];
  // writes solution into all mesh elements
  this->write_solution(mesh, 1);
  // recalculates residue
  solve(mesh);
  // reads residue from mesh
  this->read_residue(mesh, 1);

  this->Q[2] = (3.0/4.0)*this->Q[0] + (1.0/4.0)*this->Q[1] + (1.0/4.0)*dt * this->res[1];
  // writes solution into all mesh elements
  this->write_solution(mesh, 2);
  // recalculates residue
  solve(mesh);
  // reads residue from mesh
  this->read_residue(mesh, 2);

  this->Q[3] = (1.0/3.0)*this->Q[0] + (2.0/3.0)*this->Q[2] + (2.0/3.0)*dt * this->res[2];
  // writes solution into all mesh elements
  this->write_solution(mesh, 3);
  // recalculates residue
  solve(mesh);
  // reads residue from mesh
  this->read_residue(mesh, 3);
  
  L1 = mesh->get_residue_norm(0); // L1-norm
  L2 = mesh->get_residue_norm(1); // L2-norm
  Linf = mesh->get_residue_norm(2); // Linf-norm
  // vortex
  auto computational_time = dt*iter;
  auto physical_time = computational_time/this->tinf;
  auto errors = calculate_error(mesh, physical_time);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  if (this->iter % 10 == 0){
    auto residue_log = "Iter[" + std::to_string(iter) + "|" + std::to_string(physical_time) + 
                      "]: time/iter (" + std::to_string(duration.count()) +
                       " ms): L1 = " + std::to_string(std::log10(L1)) +
                       "; L2 = " + std::to_string(std::log10(L2)) + 
                       "; Linf = " + std::to_string(std::log10(Linf)) +
                       "; RhoL2 = " + std::to_string(std::log10(errors[0])) +
                       "; EntropyL2 = " + std::to_string(std::log10(errors[1])) + "\n";
    // auto residue_log = "Iter[" + std::to_string(iter) + 
    //                   "]: time/iter (" + std::to_string(duration.count()) +
    //                    " ms): L1 = " + std::to_string(std::log10(L1)) +
    //                    "; L2 = " + std::to_string(std::log10(L2)) + 
    //                    "; Linf = " + std::to_string(std::log10(Linf)) + "\n";
    //std::cout << residue_log;
    this->save_log(residue_log);
  }
    
  
  if (isinf(std::abs(log10(L2)))) throw "[Error]: Time iteration has diverged!\n";
  this->iter++;
  //std::cin.get();
}

template <>
void Time<Explicit::SSPRungeKutta>::update(
  std::shared_ptr<Static_Mesh> &background_msh,
  std::shared_ptr<Static_Mesh> &nearbody_msh,
  std::function<void(std::shared_ptr<Mesh> &)> solve,
  std::function<std::vector<double>(std::shared_ptr<Mesh> &, double)> calculate_error,
  std::function<void(std::shared_ptr<Static_Mesh> &, 
  const std::shared_ptr<Static_Mesh> &)> communicate_data)
{
  // start timer
  auto start = std::chrono::high_resolution_clock::now();

  // Strong-Stability-Preserving Runge Kutta(s,o) s-stage, o-order
  double L1 = 0.0;
  double L2 = 0.0;
  double Linf = 0.0;

  // U(0)   = U(n)
  // for i in [1, ..., s]; for k in [0, 1, ..., i-1];
  // U(i)   = U(0) + dt*SUM(cik*Residue(U(k)))
  // U(n+1) = U(s)

  auto background_msh_ = std::static_pointer_cast<Mesh>(background_msh);
  auto nearbody_msh_   = std::static_pointer_cast<Mesh>(nearbody_msh);
  std::vector<std::shared_ptr<Mesh>> meshes = {background_msh_, nearbody_msh_};
  std::future<void> m_futures;

  // Deep Copy U(n) solution into a contiguous vector
  this->read_solution(background_msh, nearbody_msh, 0);
  try {
    this->read_residue(background_msh, nearbody_msh, 0);
  } catch (const char* msg) {
    std::cerr << msg;
    throw;
  }
  

  this->Q[1] = this->Q[0] + dt * this->res[0];
  // writes solution into all mesh elements
  this->write_solution(background_msh, nearbody_msh, 1);
  // Communicate data
  communicate_data(background_msh, nearbody_msh);
  communicate_data(nearbody_msh, background_msh);
  // recalculates residue
  for (auto& mesh : meshes)
    m_futures = std::async(std::launch::async, solve, std::ref(mesh));
  m_futures.get();
  
  // reads residue from meshes
  try {
    this->read_residue(background_msh, nearbody_msh, 1);
  } catch (const char* msg) {
    std::cerr << msg;
    throw;
  }
  

  this->Q[2] = this->Q[0]*(3.0/4.0) + this->Q[1]*(1.0/4.0) + this->res[1]*((1.0/4.0)*dt);
  // writes solution into all mesh elements
  this->write_solution(background_msh, nearbody_msh, 2);
  // Communicate data
  communicate_data(background_msh, nearbody_msh);
  communicate_data(nearbody_msh, background_msh);
  // recalculates residue
  for (auto& mesh : meshes)
    m_futures = std::async(std::launch::async, solve, std::ref(mesh));
  m_futures.get();

  // reads residue from meshes
  try {
    this->read_residue(background_msh, nearbody_msh, 2);
  } catch (const char* msg) {
    std::cerr << msg;
    throw;
  }

  this->Q[3] = (1.0/3.0)*this->Q[0] + (2.0/3.0)*this->Q[2] + (2.0/3.0)*dt * this->res[2];
  // writes solution into all mesh elements
  this->write_solution(background_msh, nearbody_msh, 3);
  // Communicate data
  communicate_data(background_msh, nearbody_msh);
  communicate_data(nearbody_msh, background_msh);
  
  // recalculates residue
  for (auto& mesh : meshes)
    m_futures = std::async(std::launch::async, solve, std::ref(mesh));
  m_futures.get();

  // reads residue from mesh
  try {
    this->read_residue(background_msh, nearbody_msh, 3);
  } catch (const char* msg) {
    std::cerr << msg;
    throw;
  }
  
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  L2  = background_msh->get_residue_norm(1); // L2-norm
  L2 += nearbody_msh->get_residue_norm(1); // L2-norm

  auto computational_time = dt*iter;
  auto physical_time = computational_time/this->tinf;
  
  std::vector<double> errors = {0.0, 0.0};
  for (auto& mesh : meshes) {
    auto errors_ = calculate_error(mesh, physical_time);
    errors[0] = errors[0]+errors_[0];
    errors[1] = errors[1]+errors_[1];
  }    
  
  if (this->iter % 10 == 0)
  {
    L1  = background_msh->get_residue_norm(0); // L1-norm
    L1 += nearbody_msh->get_residue_norm(0); // L1-norm
    Linf  = background_msh->get_residue_norm(2); // Linf-norm
    Linf += nearbody_msh->get_residue_norm(2); // Linf-norm
    auto residue_log = "Iter[" + std::to_string(iter) + "|" + std::to_string(physical_time) + 
                      "]: time/iter (" + std::to_string(duration.count()) +
                       " ms): L1 = " + std::to_string(std::log10(L1)) +
                       "; L2 = " + std::to_string(std::log10(L2)) + 
                       "; Linf = " + std::to_string(std::log10(Linf)) +
                       "; RhoL2 = " + std::to_string(std::log10(errors[0])) +
                       "; EntropyL2 = " + std::to_string(std::log10(errors[1])) + "\n";
    //std::cout << residue_log;
    this->save_log(residue_log);
  }
  
  if (isinf(std::abs(log10(L2)))) throw "[Error]: Time iteration has diverged!\n";
  this->iter++;
}

// Main time step loop
template <typename Method>
void Time<Method>::loop(std::shared_ptr<Mesh> &mesh,
                        std::function<void(std::shared_ptr<Mesh> &)> solve,
                        std::function<std::vector<double>(std::shared_ptr<Mesh> &, double)> calculate_error,
                        const std::string &filename,
                        std::function<void(const std::shared_ptr<Mesh> &, const std::string &)> to_vtk)
{
  while (this->iter <= this->MAX_ITER)
  {
    try {
      this->update(mesh, solve, calculate_error);
    } catch (const char* msg) {
      std::cerr << msg;
      throw;
    }
    
    if (this->iter % 500 == 0)
    {
      std::string tstamp = std::to_string(this->iter);
      tstamp.insert(tstamp.begin(), 5 - tstamp.length(), '0');
      this->save(mesh, filename + tstamp + std::string{".vtk"}, to_vtk);
    }
  }
}

template <typename Method>
void Time<Method>::loop(
  std::shared_ptr<Static_Mesh> &background_msh,
  std::shared_ptr<Static_Mesh> &nearbody_msh,
  std::function<void(std::shared_ptr<Mesh> &)> solve,
  std::function<std::vector<double>(std::shared_ptr<Mesh> &, double)> calculate_error,
  std::function<void(std::shared_ptr<Static_Mesh> &, const std::shared_ptr<Static_Mesh> &)> communicate_data,
  const std::string &filename_bkg,
  const std::string &filename_nbd,
  std::function<void(const std::shared_ptr<Mesh> &, const std::string &)> to_vtk)
{
  auto background_msh_ = std::static_pointer_cast<Mesh>(background_msh);
  auto nearbody_msh_   = std::static_pointer_cast<Mesh>(nearbody_msh);

  while (this->iter <= this->MAX_ITER)
  {
    try {
      this->update(background_msh, nearbody_msh, solve, calculate_error, communicate_data);
    } catch (const char* msg) {
      std::cerr << msg;
      throw;
    }
    
    if (this->iter % 1000 == 0)
    {
      std::string tstamp = std::to_string(this->iter);
      tstamp.insert(tstamp.begin(), 5 - tstamp.length(), '0');
      this->save(background_msh_, filename_bkg + tstamp + std::string{".vtk"}, to_vtk);
      this->save(nearbody_msh_, filename_nbd + tstamp + std::string{".vtk"}, to_vtk);
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
void Time<Method>::save_log(const std::string & message)
{
  std::ofstream output;

  output.open(this->log_path, std::ios_base::app);
  output << message;
  output.close();
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
        // if (std::isnan(e->computational->Qsp[s_index][j]))
        //   std::cout << "Iter (" << this->iter << ") \n" << "Solution = e->computational->Qsp[" << s_index << "]["<<j<<"] = " << e->computational->Qsp[s_index][j] << "\n";
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
        // if (std::isnan(this->Q[k][index]))
        //   std::cout << "Iter (" << this->iter << ") \n" << "Solution = this->Q[" << k << "]["<<index<<"] = " << this->Q[k][index] << "\n";
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
        // if (std::isnan(e->computational->res[s_index][j]))
        //   std::cout << "Iter (" << this->iter << ") \n" << "Residue = e->computational->res[" << s_index << "]["<<j<<"] = " << e->computational->res[s_index][j] << "\n";
        
        this->res[k][index] = e->computational->res[s_index][j];
        //std::cout << "this->res[" << k << "][" << index << "] = " << this->res[k][index] << "\n";
        index++;
      }
    }
  }
  //std::cin.get();
}

template <typename Method>
void Time<Method>::read_solution(const std::shared_ptr<Static_Mesh> &background_msh, 
                                 const std::shared_ptr<Static_Mesh> &nearbody_msh, 
                                 size_t k)
{
  size_t index = 0;
  // std::cout << "\nRead BACKGROUND SOLUTION \n";
  // Background
  for (auto &e : background_msh->elems)
  {
    for (auto s_index=0; s_index<e->computational->Qsp.size(); s_index++)
    {
      for (auto j=0; j<e->computational->Qsp[s_index].size(); j++)
      {
        // std::cout << "BKG: " << e->computational->Qsp[s_index][j] << "\n";
        // if (isnan(e->computational->Qsp[s_index][j])) 
        //   throw "BackgroundMesh: Element " + std::to_string(e->id) + " SP " + std::to_string(s_index) + " has nan solution " + std::to_string(j) + "\n";
        // if (isinf(std::abs(e->computational->Qsp[s_index][j]))) 
        //   throw "BackgroundMesh: Element " + std::to_string(e->id) + " SP " + std::to_string(s_index) + " has inf solution " + std::to_string(j) + "\n";
        this->Q[k][index] = e->computational->Qsp[s_index][j];
        index++;
      }
    }
  }
  // std::cin.get();
  // std::cout << "\nRead NEAR-BODY SOLUTION \n";
  // Nearbody
  for (auto &e : nearbody_msh->elems)
  {
    for (auto s_index=0; s_index<e->computational->Qsp.size(); s_index++)
    {
      for (auto j=0; j<e->computational->Qsp[s_index].size(); j++)
      {
        // std::cout << "NRB: " << e->computational->Qsp[s_index][j] << "\n";
        // if (isnan(e->computational->Qsp[s_index][j])) 
        //   throw "NearBodyMesh: Element " + std::to_string(e->id) + " SP " + std::to_string(s_index) + " has nan solution " + std::to_string(j) + "\n";
        // if (isinf(std::abs(e->computational->Qsp[s_index][j]))) 
        //   throw "NearBodyMesh: Element " + std::to_string(e->id) + " SP " + std::to_string(s_index) + " has inf solution " + std::to_string(j) + "\n";
        this->Q[k][index] = e->computational->Qsp[s_index][j];
        index++;
      }
    }
  }
  // std::cin.get();
}

template <typename Method>
void Time<Method>::write_solution(std::shared_ptr<Static_Mesh> &background_msh,
                                  std::shared_ptr<Static_Mesh> &nearbody_msh,
                                  size_t k)
{
  size_t index = 0;
  // std::cout << "\nWrite BACKGROUND SOLUTION \n";
  // Background
  for (auto &e : background_msh->elems)
  {
    for (auto s_index=0; s_index<e->computational->Qsp.size(); s_index++)
    {
      for (auto j=0; j<e->computational->Qsp[s_index].size(); j++)
      {
        // std::cout << "wBKG: " << this->Q[k][index] << "\n";
        e->computational->Qsp[s_index][j] = this->Q[k][index];
        index++;
      }
    }
  }
  // std::cin.get();
  // std::cout << "\nWrite NEAR-BODY SOLUTION \n";
  // Nearbody
  for (auto &e : nearbody_msh->elems)
  {
    for (auto s_index=0; s_index<e->computational->Qsp.size(); s_index++)
    {
      for (auto j=0; j<e->computational->Qsp[s_index].size(); j++)
      {
        // std::cout << "wNRB: " << this->Q[k][index] << "\n";
        e->computational->Qsp[s_index][j] = this->Q[k][index];
        index++;
      }
    }
  }
  // std::cin.get();
}

template <typename Method>
void Time<Method>::read_residue(const std::shared_ptr<Static_Mesh> &background_msh,
                                const std::shared_ptr<Static_Mesh> &nearbody_msh, size_t k)
{
  size_t index = 0;
  // std::cout << "\nBACKGROUND RESIDUE \n";
  // Background
  for (auto &e : background_msh->elems)
  {
    // std::cout << "[BACKGROUND] Element: " << e->id << "\n";
    // std::cout << "               Boundary: " << e->boundary << "\n";
    // std::cout << "               Fringe: " << e->fringe << "\n";
    // std::cout << "               Nodes: \n" ;
    // std::cout << "                   0: (" << background_msh->nodes[e->nodes[0]].coords[0] << ", " <<  background_msh->nodes[e->nodes[0]].coords[1] << ")\n";
    // std::cout << "                   1: (" << background_msh->nodes[e->nodes[1]].coords[0] << ", " <<  background_msh->nodes[e->nodes[1]].coords[1] << ")\n";
    // std::cout << "                   2: (" << background_msh->nodes[e->nodes[2]].coords[0] << ", " <<  background_msh->nodes[e->nodes[2]].coords[1] << ")\n";
    // std::cout << "                   3: (" << background_msh->nodes[e->nodes[3]].coords[0] << ", " <<  background_msh->nodes[e->nodes[3]].coords[1] << ")\n";
    // std::cout << "               Edges: \n";
    // std::cout << "                      0 -> " << e->edges[0].ghost << "\n";
    // std::cout << "                      1 -> " << e->edges[1].ghost << "\n";
    // std::cout << "                      2 -> " << e->edges[2].ghost << "\n";
    // std::cout << "                      3 -> " << e->edges[3].ghost << "\n";
    // std::cout << "               Ghosts: \n";
    // if (e->edges[0].ghost != -1)
    // {
    //   std::cout << "                      0 -> 0: " << background_msh->ghosts[e->edges[0].ghost].fnodes[0].donor << "\n";
    //   std::cout << "                      0 -> 1: " << background_msh->ghosts[e->edges[0].ghost].fnodes[1].donor << "\n";
    // }
    // if (e->edges[1].ghost != -1)
    // {
    //   std::cout << "                      1 -> 0: " << background_msh->ghosts[e->edges[1].ghost].fnodes[0].donor << "\n";
    //   std::cout << "                      1 -> 1: " << background_msh->ghosts[e->edges[1].ghost].fnodes[1].donor << "\n";
    // }
    // if (e->edges[2].ghost != -1)
    // {
    //   std::cout << "                      2 -> 0: " << background_msh->ghosts[e->edges[2].ghost].fnodes[0].donor << "\n";
    //   std::cout << "                      2 -> 1: " << background_msh->ghosts[e->edges[2].ghost].fnodes[1].donor << "\n";
    // }
    // if (e->edges[3].ghost != -1)
    // {
    //   std::cout << "                      3 -> 0: " << background_msh->ghosts[e->edges[3].ghost].fnodes[0].donor << "\n";
    //   std::cout << "                      3 -> 1: " << background_msh->ghosts[e->edges[3].ghost].fnodes[1].donor << "\n";
    // }
    // std::cout << "               Solution: \n";
    // std::cout << "                   Densities: \n";
    // std::cout << "                             " << e->computational->Qsp[0][0] << "\n";
    // std::cout << "                             " << e->computational->Qsp[1][0] << "\n";
    // std::cout << "                             " << e->computational->Qsp[2][0] << "\n";
    // std::cout << "                             " << e->computational->Qsp[3][0] << "\n";
    // std::cout << "                  x-Momentum: \n";
    // std::cout << "                             " << e->computational->Qsp[0][1] << "\n";
    // std::cout << "                             " << e->computational->Qsp[1][1] << "\n";
    // std::cout << "                             " << e->computational->Qsp[2][1] << "\n";
    // std::cout << "                             " << e->computational->Qsp[3][1] << "\n";
    // std::cout << "                  y-Momentum: \n";
    // std::cout << "                             " << e->computational->Qsp[0][2] << "\n";
    // std::cout << "                             " << e->computational->Qsp[1][2] << "\n";
    // std::cout << "                             " << e->computational->Qsp[2][2] << "\n";
    // std::cout << "                             " << e->computational->Qsp[3][2] << "\n";
    // std::cout << "                      Energy: \n";
    // std::cout << "                             " << e->computational->Qsp[0][3] << "\n";
    // std::cout << "                             " << e->computational->Qsp[1][3] << "\n";
    // std::cout << "                             " << e->computational->Qsp[2][3] << "\n";
    // std::cout << "                             " << e->computational->Qsp[3][3] << "\n";
    // std::cout << "               Residue: \n";
    // std::cout << "                   Densities: \n";
    // std::cout << "                             " << e->computational->res[0][0] << "\n";
    // std::cout << "                             " << e->computational->res[1][0] << "\n";
    // std::cout << "                             " << e->computational->res[2][0] << "\n";
    // std::cout << "                             " << e->computational->res[3][0] << "\n";
    // std::cout << "                  x-Momentum: \n";
    // std::cout << "                             " << e->computational->res[0][1] << "\n";
    // std::cout << "                             " << e->computational->res[1][1] << "\n";
    // std::cout << "                             " << e->computational->res[2][1] << "\n";
    // std::cout << "                             " << e->computational->res[3][1] << "\n";
    // std::cout << "                  y-Momentum: \n";
    // std::cout << "                             " << e->computational->res[0][2] << "\n";
    // std::cout << "                             " << e->computational->res[1][2] << "\n";
    // std::cout << "                             " << e->computational->res[2][2] << "\n";
    // std::cout << "                             " << e->computational->res[3][2] << "\n";
    // std::cout << "                      Energy: \n";
    // std::cout << "                             " << e->computational->res[0][3] << "\n";
    // std::cout << "                             " << e->computational->res[1][3] << "\n";
    // std::cout << "                             " << e->computational->res[2][3] << "\n";
    // std::cout << "                             " << e->computational->res[3][3] << "\n";
    for (auto s_index=0; s_index<e->computational->res.size(); s_index++)
    {
      for (auto j=0; j<e->computational->res[s_index].size(); j++)
      {
        // std::cout << "rrBKG: e[" << e->id << "]->computational->res[" << s_index << "][" << j << "] = " << e->computational->res[s_index][j] << "\n";
        // if (isnan(e->computational->res[s_index][j])) 
        //   throw "BackgroundMesh: Element " + std::to_string(e->id) + " SP " + std::to_string(s_index) + " has nan residue\n";
        // if (isinf(std::abs(e->computational->res[s_index][j]))) 
        //   throw "BackgroundMesh: Element " + std::to_string(e->id) + " SP " + std::to_string(s_index) + " has inf residue\n";
        
        this->res[k][index] = e->computational->res[s_index][j];
        // std::cout << "rrBKG: this->res[" << k << "][" << index << "] = " << this->res[k][index] << "\n";
        index++;
      }
    }
  }
  // std::cin.get();
  // std::cout << "\nNEAR-BODY RESIDUE \n";
  // Nearbody
  for (auto &e : nearbody_msh->elems)
  {
    for (auto s_index=0; s_index<e->computational->res.size(); s_index++)
    {
      for (auto j=0; j<e->computational->res[s_index].size(); j++)
      {
        // std::cout << "rrNRB: e[" << e->id << "]->computational->res[" << s_index << "][" << j << "] = " << e->computational->res[s_index][j] << "\n";
        // if (isnan(e->computational->res[s_index][j])) 
        //   throw "NearBodyMesh: Element " + std::to_string(e->id) + " SP " + std::to_string(s_index) + " has nan residue\n";
        // if (isinf(std::abs(e->computational->res[s_index][j]))) 
        //   throw "NearBodyMesh: Element " + std::to_string(e->id) + " SP " + std::to_string(s_index) + " has inf residue\n";
        this->res[k][index] = e->computational->res[s_index][j];
        // std::cout << "rrNRB: this->res[" << k << "][" << index << "] = " << this->res[k][index] << "\n";
        index++;
      }
    }
  }
  // std::cin.get();
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
