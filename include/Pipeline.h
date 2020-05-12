#pragma once 

#include <memory>
#include <string>
#include <vector>

#include <Dict.h>
#include <Preproc.h>
#include <Posproc.h>
#include <Solver.h>


/* 
  CONTAINERS:
   - Preproc
   - Solver
   - Posproc
 */

/* Interface Pipe */
class IPipe
{
private: 
  dict params;
  std::string message;
  int type;

public:
  IPipe(){};
  virtual ~IPipe(){};

  virtual void setup(dict& params);
  virtual void setup(dict&, const std::shared_ptr<IPipe>&);
  virtual void run(void);
};

/* Templated Pipe */
//template <typename C>
template <template <typename, typename...> class ContainerClass,
          typename Type, typename... Args>
class Pipe : public IPipe
{
private: 
  ContainerClass<Type, Args...> container;
  std::shared_ptr<IPipe> left;
public:
  Pipe();
  ~Pipe();

  void setup(dict& params);
  void setup(dict& params, const std::shared_ptr<IPipe>&);
  void run(void);
};

/* pipe.cpp */
template <template <typename, typename...> class ContainerClass,
          typename Type, typename... Args>
Pipe<ContainerClass, Type, Args...>::Pipe()
{
  std::cout << "New Pipe has been created!\n";
}

template <template <typename, typename...> class ContainerClass,
          typename Type, typename... Args>
Pipe<ContainerClass, Type, Args...>::~Pipe()
{
  std::cout << "A Pipe has been deleted!\n";
}

template <template <typename, typename...> class ContainerClass,
          typename Type, typename... Args>
void Pipe<ContainerClass, Type, Args...>::setup(dict& params, const std::shared_ptr<IPipe>& left)
{
  this->left = left; // address of the last Pipe (first points to nullptr)
  this->params = params;
  this->container = ContainerClass<Type, Args...>{params};
  this->container.setup();
}


template <template <typename, typename...> class ContainerClass,
          typename Type, typename... Args>
void Pipe<ContainerClass, Type, Args...>::setup(dict& params)
{
  this->left = nullptr; // address of the last Pipe (first points to nullptr)
  this->container = ContainerClass<Type, Args...>{params};
  this->container.setup();
}


template <template <typename, typename...> class ContainerClass,
          typename Type, typename... Args>
void Pipe<ContainerClass, Type, Args...>::run(void)
{
  this->container.run();
}


class Pipeline
{
private:
  std::vector<std::shared_ptr<IPipe>> pipes;
public:
  Pipeline(std::vector<std::shared_ptr<IPipe>> vec_pipes) : pipes(vec_pipes){};
  ~Pipeline();

  void setup(std::vector<dict> dicts);
  void start(void);
};


// /* SOLVER */
// template<typename T>
// class Solver : public Pipe, public Method<T>
// {
// private:
//     /* data */
// public:
//     Solver();
//     ~Solver();
// };

// /* POSPROC */
// class Posproc : public Pipe
// {
// private:
//     /* data */
// public:
//     Posproc();
//     ~Posproc();
// };

// /* PIPELINE */
// class Pipeline
// {
// private:
//     std::vector<Pipe> pipes;
//     int verbose;
//     std::string message;
// public:
//     Pipeline(std::vector<Pipe>);
//     ~Pipeline();

//     void set(void);
//     void start(void);
//     void log(const std::string&);
//     void print(const std::string&);
// };

// void Pipeline::set(void)
// {
//     for (auto pipe : this->pipes)
//     {
//         if (pipe.params["verbose"]==1)
//         {
//             this->print(pipe.message);
//         }
//         else if(pipe.params["verbose"]==2)
//         {
//             this->log(pipe.message);
//         }
//         pipe.set();
//     }
// }

// void Pipeline::start(void)
// {
//     for (auto pipe : this->pipes)
//     {
//         if (pipe.params["verbose"]==1)
//         {
//             this->print(pipe.message);
//         }
//         else if(pipe.params["verbose"]==2)
//         {
//             this->log(pipe.message);
//         }
//         pipe.run();
//     }
// }

// void Pipeline::log(const std::string& message)
// {
//     //std::cout << ""
// }
