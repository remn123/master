 
#include <memory>
#include <string>
#include <vector>

#include <Dict.h>
// #include <Preproc.h>
// #include <Solver.h>
// #include <Posproc.h>

// typedef std::unordered_map<const std::string, const void> dict;

/* Interface Pipe */
class IPipe
{
private: 
  dict params;
  std::string message;
  int type;
public:
  IPipe();
  virtual ~IPipe();

  virtual void setup(dict& params);
  virtual void run(void);
};

/* Templated Pipe */
template <typename C>
class Pipe : public IPipe
{
private: 
  C container;
public:
  Pipe();
  ~Pipe();

  void setup(dict& params);
  void run(void);
};

/* pipe.cpp */
template <typename C>
void Pipe<C>::setup(dict& params)
{
  this->container = C{params};
  this->container.setup();
}

template <typename C>
void Pipe<C>::run(void)
{
  this->container.run();
}

/* 
  CONTAINERS:
   - Preproc
   - Solver
   - Posproc
 */

/* Preproc */
class Preproc
{
private:
  
public:
  Preproc(dict&);
  ~Preproc();

  void setup(dict& params);
  void run(void);
};

/* Solver */
template <typename S, typename T>
class Solver
{
private:
  S space;
  T time;
public:
  Solver(dict&);
  ~Solver();

  void setup(dict& params);
  void run(void);
};

/* Posproc */
class Posproc
{
private:
  
public:
  Posproc(dict&);
  ~Posproc();

  void setup(dict& params);
  void run(void);
};


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
