#ifndef POLY_H
#define POLY_H

#include <iostream>
#include <vector>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>

// Poly.h
class Poly
{
public:
  std::vector<double> nodes;
  inline static double MAX_ERROR = 1e-5;
  inline static std::vector<double> memo{1};

public:
  Poly(){std::cout << "Poly Interface is Alive; \n";};
  virtual ~Poly() = default;
    
  virtual void setup(unsigned int)=0;
  virtual std::vector<double> get_nodes(void)=0;
  virtual double get_node(unsigned int)=0;
  unsigned long factorial(unsigned int); // with memoization
  void newton_raphson(double, double, unsigned int, unsigned int, double);
  double Pn(double, double, unsigned int, double); // (alpha, beta, n, x)
  double dPn(double, double, unsigned int, double); // (alpha, beta, n, x)
};

//std::vector<double> Poly::memo{};  // static memo definition


// Chebyshev
class Chebyshev : public Poly
{
public:
  Chebyshev(){std::cout << "Chebyshev is Alive; \n";};
  ~Chebyshev(){std::cout << "Chebyshev is Dead; \n";};
    
  void setup(unsigned int);
  std::vector<double> get_nodes(void);
  double get_node(unsigned int);
};


// Gauss-Legendre 
class GL : public Poly
{
public:
  GL(){std::cout << "GL is Alive; \n";};
  ~GL(){std::cout << "GL is Dead; \n";};
    
  void setup(unsigned int);
  std::vector<double> get_nodes(void);
  double get_node(unsigned int);
};

// Gauss-Legendre-Lobatto
class GLL : public Poly
{
public:
  GLL(){std::cout << "GLL is Alive; \n";};
  ~GLL(){std::cout << "GLL is Dead; \n";};

  void setup(unsigned int);    
  std::vector<double> get_nodes(void);
  double get_node(unsigned int);
};

#endif
