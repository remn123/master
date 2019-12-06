#ifndef HELPERS_H
#define HELPERS_H

#include <iostream>
#include <vector>
#include <math.h>
//#include <Poly.h>


// Helpers.h
template <typename P>
class Helpers
{
private:
  static P pol;
public:

  Helpers()=delete; // not instatiable
  virtual ~Helpers(){std::cout<<"Helpers is Dead;\n";};

  static void init(void);
  static void set_nodes(int);
  static void set_nodes(std::vector<double>);
  static std::vector<double> get_nodes(void);
  static double get_node(int);
  static void print_nodes(void);
  static void delete_nodes(void);
};


template <class P> P Helpers<P>::pol;

// Helpers.cpp
template <typename P>
void Helpers<P>::print_nodes(void)
{
  int i=0;
  for (auto n : Helpers<P>::pol.get_nodes())
  {
    i++;
    std::cout << "Node (" << i << "): " << n << "\n";
  }
}

template <typename P>
std::vector<double> Helpers<P>::get_nodes(void)
{
  return Helpers<P>::pol.get_nodes();
}

template <typename P>
double Helpers<P>::get_node(int k)
{
  return Helpers<P>::pol.get_node(k);
}

template <typename P>
void Helpers<P>::set_nodes(int n)
{
  if (n<1)
  {
    std::cout << "n must be greater or equal to 1" << "\n";
  }
  else
  {
    Helpers<P>::pol.setup(n);
  }
}

template <typename P>
void Helpers<P>::set_nodes(std::vector<double> nodes)
{
  if (nodes.size()<1)
  {
    std::cout << "Node vector is empty!" << "\n";
  }
  else
  {
    Helpers<P>::pol.setup(nodes);
  }
}

template <typename P>
void Helpers<P>::init(void)
{
  // stuff
  std::cout << "Initializing Helper functions!" << "\n";
}

template <typename P>
void Helpers<P>::delete_nodes(void)
{
  std::cout << "Freeing memory!" << "\n";
  Helpers<P>::pol.delete_nodes();
}


#endif // HELPERS_H


