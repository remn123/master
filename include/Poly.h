#ifndef POLY_H
#define POLY_H

#include <iostream>
#include <vector>


// Poly.h
class Poly
{
public:
    std::vector<double> nodes;
    double MAX_ERROR;
    inline static std::vector<double> memo{1};

public:
    Poly(){std::cout << "Poly Interface is Alive; \n"; MAX_ERROR = 1e-16;};
    virtual ~Poly() = default;
    
    virtual void setup(int)=0;
    virtual std::vector<double> get_nodes(void)=0;
    virtual double get_node(int)=0;
    long factorial(int); // with memoization
    void newton_raphson(double, double, int, double);
    double Pn(double, double, int, double); // (alpha, beta, n, x)
    double dPn(double, double,int, double); // (alpha, beta, n, x)
};

//std::vector<double> Poly::memo{};  // static memo definition


// Chebyshev
class Chebyshev : public Poly
{
public:
    Chebyshev(){std::cout << "Chebyshev is Alive; \n";};
    ~Chebyshev(){std::cout << "Chebyshev is Dead; \n";};
    
    void setup(int);
    std::vector<double> get_nodes(void);
    double get_node(int);
};


// Gauss-Legendre 
class GL : public Poly
{
public:
    GL(){std::cout << "GL is Alive; \n";};
    ~GL(){std::cout << "GL is Dead; \n";};
    
    void setup(int);
    std::vector<double> get_nodes(void);
    double get_node(int);
};




// Gauss-Legendre-Lobatto
class GLL : public Poly
{
public:
    GLL(){std::cout << "GLL is Alive; \n";};
    ~GLL(){std::cout << "GLL is Dead; \n";};
    
    void setup(int);    
    std::vector<double> get_nodes(void);
    double get_node(int);
};

#endif
