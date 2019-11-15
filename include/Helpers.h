#ifndef HELPERS_H
#define HELPERS_H

#include <iostream>
#include <vector>
#include <math.h>

// Poly.h
class Poly
{
public:
    std::vector<double> nodes;
    double MAX_ERROR;
    static std::vector<double> memo;

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

std::vector<double> Poly::memo{};  // static memo definition

// Poly.cpp
void Poly::newton_raphson(double a, double b, int n, double guess)
{
    double xn = 0.0;
    double x0 = guess;
    double dx = 1.0;

    while (dx > this->MAX_ERROR)
    {
        xn = x0 - this->Pn(a,b,n,x0)/this->dPn(a,b,n,x0);
        dx = abs(xn-x0);
        x0 = xn;
    }
    this->nodes.push_back(xn);
}

double Poly::Pn(double a, double b, int n, double x)
{
    if (n == 0)
    {
        return 1.0;
    }
    else if (n == 1)
    {
        return 0.5*(a-b+(a+b+2.0)*x);
    }
    // for n==2 we will use the recurrence relation
    double a1 = 2.0*n*(n+a+b)*(2.0*(n-1.0)+a+b);
    double a2 = (2.0*(n-1.0)+a+b+1.0)*(a*a-b*b);
    double a3 = (2.0*(n-1.0)+a+b)*(2.0*(n-1.0)+a+b+1.0)*(2.0*n+a+b);
    double a4 = 2.0*(n-1.0+a)*(n-1.0+b)*(2.0*n+a+b);

    return (1.0/a1)*((a2+a3*x)*Pn(a,b,n-1,x)-a4*Pn(a,b,n-2,x));
}

double Poly::dPn(double a, double b, int n, double x)
{
    if (n == 0)
    {
        return 0.0;
    }
    else if (n == 1)
    {
        return 0.5*(a+b+2.0);
    }
    // for n==2 we will use the recurrence relation
    double b1 = (2.0*n+a+b)*(1.0-x*x);
    double b2 = n*(a-b-(2.0*n+a+b)*x);
    double b3 = 2.0*(n+a)*(n+b);

    return (b2/b1)*Pn(a,b,n,x) + (b3/b1)*Pn(a,b,n-1,x); 
}

long Poly::factorial(int n)
{
    if(n==0)
    {
        if(this->memo.size() <= 1)
        {
            this->memo[0] = 1;
        }
        
        return this->memo[n];
    }
    else if (this->memo.size() > n)
    {
        std::cout << "Ahá MEMO for n = " << n << "!!!!!\n";
        std::cout << "n! = " << this->memo[n] << "\n";
        return this->memo[n];
    }
    this->memo.resize(n+1);
    this->memo[n] = n*this->factorial(n-1);
    return this->memo[n];
}

// ------------------------------- POLY ------------------------------


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

void Chebyshev::setup(int n)
{
    std::cout << "Setting up Chebyshev Polynomials" << "\n";
    this->nodes.reserve(n);

    double xk = 0.0;
    for (auto k=0; k<n; k++)
    {
        xk = cos((2.0*k+1.0)*M_PI/(2.0*n));
        this->nodes.push_back(xk);
    }
    //this->nodes = a+1;
}

std::vector<double> Chebyshev::get_nodes(void)
{
    return this->nodes;
}
double Chebyshev::get_node(int k)
{
    return this->nodes[k];
}
// ------------------------------- CHEBYSHEV ------------------------------

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


void GL::setup(int n)
{
    std::cout << "Setting up Gauss-Legendre Polynomials" << "\n";
    double guess_k;
    Chebyshev cheb{};
    cheb.setup(n);

    this->nodes.reserve(n);
    for (auto k=0; k<n; k++)
    {
        guess_k = cheb.get_node(k); 
        newton_raphson(0.0, 0.0, k, guess_k); // Legendre Polynomials
    }
    //this->nodes = a;
}

std::vector<double> GL::get_nodes(void)
{
    return this->nodes;
}


double GL::get_node(int k)
{
    return this->nodes[k];
}
// ------------------------------- GL ------------------------------


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


void GLL::setup(int n)
{
    std::cout << "Setting up Gauss-Legendre-Lobatto Polynomials" << "\n";
    
    double guess_k;
    Chebyshev cheb{};
    cheb.setup(n);

    this->nodes.reserve(n);
    this->nodes[0] = -1.0;
    for (auto k=1; k<n-1; k++)
    {
        guess_k = cheb.get_node(k); 
        newton_raphson(1.0, 1.0, k, guess_k);
    }
    this->nodes[n-1] = 1.0;
}

std::vector<double> GLL::get_nodes(void)
{
    return this->nodes;
}

double GLL::get_node(int k)
{
    return this->nodes[k];
}
// ------------------------------- GLL ------------------------------




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
    static std::vector<double> get_nodes(void);
    static void print_nodes(void);
};

template <class P> P Helpers<P>::pol;         // template definition


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
void Helpers<P>::init(void)
{
    // stuff
    std::cout << "Initializing Helper functions!" << "\n";
}

#endif // HELPERS_H


