#include <Helpers.h>
#include <Poly.h>

// Poly.cpp
void Poly::newton_raphson(double a, double b, unsigned int n, double guess)
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

double Poly::Pn(double a, double b, unsigned int n, double x)
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

double Poly::dPn(double a, double b, unsigned int n, double x)
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

unsigned long Poly::factorial(unsigned int n)
{
  if(n==0)
  {
    return this->memo[n];
  }
  else if (this->memo.size() > n)
  {
    std::cout << "Ahá MEMO for n = " << n << "!!!!!\n";
    std::cout << "n! = " << this->memo[n] << "\n";
    return this->memo[n];
  }
  
  long tmp = n*this->factorial(n-1);
  
  this->memo.resize(n+1);
  this->memo[n] = tmp;
  return this->memo[n];
}

void Poly::delete_nodes(void)
{
  this->nodes.clear();
}
// ------------------------------- POLY ------------------------------//


// CHebyshev
void Chebyshev::setup(unsigned int n)
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
double Chebyshev::get_node(unsigned int k)
{
  return this->nodes[k];
}
// ------------------------------- CHEBYSHEV ------------------------------

// Gauss-Legendre
void GL::setup(unsigned int n)
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


double GL::get_node(unsigned int k)
{
  return this->nodes[k];
}
// ------------------------------- GL ------------------------------


// Guass-Legendre-Lobatto
void GLL::setup(unsigned int n)
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

double GLL::get_node(unsigned int k)
{
  return this->nodes[k];
}
// ------------------------------- GLL ------------------------------


// Lagrange
void Lagrange::setup(std::vector<double> nodes)
{
  std::cout << "Setting up Lagrange Polynomials" << "\n";
  this->nodes = nodes;
}

std::vector<double> Lagrange::get_nodes(void)
{
  return this->nodes;
}
double Lagrange::get_node(unsigned int k)
{
  return this->nodes[k];
}

double Lagrange::Pn(unsigned int i, double x)
{
  double xi=0.0, xk=0.0, prod=1.0;
  unsigned int n = this->nodes.size();
  if (n > 0)
  {
    for (auto k=0; k<n; k++)
    {
      if (i != k)
      {
        xi = this->nodes[i];
        xk = this->nodes[k];

        prod *= (x - xk)/(xi - xk);
      }
    }
  }
  return prod;
}

double Lagrange::dPn(unsigned int i, double x)
{
  double xi=0.0, xk=0.0, xj=0.0, prod=1.0, sum=0.0;
  unsigned int n = this->nodes.size();
  if (n > 0)
  {
    for (auto j=0; j<n; j++)
    {
      xj = this->nodes[j];
      for (auto k=0; k<n; k++)
      {
        if (k != i && k != j)
        {
          xi = this->nodes[i];
          xk = this->nodes[k];

          prod *= (x - xk)/(xi - xk);
        }
      }
      if (j != i)
      {
        sum += (-xj)/(xi - xj);
      }
    }
  }
   return sum;
}
// ------------------------------- LAGRANGE ------------------------------



