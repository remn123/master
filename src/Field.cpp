#pragma once 
#include <algorithm>
#include <assert.h>
#include <math.h>

#include <Ghost.h>
#include <Field.h>
#include <Numericals.h>

/* Declaring unitary field */
std::vector<double> FIELDS::DEFAULT_FIELD_MAPPING (const Node& n)
{
  std::vector<double> vec(4, 0.0); 

  vec[0] = 1.0; 
  vec[1] = 0.0; 
  vec[2] = 0.0; 
  vec[3] = 1.0/(1.4*0.4); 
  return vec;
};


std::vector<double> FIELDS::GAUSSIAN_FIELD_MAPPING (const Node& n)
{
  std::vector<double> vec(4, 0.0);
  
  // Coordinates
  double x = n.coords[0];
  double y = n.coords[1];
  double mux=0.5, muy=0.5;
  double sigmax=0.1, sigmay=0.1;

  double rho = 1.0;
  double u = (1.0/(sigmax*std::sqrt(2.0*M_PI)))*std::exp(-(x-mux)*(x-mux)/(2.0*sigmax));
  double v = (1.0/(sigmay*std::sqrt(2.0*M_PI)))*std::exp(-(y-muy)*(y-muy)/(2.0*sigmay));
  double E = 1.0;

  vec[0] = rho; 
  vec[1] = rho*u; 
  vec[2] = rho*v; 
  vec[3] = E; 
  
  return vec;
};

std::vector<double> FIELDS::LINEAR_FIELD_MAPPING (const Node& n)
{
  std::vector<double> vec(4, 0.0);
  
  // Coordinates
  double x = n.coords[0];
  double y = n.coords[1];
  double ax=0.5, ay=0.5;
  double bx=0.0, by=0.0;

  double rho = 1.0;
  double u = ax*x + bx;
  double v = ay*y + by;
  double E = 1.0;

  vec[0] = rho; 
  vec[1] = rho*u; 
  vec[2] = rho*v; 
  vec[3] = E; 
  
  return vec;
};


/* CYLINDER FLOW */
std::vector<double> FIELDS::CYLINDER_FIELD_MAPPING (const Node& n)
{
  std::vector<double> vec(4, 0.0);
  
  // Coordinates
  double x = n.coords[0];
  double y = n.coords[1];

  double epsilon = 1E-2; 
  // Center 
  double xc=0.5, yc=0.5;
  // Outer Circle
  double R=0.5;
  // Inner Circle
  double r=0.5/16.0;

  double radius = std::sqrt(std::pow((x-xc), 2.0) + std::pow((y-yc), 2.0));

  double rho=1.0, u=0.2, v=0.0, E=1.0, p=1.0/1.4;
  double a=1.0, umag=1.0, q=1.0;

  // Problem specifications
  double gamma=1.4;
  double Mach=Ghost::Mach;
  
  //Check if it is around the Outer Boundary
  // if (std::abs(R - radius) < epsilon)
  // {
  //   rho = 1.0;
  //   u = Mach;
  //   v = 0.0;

  //   q = u*u + v*v;
  //   umag=std::sqrt(q);
  //   p = 1.0/gamma;
  //   E = p/(gamma - 1.0) + 0.5 * rho * (u * u + v * v);
    
  // } 
  // // Check if it is around the Inner Boundary
  // else if (std::abs(radius - r) < epsilon)
  // {
  //   rho = 1.0;
  //   u = Mach;
  //   v = 0.0;
  //   E = p/(gamma - 1.0) + 0.5 * rho * (u * u + v * v);
  // }
  // E = p/(gamma - 1.0) + 0.5 * rho * (u * u + v * v);

  rho = 1.0;
  u = Mach;
  v = 0.0;

  q = u*u + v*v;
  umag=std::sqrt(q);
  p = 1.0/gamma;
  E = p/(gamma - 1.0) + 0.5 * rho * (u * u + v * v);

  vec[0] = rho; 
  vec[1] = rho*u; 
  vec[2] = rho*v; 
  vec[3] = E; 
  
  return vec;
};

/* IMPLOSION FLOW */
std::vector<double> FIELDS::IMPLOSION_FIELD_MAPPING (const Node& n)
{
  std::vector<double> vec(4, 0.0);
  
  // Coordinates
  double x = n.coords[0];
  double y = n.coords[1];

  double rho = 1.0;
  double u = 0.0;
  double v = 0.0;
  double p = 0.0;
  double E = 1.0;
  
  double L = 0.5;

  // For x+y > 0.15, the initial density and pressure is 1.0, otherwise ρ = 0.125 and P = 0.14
  if (x + y > L)
  {
    p = 1.0;
    rho = 1.0;
  }
  else
  {
    rho = 0.125;
    p = 0.14;
  }

  E = p/0.4;

  vec[0] = rho; 
  vec[1] = rho*u; 
  vec[2] = rho*v; 
  vec[3] = E; 

  return vec;
};

/* RINGLEB FLOW */
std::vector<double> FIELDS::RINGLEB_FIELD_MAPPING (const Node& n)
{
  std::vector<double> vec(4, 0.0);
  boost::numeric::ublas::vector<double> guess(1);
  boost::numeric::ublas::vector<double> u1(1), u2(1), ul(1), uh(1);
  boost::numeric::ublas::vector<double> fl(1), fh(1);
  boost::numeric::ublas::vector<double> velocity(1);
  boost::numeric::ublas::vector<double> f_guess(1);
  boost::numeric::ublas::matrix<double> df_guess(1,1);

  // Coordinates
  double x = n.coords[0];
  double y = n.coords[1];
  //std::cout << "Solving for (" << x << "; " << y << ")\n";
  double gamma = 1.4;
  double q = std::sqrt(x*x + y*y); // initial guess

  // lambda function for ringleb flow x,y, V relation
  auto ringleb_f = [&](const boost::numeric::ublas::vector<double>& u, 
                       boost::numeric::ublas::vector<double>& res) -> bool
  {
    //double V = std::max(0.5, std::min(u(0), 1.5));
    double V = u(0);
    double b = std::sqrt(1.0-0.2*std::pow(V, 2.0));
    double rho = std::pow(b, 5.0);
    // std::cout << "V   = " << V << std::endl;
    // std::cout << "b   = " << b << std::endl;
    // std::cout << "rho = " << rho << std::endl;
    //assert (rho > 0.0);
    double L = (1.0/b) + (1.0/(3.0*std::pow(b, 3.0))) + (1.0/(5.0*std::pow(b, 5.0))) - (0.5*(std::log((1.0+b)/(1.0-b))));

    res(0) = std::pow((x+L/2.0), 2.0) + std::pow(y, 2.0) - std::pow(1.0/(2.0*rho*std::pow(V, 2.0)), 2.0);

    return true;
  };

  // lambda function for ringleb flow x,y, V derivative relation
  auto ringleb_df = [&](const boost::numeric::ublas::vector<double>& u, 
                        boost::numeric::ublas::matrix<double>& res) -> bool
  {
    //double V = std::max(0.5, std::min(u(0), 1.5));
    double V = u(0);
    double b = std::sqrt(1.0-0.2*std::pow(V, 2.0));
    double rho = std::pow(b, 5.0);
    double L = (1.0/b) + (1.0/(3.0*std::pow(b, 3.0))) + (1.0/(5.0*std::pow(b, 5.0))) - (0.5*(std::log((1.0+b)/(1.0-b))));

    double db_dV = -0.2*V/b;
    double drho_dV = 5.0*db_dV*std::pow(b, 4.0);

    double dl1 = -db_dV*std::pow(b, -2.0);
    double dl2 = -db_dV*std::pow(b, -4.0);
    double dl3 = -db_dV*std::pow(b, -6.0);
    double dl4 = -db_dV/((1.0-b)*(1.0+b));
    double dL_dV = dl1 + dl2 + dl3 + dl4;

    double dp1 = (x+L/2.0)*dL_dV;
    double dp2 = (4.0*V*(drho_dV*V + 2.0*rho))/std::pow(2.0*rho*std::pow(V, 2.0), 3.0);
    res(0, 0) = dp1 + dp2;

    return true;
  };

  // Finding
  std::vector<double> guesses = {0.49, 0.75, 1.0, 1.25, 1.5};
  bool solved = false;
  double cur_g = 0.0, cur_fg = 0.0;
  double last_guess=0.0;
  int counter=0;
  for (auto g : guesses)
  {
    guess(0) = g;
    solved = ringleb_f(guess, f_guess);
    if (!solved) throw "ERROR: unvalid guess for ringleb function";
    if (counter>0) 
    {
      if (cur_fg*f_guess(0)<0.0) 
      {
        // std::cout << "last fg: " << cur_fg << "; cur_fg = " << f_guess(0) << "\n";
        // std::cout << "last guess: " << cur_g << "; current guess = " << g << "\n";
        //guess(0)=-(cur_fg*g-f_guess(0)*cur_g)/(f_guess(0)-cur_fg);
        u1(0) = cur_g;
        u2(0) = g;
        fl(0) = cur_fg;
        fh(0) = f_guess(0);
        last_guess=g;
        break;
      }
    }
    cur_fg = f_guess(0);
    cur_g = g;
    counter++;
  }
  velocity(0) = 0.0;

  solved = false;
  double ds = 1E-5;
  double b = 0.0;
  double rho = 0.0;
  //q = 0.5;
  q = u1(0);
  while (!solved)
  {
    // if (fl(0) < 0.0) 
    // {
    //   ul(0)=u1(0);
    //   uh(0)=u2(0);
    // } else {
    //   uh(0)=u1(0);
    //   ul(0)=u2(0);
    // }
    guess(0) = q;
    solved = newton_raphson(guess, ringleb_f, ringleb_df, velocity);
    b = std::sqrt(1.0-0.2*std::pow(velocity(0), 2.0));
    rho = std::pow(b, 5.0);
    // if (abs(x+0.277229)<1E-5 && abs(y-0.814547)<1E-5)
    // {
    //   solved = ringleb_df(guess, df_guess);
    //   solved = ringleb_f(guess, f_guess);
    //   std::cout << "q  = " << q << "\n";
    //   std::cout << "u  = " << guess << "\n";
    //   std::cout << "f  = " << f_guess << "\n";
    //   std::cout << "df = " << df_guess << "\n";
    //   std::cout << "sol  = " << velocity << "\n";
    //   std::cin.get();
    // }
    if (rho <= 0.0 || velocity(0) > 2.236 || velocity(0) < 0.0) solved=false;
    if (q > last_guess) throw "ERROR: domain limit has been reached!";
    q += ds;
  }
  
  assert (solved==true);
  q = velocity(0);
  //std::cout << "(x, y, V) = (" << x << ", "<< y << ", " << q << ")\n";
  b = std::sqrt(1.0-0.2*std::pow(q, 2.0));
  double L = (1.0/b) + (1.0/(3.0*std::pow(b, 3.0))) + (1.0/(5.0*std::pow(b, 5.0))) - (0.5*(std::log((1.0+b)/(1.0-b))));
  rho = std::pow(b, 5.0);
  double k = 1.0/std::sqrt(1.0/(2.0*q*q) + rho*(x+L/2.0));
  double theta = std::asin(q/k);
  assert (rho > 0.0);

  double u = q*std::cos(theta);
  if (y < 0)
    u = -u;
  double v = -q*std::sin(theta);
  double p = (1.0/gamma)*std::pow(b, (2.0*gamma/(gamma-1.0)));
  double E = p/(gamma-1.0) + rho*(std::pow(q, 2.0))/2.0;

  vec[0] = rho; 
  vec[1] = rho*u; 
  vec[2] = rho*v; 
  vec[3] = E; 
  
  return vec;
}


std::vector<double> FIELDS::VORTEX_FIELD_MAPPING (const Node& n)
{
  std::vector<double> vec(4, 0.0);
  
  // Coordinates
  double x = n.coords[0];
  double y = n.coords[1];

  double du, dv, dT;
  double rho, T, p;
  double gamma = 1.4;
  double beta, r, u, v, q, E, R;
  double xc = 0.05;
  double yc = 0.05;

  beta = 1.0/50.0;
  R = 0.005;

  // r = std::sqrt(std::pow((x-xc), 2.0) + std::pow((y-yc), 2.0))/R;
  // u = Ghost::U*(1.0-beta*(y-yc)*(std::exp(-std::pow(r, 2.0)/2.0))/R);
  // v = Ghost::U*(beta*(x-xc)*(std::exp(-std::pow(r, 2.0)/2.0))/R);
  // T = Ghost::T - std::pow(Ghost::U*beta, 2.0)*(std::exp(-std::pow(r, 2.0)))/(2.0*Ghost::Cp);

  // rho = Ghost::rho*std::pow((T/Ghost::T), 1.0/(gamma-1.0));
  // p = rho*Ghost::R*T;
  
  double alpha = 0.0;
  double sigma = 1.0;
  auto f = (-0.5/std::pow(sigma, 2.0))*(std::pow((x-xc)/R, 2.0) + std::pow((y-yc)/R, 2.0));
  auto Omega = beta*std::exp(f);
  // auto a_inf = std::sqrt(gamma*Ghost::p/std::pow(Ghost::T,  1.0/(gamma-1.0)));
  du = -(y-yc)*Omega/R;
  dv = +(x-xc)*Omega/R;
  dT = -(gamma-1.0)*std::pow(Omega, 2.0)/2.0;
  rho = std::pow(1.0 + dT, 1.0/(gamma-1.0));
  u = Ghost::Mach*1.0 + du;
  v = dv;
  p = ((1.0/gamma)*std::pow(1.0 + dT, gamma/(gamma-1.0)));
  //p = (rho*std::pow(1.0 + dT, 1.0/(gamma-1.0)))/gamma;
  q = std::sqrt(std::pow(u, 2.0)+std::pow(v, 2.0));
  E = p/(gamma-1.0) + rho*(std::pow(q, 2.0))/2.0;

  vec[0] = rho; 
  vec[1] = rho*u; 
  vec[2] = rho*v; 
  vec[3] = E; 

  return vec;
};

