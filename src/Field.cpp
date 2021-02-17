#include <Field.h>


/* Declaring unitary field */
std::vector<double> FIELDS::DEFAULT_FIELD_MAPPING (const Node& n)
{
  std::vector<double> vec(4, 1.0);
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

  double pho = 1.0;
  double u = (1.0/(sigmax*sqrt(2.0*M_PI)))*exp(-(x-mux)*(x-mux)/(2.0*sigmax));
  double v = (1.0/(sigmay*sqrt(2.0*M_PI)))*exp(-(y-muy)*(y-muy)/(2.0*sigmay));
  double E = 1.0;

  vec[0] = pho; 
  vec[1] = pho*u; 
  vec[2] = pho*v; 
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

  double pho = 1.0;
  double u = ax*x + bx;
  double v = ay*y + by;
  double E = 1.0;

  vec[0] = pho; 
  vec[1] = pho*u; 
  vec[2] = pho*v; 
  vec[3] = E; 
  
  return vec;
};
