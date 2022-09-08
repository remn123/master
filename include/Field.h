#include <vector>
#include <math.h>  

#include <Node.h>

namespace FIELDS {
  std::vector<double> DEFAULT_FIELD_MAPPING (const Node&);
  std::vector<double> LINEAR_FIELD_MAPPING (const Node&);
  std::vector<double> GAUSSIAN_FIELD_MAPPING (const Node&);
  std::vector<double> CYLINDER_FIELD_MAPPING (const Node&);
  std::vector<double> IMPLOSION_FIELD_MAPPING (const Node&);
  std::vector<double> RINGLEB_FIELD_MAPPING (const Node&);
  std::vector<double> VORTEX_FIELD_MAPPING (const Node&);
  std::vector<double> VORTEX_T_FIELD_MAPPING (const Node&);
}
