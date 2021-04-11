#include <vector>
#include <math.h>  

#include <Node.h>

namespace FIELDS {
  std::vector<double> DEFAULT_FIELD_MAPPING (const Node&);
  std::vector<double> LINEAR_FIELD_MAPPING (const Node&);
  std::vector<double> GAUSSIAN_FIELD_MAPPING (const Node&);
  std::vector<double> RINGLEB_FIELD_MAPPING (const Node&);
}
