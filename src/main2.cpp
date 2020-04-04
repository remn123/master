#include <iostream>
#include <filesystem>
#include <string>
#include <map>
#include <memory>
//#include <crtdbg.h>

#include <Mesh.h>
#include <Solver.h>

namespace fs = std::filesystem;

int main()
{
  std::shared_ptr<Mesh> mesh = std::make_shared(Mesh{2});
  
  int order = 2;
  std::shared_ptr<SD<Euler>> sd = std::make_shared(SD<Euler>{order});

  fs::path cur_path = fs::current_path();

  mesh->read_gmsh((cur_path / ".." / "resources" / "mesh1.msh").string());
  
  // 1) Setup solver
  //    Calculates all flux and solution nodes and initial values
  sd->create_nodes();
  
  for (auto& e : mesh->elems)
  {
    sd->initialize_properties(e, mesh->nodes[e->nodes]);
  }
  // 2) Main loop
  while (1)
  {
    for (auto& e : mesh->elems)
    {
      sd->interpolate_sp2fp(e);
      sd->calculate_fluxes(e);
      sd->bondary_condition(e);
      
      sd->riemann_solver(e);
      sd->interpolate_fp2sp(e);
      sd->residue(e);
    }

    for (auto& e : mesh->elems)
    {
      sd.solve(e);
    }
  }


  std::cin.get();

  return 1;
}
