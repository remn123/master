#include <iostream>
#include <experimental/filesystem>
#include <string>
#include <map>
//#include <crtdbg.h>

#include <Mesh.h>
#include <Solver.h>

namespace fs = std::experimental::filesystem;

int main()
{
  Mesh mesh {2};
  
  int order = 2;
  SD sd {order};

  fs::path cur_path = fs::current_path();

  mesh.read_gmsh((cur_path / ".." / "resources" / "mesh1.msh").string());
  
  // 1) Setup solver
  //    Calculates all flux and solution nodes
  sd.setup();

  // 2) Main loop
  while (1)
  {
    for (auto& e : mesh.elems)
    {
      sd.solve(e);
    }
  }


  std::cin.get();

  return 1;
}
