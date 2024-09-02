#include <iostream>
#include <filesystem>
#include <string>
#include <map>

#include <Mesh.h>

namespace fs = std::filesystem;

int main()
{
  Static_Mesh mesh1 {2};
  Mesh mesh2 {2};
  
  fs::path cur_path = fs::current_path();

  mesh1.read_gmsh((cur_path / ".." / "resources" / "mesh1.msh").string());
  mesh2.read_gmsh((cur_path / ".." / "resources" / "mesh2.msh").string());

  mesh1.print_element_by_id(0);
  mesh1.print_element_by_id(1);
  
  // Creating kd-tree
  mesh1.create_kdtree();
  mesh1.to_graphviz();
	
  std::cout << "Marking Fringes" << std::endl;
  
  long elm_id;

  for (auto& e : mesh2.elems)
  {
    if(e->boundary==1)
    {
      for (auto& ed : e->edges)
      {
        if(ed.boundary==1)
        {
          for (auto& n : ed.nodes)
          {
            elm_id = mesh1.mark_fringes(mesh2.nodes[n]);
          }
        }
      }
    }
  }
	
  mesh1.to_vtk((cur_path / ".." / "results" / "mesh3.vtk").string());
  mesh2.to_vtk((cur_path / ".." / "results" / "mesh2.vtk").string());

  std::cout << " PRINTING KDTREE " << '\n';
  long height = 0;
  std::string before = {""};
  int flag_right=0;
  mesh1.print_tree(mesh1.root, height, before, flag_right);

  std::cin.get();

  return 1;
}
