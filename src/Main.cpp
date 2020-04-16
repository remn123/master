#include <iostream>
#include <filesystem>
#include <string>
#include <map>
//#include <crtdbg.h>

#include <Mesh.h>

namespace fs = std::filesystem;

int main()
{
  Static_Mesh mesh1 {2};
  Mesh mesh2 {2};
  //Mesh mesh1 = Mesh(2);
  fs::path cur_path = fs::current_path();

  mesh1.read_gmsh((cur_path / ".." / "resources" / "mesh1.msh").string());
  mesh2.read_gmsh((cur_path / ".." / "resources" / "mesh2.msh").string());

  mesh1.print_element_by_id(0);
  mesh1.print_element_by_id(1);
  //mesh1.print_element_by_id(80);
  // Creating kd-tree
  mesh1.createKDtree();
  mesh1.to_graphviz();
	
  // Instantiating a new node
  //Node n2{ 0.4, 0.3, 0.0 };

  std::cout << "Marking Fringes" << std::endl;
  //long elm_id = mesh1.mark_fringes(n2);
  //std::cout << "Element " << elm_id << " is a fringe!" << std::endl;
	
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
            //std::cout << "Element " << elm_id << " is a fringe!" << std::endl;
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

  /*for (auto& elm : mesh1.elems)
  {
    for (auto& e : elm->edges)
    {
      std::cout << "Element(" << elm->id << "/" << e.left << ")->neighbor = " << e.right << std::endl;
    }
    std::cin.get();
  }*/

  /*for (auto& elm : mesh1.elems)
  {
    for (auto& e : elm->edges)
    {
      for (auto& v : e.nodes)
      {
	std::cout << v << " ";
      }
      std::cout << " -> " << e.id << std::endl;
    }
  }*/
  /*for (auto& e : mesh1.elems[0]->edges_map)
  {
    for (auto& v : e.first)
    {
      std::cout << v << " ";
    }
    std::cout << " -> " << e.second << std::endl;
  }*/


  /*std::cout << "Elements from each node " << std::endl;
  for (auto& n : mesh1.nodes)
  {
    std::cout << "Node(" << n.id << ") is on elements: ";
    for (auto& e : n.elems)
    {
      std::cout << e << " ";
    }
    std::cout << std::endl;
  }*/

  std::cin.get();

  return 1;
}
