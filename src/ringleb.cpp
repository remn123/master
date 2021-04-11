#include <iostream>
#include <filesystem>
#include <functional>
#include <string>
#include <map>
#include <memory>
//#include <crtdbg.h>

#include <Dummies.h>
#include <Mesh.h>
#include <SD.h>
//#include <Solver.h>
#include <Time.h>

namespace fs = std::filesystem;

// explicit instantiations
// template class SD<Euler>;
// template class SD<NavierStokes>;

// template class Time<Explicit::SSPRungeKutta>;

/* Declaring analytical field solution */

int main()
{
  fs::path cur_path = fs::current_path();
  auto mesh = std::make_shared<Mesh>(2);

  mesh->read_gmsh((cur_path.parent_path() / "resources" / "ringleb_v0.msh").string());
  

  int order = 3;
  auto sd = std::make_shared<SD<Euler>>(order, 2);
  /*
    1) Setup (all element in Mesh)
      1.1) Calculate solution and fluxes points
      1.2) Initialize solution and fluxes
  */
  //auto field = [](const Node& n){return FIELDS::RINGLEB_FIELD_MAPPING(n);};
  Ghost::analytical_solution = FIELDS::RINGLEB_FIELD_MAPPING;

  sd->setup(
    mesh, 
    FIELDS::RINGLEB_FIELD_MAPPING
  );

  auto filename = (cur_path.parent_path() / "results" / "animation" / "pp_mesh_ringleb_v0_").string();
  std::string tstamp = std::to_string(0);
  tstamp.insert(tstamp.begin(), 5 - tstamp.length(), '0');
  sd->to_vtk(mesh, filename + tstamp + std::string{".vtk"});
  

  // std::cout << "this->Ngh = " << mesh->Ngh << "\n";
  // for (auto &g : mesh->ghosts)
  // {
  //   std::cout << "g.id = " << g.id << "\n";
  //   std::cout << "g.edg_id = " << g.edg_id << "\n";
  //   std::cout << "g.elm_id = " << g.elm_id << "\n";
  //   std::cout << "g.group = " << g.group << "\n";
  //   std::cout << "g.lr_edge = " << g.lr_edge << "\n";
  //   std::cout << "g.type = " << g.type << "\n";
  //   for (auto &f : g.fnodes)
  //   {
  //     std::cout << "f[" << f.id << "].coords = (" << f.coords[0] << ", " << f.coords[1] << ")\n";  
  //     std::cout << "f[" << f.id << "].local = " << f.local << "\n";  
  //     std::cout << "f[" << f.id << "].right = " << f.right << "\n\n";  
  //   }
  // }
  // std::cin.get();

  /*
    2) Solver Loop (for each element in Mesh)
       2.1) Apply Boundary Conditions on interfaces (flux points)
       2.2) Interpolate Q and dQ from SPs to FPs
       2.3) Calculate Fluxes at internal FPs
       2.4) Calculate Fluxes at faces FPs (riemann solver and BR2)
       2.5) Interpolate Fluxes derivatives from FPs to SPs
       2.6) Calculate Residue
  */

  sd->solve(mesh);

  /*
    3) Time Marching Loop
      3.1) Calculate residue norm
      3.2) Check if it's already converged
      3.3) (if not) Apply time iteration then go to (2)
  */
  double CFL = 0.01;
  long MAX_ITER = 1E+4;
  int rk_order = 4;
  int stages = 5;
  int size = mesh->Nel * (order * order)*4; // overall number of solution points

  auto time = std::make_shared<Time<Explicit::SSPRungeKutta>>(CFL,
                                                              MAX_ITER,
                                                              stages,
                                                              rk_order,
                                                              size);
  
  
  time->loop(
    mesh, 
    [&sd](std::shared_ptr<Mesh> & m){sd->solve(m);},
    filename,
    [&sd](const std::shared_ptr<Mesh> & m, const std::string & f){sd->to_vtk(m, f);}
  );

  
  time->save(
    mesh, 
    filename, 
    [&sd](const std::shared_ptr<Mesh> & m, const std::string & f){sd->to_vtk(m, f);}
  );

  return 0;
}
