#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <sstream>
#include <sys/stat.h>
#include <iterator>
#include <algorithm>
#include <functional>
#include <random>
#include <limits>

#include <Mesh.h>
#include <Ghost.h>

int Static_Mesh::dim = 0;
int Mesh::number_meshes = 0;
long Element::num_elems = -1;
long Element::num_faces = -1;
long Element::num_edges = -1;
std::unordered_map<std::vector<long>, std::vector<long>> Element::faces_map{};
std::unordered_map<std::vector<long>, std::vector<long>> Element::edges_map{};
//long Vertice::num_nodes = -1;

//Mesh
Mesh::Mesh(int d)
{
  number_meshes++;
  this->id = number_meshes;
  this->dimension = d;
  std::cout << "Mesh has been created!" << std::endl;
}

Mesh::~Mesh(void)
{
  std::cout << "Mesh has been deleted!" << std::endl;
}

// enumerate gmsh blocks types
enum blocks
{
  MESH_FORMAT = 1,
  eMESH_FORMAT,
  NODES,
  eNODES,
  ELEMENTS,
  eELEMENTS,
  PHYSICALNAMES,
  ePHYSICALNAMES,
  NODEDATA,
  eNODEDATA,
  ELEMENTDATA,
  eELEMENTDATA,
  ELEMENTNODEDATA,
  eELEMENTNODEDATA
};

const int blocks_arr[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

std::ostream &operator<<(std::ostream &out, const blocks block)
{
  return std::cout << blocks_arr[block];
}

void Mesh::update_element_neighbors(void)
{
  // scan all elements
  for (auto &e : this->elems)
  {
    // scan all elements' edges
    for (auto &ed1 : e->edges)
    {
      // if there is a neighbor, update neighbor's edges
      if (ed1.right >= 0)
      {
        for (auto &ed2 : this->elems[ed1.right]->edges)
        {
          if (ed2.id == ed1.id && ed2.right < 0 && e->id != this->elems[ed1.right]->id)
          {
            ed2.right = ed1.left;
          }
        }
      }
    }
  }
}

//void Mesh::calculate_jacobians(std::shared_ptr<Element>& e)
//{
//  std::vector<Node> enodes;
//  enodes.reserve(e->nodes.size());
//
//  for (auto& n_id : e->nodes)
//  {
//    enodes.push_back(this->nodes[n_id]);
//  }
//  e->calculate_jacobian(enodes);
//  enodes.clear();
//}

void Mesh::mark_boundaries(void)
{
  int local_id = 0;
  // scan all elements
  for (auto &e : this->elems)
  {
    e->boundary = 0;
    local_id = 0;
    // scan all elements' edges
    for (auto &ed : e->edges)
    {
      // if there is a negative neighbor, mark element as boundary
      if (ed.right < 0)
      {
        ed.boundary = 1;
        e->boundary = 1;
        this->ghosts.emplace_back(Ghost{e->id, ed.id, local_id, ed.type, ed.group});
        ed.ghost = Ghost::num_ghosts;
        this->Ngh = Ghost::num_ghosts;
        //break;
      }
      local_id++;
    }
  }
}

void Mesh::to_vtk(const std::string &filename)
{
  std::ofstream output;

  output.open(filename);
  output << "# vtk DataFile Version 3.0" << std::endl;
  output << "My Mesh" << std::endl;
  output << "ASCII" << std::endl;
  output << "DATASET UNSTRUCTURED_GRID" << std::endl;
  output << "POINTS " << this->N << " float" << std::endl;
  //output << "POINT_DATA " << this->N << std::endl;

  // writing nodes
  for (auto &n : this->nodes)
  {
    //output << n.id+1 << " " << n.coords[0] << " " << n.coords[1] << " " << n.coords[2];
    output << n.coords[0] << " " << n.coords[1] << " " << n.coords[2];
    output << std::endl;
  }
  output << std::endl;

  auto nodes_size = 0;

  for (auto &e : this->elems)
  {
    nodes_size += 4;
  }

  output << "CELLS " << this->Nel << " " << this->Nel + nodes_size << std::endl;

  // writing elements
  for (auto &e : this->elems)
  {
    output << e->nodes.size();
    for (auto &enode : e->nodes)
    {
      output << " " << enode;
    }
    output << std::endl;
  }
  output << std::endl;

  output << "CELL_TYPES " << this->Nel << std::endl;

  // writing elements type
  for (auto &e : this->elems)
  {
    // temporary "fix"
    output << 9 << std::endl;
  }
  output << std::endl;

  output << "CELL_DATA " << this->Nel << std::endl;
  output << "SCALARS boundary int 1" << std::endl;
  output << "LOOKUP_TABLE default" << std::endl;

  // writing element kdtree boundary flag
  for (auto &e : this->elems)
  {
    output << e->boundary << std::endl;
  }

  output << "SCALARS fringe int 1" << std::endl;
  output << "LOOKUP_TABLE default" << std::endl;

  // writing element kdtree boundary flag
  for (auto &e : this->elems)
  {
    output << e->fringe << std::endl;
  }
  output.close();
}

void Mesh::read_gmsh(const std::string &filename)
{
  struct stat buffer;
  // check if file exists
  if (stat(filename.c_str(), &buffer) == 0)
  {
    std::cout << "File <" + filename + "> was found!" << std::endl;

    std::ifstream gmshFile(filename);
    std::string line;

    blocks _BLOCK_ = blocks::MESH_FORMAT;

    std::map<std::string, blocks> blck_map = {
        {"$MeshFormat", blocks::MESH_FORMAT},
        {"$EndMeshFormat", blocks::eMESH_FORMAT},
        {"$Nodes", blocks::NODES},
        {"$EndNodes", blocks::eNODES},
        {"$Elements", blocks::ELEMENTS},
        {"$EndElements", blocks::eELEMENTS},
        {"$NodeData", blocks::NODEDATA},
        {"$EndNodeData", blocks::eNODEDATA},
        {"$ElementData", blocks::ELEMENTDATA},
        {"$EndElementData", blocks::eELEMENTDATA},
        {"$ElementNodeData", blocks::ELEMENTNODEDATA},
        {"$EndElementNodeData", blocks::eELEMENTNODEDATA}};

    while (!gmshFile.eof())
    {
      std::getline(gmshFile, line);

      auto search = blck_map.find(line);
      //auto search2 = blck_map.find("$MeshFormat\n");

      // std::cout << "line.length = " << line.length() << std::endl;
      // std::cout << "line= " << line << std::endl;
      // for (auto &a : line)
      //   std::cout << a << "\n";
      // std::cout << "FIM\n";
      // std::cout << "search->first.length() = " << search->first.length() << std::endl;
      //std::cout << "search2->first = " << search2->first << " search2->first.length() " << search2->first.length() << std::endl;
      if (search != blck_map.end())
      {
        //std::cout << "Found " << search->first << " " << search->second << std::endl;
        _BLOCK_ = search->second;
      }
      // std::cin.get();
      std::istringstream check1(line);

      //std::vector<std::string> parser;
      std::string value;

      // Getting into the Node Block
      if (_BLOCK_ == blocks::NODES)
      {
        std::vector<std::string> parser(std::istream_iterator<std::string>{check1},
                                        std::istream_iterator<std::string>());

        if (parser.size() > 2)
        {
          // std::cout << "[Nodes]: " << line << std::endl;
          // std::cout << "That's the parser, dude: " << std::endl;
          // for (const auto &v : parser)
          //   std::cout << v << ", ";
          // std::cout << std::endl;

          // [0] node_ID
          // [1] x
          // [2] y
          // [3] z

          //long node_ID = std::stol(parser[0]);
          auto first = parser.begin() + 1;
          auto last = parser.end();

          std::vector<std::string> newVec{first, last};
          // std::cout << "Creating mesh vertice..." << std::endl;
          // for (auto &n : newVec)
          //   std::cout << n << " ";
          // std::cout << std::endl;
          // std::cout << "N Nodes = " << this->N << " - this->nodes.size = " << this->nodes.size() << std::endl;
          this->nodes.emplace_back(Vertice{newVec});
          //std::cout << "Created mesh vertice!" << std::endl;

          // for (auto &n : this->nodes)
          //   std::cout << n.id << " ";
          // std::cout << std::endl;
          //newVec.clear();
          //parser.clear();
        }
        else if (line.find("$Nodes") == std::string::npos &&
                 line.find("$EndNodes") == std::string::npos)
        {
          //std::cout << "Reserving N Nodes..." << std::endl;
          this->N = std::stol(parser[0]);
          this->nodes.reserve(this->N);
          //std::cout << "N Nodes = " << this->N << " - this->nodes.size = " << this->nodes.size() << std::endl;
          // this->nodes.resize(this->N);
        }
      }
      // Getting into the Element Block
      else if (_BLOCK_ == blocks::ELEMENTS)
      {
        std::vector<std::string> parser(std::istream_iterator<std::string>{check1},
                                        std::istream_iterator<std::string>());

        if (parser.size() > 2)
        {
          //std::cout << "[Elements]: " << line << std::endl;
          //std::cout << "That's the parser, dude: " << std::endl;
          //for (const auto& v : parser)
          //	std::cout << v << ", ";
          //std::cout << std::endl;

          // [0] elem_ID
          // [1] elem_type
          // [2] number_of_tags (k)
          // [3..2+k] tags(k)
          // [3+k...vector.size()-1] nodes

          //long elem_ID = std::stol(parser[0]);
          int elem_type = std::stoi(parser[1]);
          int num_tags = std::stoi(parser[2]);
          int physical_tag = std::stoi(parser[3]);
          int elementary_tag = std::stoi(parser[4]);
          auto first = parser.begin() + 2 + num_tags + 1;
          auto last = parser.end();
          std::vector<std::string> newVec(first, last);

          if (elm_type(elem_type) != elm_type::NODE2_LINE &&
              elm_type(elem_type) != elm_type::NODE3_O2_LINE &&
              elm_type(elem_type) != elm_type::NODE1_POINT &&
              elm_type(elem_type) != elm_type::NODE4_O3_EDG &&
              elm_type(elem_type) != elm_type::NODE5_O4_EDG &&
              elm_type(elem_type) != elm_type::NODE6_O5_EDG)
          {
            switch (elem_type)
            {
            case 2:
              elems.emplace_back((std::make_shared<Triangle>(newVec)));
              break;
            case 3:
              elems.emplace_back((std::make_shared<Quadrangle>(newVec)));
              break;
            case 4:
              elems.emplace_back((std::make_shared<Tetrahedron>(newVec)));
              break;
            case 5:
              elems.emplace_back((std::make_shared<Hexahedron>(newVec)));
              break;
            case 6:
              elems.emplace_back((std::make_shared<Prism>(newVec)));
              break;
            case 7:
              elems.emplace_back((std::make_shared<Pyramid>(newVec)));
              break;
            case 9:
              elems.emplace_back((std::make_shared<Triangle>(newVec)));
              break;
            case 10:
              elems.emplace_back((std::make_shared<Quadrangle>(newVec)));
              break;
            case 11:
              elems.emplace_back((std::make_shared<Tetrahedron>(newVec)));
              break;
            case 12:
              elems.emplace_back((std::make_shared<Hexahedron>(newVec)));
              break;
            case 13:
              elems.emplace_back((std::make_shared<Prism>(newVec)));
              break;
            case 14:
              elems.emplace_back((std::make_shared<Pyramid>(newVec)));
              break;
            case 16:
              elems.emplace_back((std::make_shared<Quadrangle>(newVec)));
              break;
            case 17:
              elems.emplace_back((std::make_shared<Hexahedron>(newVec)));
              break;
            case 18:
              elems.emplace_back((std::make_shared<Prism>(newVec)));
              break;
            case 19:
              elems.emplace_back((std::make_shared<Pyramid>(newVec)));
              break;
            case 20:
              elems.emplace_back((std::make_shared<Triangle>(newVec)));
              break;
            case 21:
              elems.emplace_back((std::make_shared<Triangle>(newVec)));
              break;
            case 22:
              elems.emplace_back((std::make_shared<Triangle>(newVec)));
              break;
            case 23:
              elems.emplace_back((std::make_shared<Triangle>(newVec)));
              break;
            case 24:
              elems.emplace_back((std::make_shared<Triangle>(newVec)));
              break;
            case 25:
              elems.emplace_back((std::make_shared<Triangle>(newVec)));
              break;
            case 29:
              elems.emplace_back((std::make_shared<Tetrahedron>(newVec)));
              break;
            case 30:
              elems.emplace_back((std::make_shared<Tetrahedron>(newVec)));
              break;
            case 31:
              elems.emplace_back((std::make_shared<Tetrahedron>(newVec)));
              break;
            }
            this->append_elem_to_nodes(elems.back());
            newVec.clear();
          }
          else if (elm_type(elem_type) == elm_type::NODE2_LINE)
          {
            this->append_boundary_face(physical_tag, newVec);
          }

          parser.clear();
        }
        else if (line.find("$Elements") == std::string::npos &&
                 line.find("$EndElements") == std::string::npos)
        {
          this->Nel = std::stol(parser[0]);
          this->elems.reserve(this->Nel);
        }
      }
    }
    gmshFile.close();
    this->Nel = Element::num_elems + 1;
    this->Ned = Element::num_edges + 1;
    this->update_element_neighbors();
    this->mark_boundaries();
    this->Ngh = Ghost::num_ghosts;
    // reset num_elems, num_edges and num_ghosts
    Element::num_elems = -1;
    Element::num_faces = -1;
    Element::num_edges = -1;
    Ghost::num_ghosts = 0;
    Vertice::num_nodes = -1;
    Element::edges_map.clear();

    std::cout << "File <" + filename + "> has been read." << std::endl;
  }
  else
  {
    std::cout << "File <" + filename + "> not found!" << std::endl;
  }
}

int Mesh::get_id(void)
{
  return id;
}

int Mesh::get_dimension(void)
{
  return dimension;
}

long Mesh::get_number_nodes(void)
{
  return N;
}

long Mesh::get_number_elements(void)
{
  return Nel;
}

void Mesh::print_node_id(long id)
{
  std::cout << "Vertice(" << id << ").id = " << this->nodes[id].id << std::endl;
}

void Mesh::print_element_by_id(long elem_id)
{
  this->elems[elem_id]->print_nodes();
}

void Mesh::print_node_by_id(long node_id)
{
  this->nodes[node_id].print_coords();
}

double Mesh::get_area(const std::vector<Vertice> &vec)
{
  size_t size = 0;
  double crossp = 0.0;
  double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;

  while (size < vec.size())
  {
    if (size < vec.size() - 1)
    {
      x1 = vec[size].coords[0];
      x2 = vec[size + 1].coords[0];
      y1 = vec[size].coords[1];
      y2 = vec[size + 1].coords[1];
    }
    else
    {
      x1 = vec[size].coords[0];
      x2 = vec[0].coords[0];
      y1 = vec[size].coords[1];
      y2 = vec[0].coords[1];
    }
    crossp += x1 * y2 - y1 * x2;
    size++;
  }
  std::cout << "Area = " << crossp << std::endl;
  return crossp;
}

double Mesh::get_area(const std::vector<long> &vec, const Vertice &n2)
{
  size_t size = 0;
  double crossp = 0.0;
  double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;

  while (size <= vec.size())
  {
    if (size < vec.size() - 1)
    {
      x1 = this->nodes[vec[size]].coords[0];
      x2 = this->nodes[vec[size + 1]].coords[0];
      y1 = this->nodes[vec[size]].coords[1];
      y2 = this->nodes[vec[size + 1]].coords[1];
    }
    else if (size == vec.size() - 1)
    {
      x1 = this->nodes[vec[size]].coords[0];
      x2 = n2.coords[0];
      y1 = this->nodes[vec[size]].coords[1];
      y2 = n2.coords[1];
    }
    else if (size == vec.size())
    {
      x1 = n2.coords[0];
      x2 = this->nodes[vec[0]].coords[0];
      y1 = n2.coords[1];
      y2 = this->nodes[vec[0]].coords[1];
    }
    crossp += x1 * y2 - y1 * x2;
    size++;
  }

  return crossp;
}

double Mesh::get_volume(const std::vector<Vertice> &vec)
{
  // pass
  return vec[0].coords[0];
}

void Mesh::append_elem_to_nodes(const std::shared_ptr<Element> &e)
{
  /* append the element id to all nodes that this belongs to this element */

  for (auto &n_id : e->nodes)
  {
    //std::cout << "Appending " << e->id << " to node " << n_id << std::endl;
    this->nodes[n_id].elems.push_back(e->id);
  }
  //std::cout << "Appending done!" << std::endl;
}

void Mesh::append_boundary_face(const int physical_tag, const std::vector<std::string> &node_list)
{
  std::vector<long> edges_nodes = {};
  edges_nodes.reserve(2);

  // transforming vector<string> to vector<int>
  // subtract 1 from the node_id because the first id for GMSH is 1, but as I am using the node id to consult it on a vector
  // I need it to be 0-started.
  std::transform(node_list.begin(), node_list.end(), std::back_inserter(edges_nodes),
                 [](const std::string &str) { return std::stol(str) - 1; });

  // Sort vector so when we search for this edge, we have an easy match
  std::sort(edges_nodes.begin(), edges_nodes.end());

  // Add new boundary edge to the map
  this->bc_map.emplace(std::make_pair(edges_nodes, physical_tag));
}

// ------------------------------------------------------------------------ //

// Static Mesh
Static_Mesh::Static_Mesh(int d) : Mesh(d)
{
  std::cout << "Static_Mesh has been created!" << std::endl;
}

Static_Mesh::~Static_Mesh(void)
{
  std::cout << "Static_Mesh with root at node(" << this->nodes[this->root].id << ")has been deleted!" << std::endl;
}

void Static_Mesh::create_kdtree(void)
{
  std::cout << "Start KD-Tree creation..." << std::endl;
  std::vector<long> nodes_ptr;
  long k = -1;
  //std::vector<Node *> nodes_ptr;

  for (auto &n : this->nodes)
  {
    //std::cout << "Node(" << n.id << "): " << " n->coords ={" << n.coords[0] << " " << n.coords[1] << " " << n.coords[2] << std::endl;
    //std::cout << "Node Address: " << &n << std::endl;
    nodes_ptr.emplace_back(n.id);
  }

  //std::cout << "1) Size of v_nodes = " << nodes_ptr.size() << " and Size of nodes = " << this->nodes.size() << std::endl;
  this->get_pivot(this->root, nodes_ptr);
  //std::cout << " Root has been assigned" << std::endl;
  this->build_kdtree(this->root, nodes_ptr, k);

  std::cout << "KD-Tree has been created with root at " << this->nodes[this->root].id << std::endl;
  nodes_ptr.clear();
}

//void Static_Mesh::print_tree(void)
//{
//
//}

// Status
Status::~Status(void)
{
}

void Static_Mesh::build_kdtree(long &root_, std::vector<long> &v_nodes, long k)
//void Static_Mesh::build_kdtree(Node* root_, std::vector<Node*>& v_nodes)
{
  /*
    First I will get a random sample of my node vector.
    Then I order this sample by the layer dimension (x, y or z)
    get the median and return it as pivot.
  */
  //std::cout << "Building KDtree..." << std::endl;

  //for (auto& n : v_nodes)
  //{
  //	std::cout << "Build -> Node(" << this->nodes[n].id << "): " << " n->coords = {" << this->nodes[n].coords[0] << " " << this->nodes[n].coords[1] << " " << this->nodes[n].coords[2] << "}" << std::endl;
  //}

  k++;
  k = k % this->dimension;

  if (v_nodes.size() > 1)
  {
    std::vector<long> left_nodes;
    std::vector<long> right_nodes;
    //std::vector<Node *> left_nodes;
    //std::vector<Node *> right_nodes;

    // Split Nodes which have the coordinate dim less or greater than the root (pivot)
    for (auto &n : v_nodes)
    {
      if (n != root_)
      {
        if (this->nodes[n].coords[k] > this->nodes[root_].coords[k])
        {
          right_nodes.emplace_back(n);
        }
        //if (this->nodes[n].coords[dim] <= this->nodes[root_].coords[dim])
        else
        {
          left_nodes.emplace_back(n);
        }
      }
    }

    //dim++;
    //dim = dim % this->dimension; // cyclic layer dim
    //std::cout << "Dim = " << dim << std::endl;
    dim = (k + 1) % this->dimension;
    // Find out left_nodes pivot and traverse it
    this->get_pivot(this->nodes[root_].left, left_nodes);
    if (this->nodes[root_].left > -1)
    {
      this->build_kdtree(this->nodes[root_].left, left_nodes, k);
      //std::cout << "Node(" << this->nodes[root_].id << ")[" << this->nodes[root_].coords[dim] << "]->left = " << "Node(" << this->nodes[this->nodes[root_].left].id << ")[" << this->nodes[this->nodes[root_].left].coords[dim] << "]" << std::endl;
    }
    left_nodes.clear();

    // Find out right_nodes pivot and traverse it
    this->get_pivot(this->nodes[root_].right, right_nodes);
    if (this->nodes[root_].right > -1)
    {
      this->build_kdtree(this->nodes[root_].right, right_nodes, k);
      //std::cout << "Node(" << this->nodes[root_].id << ")[" << this->nodes[root_].coords[dim] << "]->right = " << "Node(" << this->nodes[this->nodes[root_].right].id << ")[" << this->nodes[this->nodes[root_].right].coords[dim] << "]" << std::endl;
    }
    right_nodes.clear();
  }
}

void Static_Mesh::get_pivot(long &root_, const std::vector<long> &v_nodes)
{
  /*
    1) Sample v_nodes
    2) Sort the sampling
    3) Get the median
    4) Return median as pivot
  */

  //std::cout << "Size of v_nodes = " << v_nodes.size() << std::endl;

  if (v_nodes.size() == 0)
  {
    root_ = -1;
  }
  else
  {
    std::vector<long> sampled_nodes;
    long k = static_cast<long>(v_nodes.size() * 5 / 100); // 5% of population

    if (k < 2)
    {
      k = static_cast<long>(v_nodes.size());
    }

    std::sample(v_nodes.begin(), v_nodes.end(), std::back_inserter(sampled_nodes), k, std::mt19937{std::random_device{}()});

    size_t size = sampled_nodes.size();

    std::sort(sampled_nodes.begin(), sampled_nodes.end(), [this](long n1, long n2) {
      return this->nodes[n1].coords[dim] < this->nodes[n2].coords[dim];
    });

    if (size % 2 == 0)
    {
      //std::cout << "Dim = " << dim << "; Median(" << this->nodes[sampled_nodes[size / 2]].id << ") -> " << this->nodes[sampled_nodes[size / 2]].coords[dim] << std::endl;
      root_ = sampled_nodes[size / 2];
    }
    //else if (size == i1)
    //{
    //std::cout << "Dim = " << dim << "; Median(" << this->nodes[sampled_nodes[0]].id << ") -> " << this->nodes[sampled_nodes[0]].coords[dim] << std::endl;
    //  root_ = sampled_nodes[0];
    //}
    else
    {
      //std::cout << "Dim = " << dim << "; Median(" << this->nodes[sampled_nodes[(size + 1) / 2]].id << ") -> " << this->nodes[sampled_nodes[(size + 1) / 2]].coords[dim] << std::endl;
      root_ = sampled_nodes[(size + 1) / 2 - 1];
    }
    sampled_nodes.clear();
  }
}

long Static_Mesh::mark_fringes(const Vertice &n2)
{
  int k = -1;
  double min_dist = std::numeric_limits<double>::max();
  //std::cout << min_dist << std::endl;
  long closest = -1;
  closest = this->_search(n2, this->root, k, min_dist, this->root);
  std::cout << "Test Node: " << n2.id << std::endl;
  std::cout << "Closest: " << closest << std::endl;
  long find_elm = this->_find_element(n2, this->nodes[closest].elems);
  std::cout << "Find fringe at " << find_elm << std::endl;
  return find_elm;
}

void Static_Mesh::_get_kneighbors(long &sroot, std::vector<long> &kneighbors)
{
  if (sroot != -1)
  {
    //std::cout << "kNEIGHBORS(" << &kneighbors << ")\n";
    kneighbors.push_back(sroot);
    this->_get_kneighbors(this->nodes[sroot].left, kneighbors);
    this->_get_kneighbors(this->nodes[sroot].right, kneighbors);
  }
}

long Static_Mesh::_get_closest(const Vertice &n2, std::vector<long> &kneighbors)
{
  double min_dist = -1.0, dist = 0.0;
  double x1, y1, z1;
  double x2, y2, z2;
  long closest;

  x2 = n2.coords[0];
  y2 = n2.coords[1];
  z2 = n2.coords[2];

  for (auto &n : kneighbors)
  {
    x1 = this->nodes[n].coords[0];
    y1 = this->nodes[n].coords[1];
    z1 = this->nodes[n].coords[2];

    dist = sqrt(pow((x1 - x2), 2.0) + pow((y1 - y2), 2.0) + pow((z1 - z2), 2.0));

    if (min_dist > dist || min_dist < 0)
    {
      min_dist = dist;
      closest = n;
    }
  }

  return closest;
}

double Static_Mesh::_get_distance(const Vertice &n1, const Vertice &n2)
{
  double x1, y1, z1;
  double x2, y2, z2;

  x1 = n1.coords[0];
  y1 = n1.coords[1];
  z1 = n1.coords[2];

  x2 = n2.coords[0];
  y2 = n2.coords[1];
  z2 = n2.coords[2];

  return sqrt(pow((x1 - x2), 2.0) + pow((y1 - y2), 2.0) + pow((z1 - z2), 2.0));
}

long Static_Mesh::_search(const Vertice &n2, long &sroot, int k, double min_dist, long closest)
{
  k++;
  k = k % this->dimension;

  double max_k_axis = 0.0, min_k_axis = 0.0;
  long new_closest = -1;

  // First, calculate the distance between the current root and test node (radius)
  double radius = this->_get_distance(this->nodes[sroot], n2);
  min_dist = this->_get_distance(this->nodes[closest], n2);

  if (radius <= min_dist)
  {
    min_dist = radius;
    new_closest = sroot;

    // Circunference Limits in the k direction
    max_k_axis = n2.coords[k] + radius;
    min_k_axis = n2.coords[k] - radius;
  }
  else
  {
    new_closest = closest;
    max_k_axis = n2.coords[k] + min_dist;
    min_k_axis = n2.coords[k] - min_dist;
  }

  if (n2.coords[k] <= this->nodes[sroot].coords[k] || min_k_axis < this->nodes[sroot].coords[k])
  {
    if (this->nodes[sroot].left >= 0)
    {
      new_closest = this->_search(n2, this->nodes[sroot].left, k, min_dist, new_closest);
      min_dist = this->_get_distance(this->nodes[new_closest], n2);
      max_k_axis = n2.coords[k] + min_dist;
    }
  }

  if (n2.coords[k] > this->nodes[sroot].coords[k] || max_k_axis > this->nodes[sroot].coords[k])
  {
    if (this->nodes[sroot].right >= 0)
    {
      new_closest = this->_search(n2, this->nodes[sroot].right, k, min_dist, new_closest);
    }
  }

  return new_closest;
}

//long Static_Mesh::_search(const Node& n2, long& sroot, int& k, long& height)
//{
//	k++;
//	k = k % this->dimension;
//
//        double max_height = 3.0;

//	height++;
//	if (n2.coords[k] <= this->nodes[sroot].coords[k])
//	{

//		if (this->nodes[sroot].left < 0 || height > (long)(log2(this->N)-max_height))
//		{
//			std::vector<long> kneighbors;
//std::cout << "FOUND NODE: " << sroot << std::endl;

//			this->_get_kneighbors(sroot, kneighbors);

//std::cout << "WITH NEIGHBORS(" << &kneighbors << "): " ;
//for (auto& n : kneighbors)
//	std::cout << n << " ";
//std::cout << std::endl;

//			long closest = this->_get_closest(n2, kneighbors);
//std::cout << "CLOSEST NODE: " << closest << std::endl;
//std::cin.get();
//			return this->_find_element(n2, this->nodes[closest].elems);
//		}
//		else return this->_search(n2, this->nodes[sroot].left, k, height);
//	}
//	else
//	{
//k++;
//k = k % this->dimension;
//		if (this->nodes[sroot].right < 0 || height > (long)(log2(this->N)-max_height))
//		{
//			std::vector<long> kneighbors;
//std::cout << "FOUND NODE: " << sroot << std::endl;

//			this->_get_kneighbors(sroot, kneighbors);

//std::cout << "WITH NEIGHBORS(" << &kneighbors << "): " ;
//for (auto& n : kneighbors)
//	std::cout << n << " ";
//std::cout << std::endl;

//			long closest = this->_get_closest(n2, kneighbors);
//std::cout << "CLOSEST NODE: " << closest << std::endl;
//std::cin.get();
//			return this->_find_element(n2, this->nodes[closest].elems);
//		}
//		else return this->_search(n2, this->nodes[sroot].right, k, height);
//	}
//	return -2;
//}

long Static_Mesh::_find_element(const Vertice &n2, std::vector<long> &elm_ids)
{
  /*
  Inside the element:
   _________
  |         |
  |    0----|-------> INF
  |_________|

  Outside the element:
   _________
  |         |
  |         | 0-----> INF
  |_________|

   _________
  |         |
  0---|---------|-------> INF
  |_________|


  Edge type I:
  node_A ---- node_B : node_A[1] & node_B[1]

  Algorithm:
  1) Check if the new node is colinear wrt the edge (or lies on the face's plane);
  2) Find out the "edge type";
  3) Validade if the "infinity line segment" may pierce the edge/face and aggregate it on the pierce counter;
  4) Check if the number of pierces is odd (inside) or even (outside) and return;
  */

  double area, min_area;
  //double volume;
  long fringed_elm_id = -1;

  // find which element the node n2 is covered by
  for (auto &e_id : elm_ids)
  {
    if (this->dimension == 2)
    {
      area = 0.0;
      min_area = 10.0;
      for (auto &ed : this->elems[e_id]->edges)
      {
        //auto n = ed.nodes;
        //auto last_node = n2;
        //n.push_back(last_node);

        area = this->get_area(ed.nodes, n2);
        //n.clear();
        min_area = (min_area >= area) ? area : min_area;
        // if the triangle formed by an edge and the n2 point has negative area,
        // it means that n2 is outside of the element.
        if (min_area <= 0.0)
        {
          break;
        }
      }

      if (min_area >= 0.0)
      {
        // if(e_id==2 || e_id==21)
        // {
        // 	std::cout << "Point(" << n2.coords[0] << "," << n2.coords[1] << ") is inside of the element (" << e_id << ")" << std::endl;
        // 	this->print_element_by_id(e_id);
        // }

        this->elems[e_id]->fringe = 1; // mark as fringed
        fringed_elm_id = e_id;
      }
    }
    else if (this->dimension == 3)
    {
      /*volume = 0.0;
      for (auto& face : this->elems[e_id]->faces)
      {
        //auto n = face.nodes;
        //auto last_node = n2;
        //n.push_back(last_node);

        //volume = this->get_volume(face.nodes);
        //n.clear();
        // if the volume formed by a face and the n2 point has negative value,
        // it means that n2 is outside of the element.
        if (volume < 0.0)
        {
          this->elems[e_id]->fringe = 1; // mark as fringed
          fringed_elm_id = e_id;
          break;
        }
      }*/
    }

    if (fringed_elm_id != -1)
    {
      break;
    }
  }
  return fringed_elm_id;
}

/*
Static_Mesh::to_graphviz
prints out a digraph based on graphviz pattern
it will become possible to visualize the kdtree
*/

void Static_Mesh::to_graphviz(void)
{
  std::ofstream output;

  output.open("./kdtree_digraph.dot");
  output << "digraph D {\n";
  output << "node [shape=box]\n";
  output << "graph [ordering=\"out\"]\n";

  for (auto &n : this->nodes)
  {
    output << n.id << " [label=\"" << n.id << " (" << n.coords[0] << "," << n.coords[1] << "," << n.coords[2] << ")\"]" << std::endl;
  }

  for (auto &n : this->nodes)
  {
    if (n.left > -1)
    {
      output << "     " << n.id << " -> L" << this->nodes[n.left].id << ";" << std::endl;
    }
    if (n.right > -1)
    {
      output << "     " << n.id << " -> R" << this->nodes[n.right].id << ";" << std::endl;
    }
  }
  output << "}";
  output.close();
}

/*
1
|-2 
| --3
|   |-6
|   --7
| 
--5
*/
void Static_Mesh::print_tree(long &root_, long &height, std::string before, int &flag_right)
{
  int i = 0;
  std::string space{" "};
  std::string down{"\u2502"};
  std::string left{"\u251c"};  // |-
  std::string right{"\u2514"}; // |_

  if (root_ != -1)
  {
    std::cout << root_ << " [";
    std::cout << this->nodes[root_].coords[0] << ", ";
    std::cout << this->nodes[root_].coords[1] << ", ";
    std::cout << this->nodes[root_].coords[2] << "]" << std::endl;

    if (this->nodes[root_].left != -1)
    {
      if (height > 0)
      {
        if (flag_right == 0)
          before += (down + space);
        if (flag_right == 1)
          before += (space + space);
      }
      if (this->nodes[root_].right != -1)
      {
        flag_right = 0;
        std::cout << before + left + " S:";
      }
      else
      {
        flag_right = 1;
        std::cout << before + right + " S:";
      }
      height++;
      this->print_tree(this->nodes[root_].left, height, before, flag_right);
    }

    if (this->nodes[root_].right != -1)
    {
      flag_right = 1;
      std::cout << before + right + " N:";
      height++;
      this->print_tree(this->nodes[root_].right, height, before, flag_right);
    }
  }
  before = before.substr(0, before.length() - 2);
  height--;
}

/*

build_kdtree(const std::vector<Node>& root, std::vector<std::shared_ptr<Node>>)

vector:	1(0,0) 7(1,0) 9(0.5,0.5) 4(0.2,0) 2(0.7,0.8) 5(0.9,1) 3(0.8,0.2) 8(1,1)

sample: 5(0.9,1) 3(0.8,0.2) 8(1,1)

order by dimLayer: 3(0.8,0.2) 5(0.9,1) 8(1,1)

median: 5(0.9,1)

save address into the root: root = &median

split:
1) <  0.9: 1 9 4 2 3
2) >= 0.9: 7 5 8

build_kdtree(root.left,  <1 9 4 2 3>, dimLayer);
build_kdtree(root.right, <7 5 8>, dimLayer);


*/

/*
How to select a subvector<Node> without copying data?

1) vector<shared_ptr<Node>> subvec that will hold the address of all nodes of interest for each run of my build_kdtree
- it might have less memory usage?
2) vector<int> that will hold all node indices to access the nodes of interest within my vector<Node>
- I will try this one first and see if it's ok


*/
