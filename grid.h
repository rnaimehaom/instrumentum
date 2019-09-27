#include "node.h"
#include "molecule.h"

#ifndef _gridh
#define _gridh

typedef struct {
  int atom;
  int node;
  int locale;
} State;

/// A class representing the entire grid of nodes on which molecules will be built. 
class Grid {
 private:
  /// This integer property represents the number of nodes extending 
  /// in the positive x direction. 
  int D1 = 0;
  /// This integer property represents the number of nodes extending 
  /// in the positive y direction. 
  int D2 = 0;
  /// This integer property represents the number of nodes extending 
  /// in the positive z direction. 
  int D3 = 0;
  /// The total number of nodes in the grid and thus 
  /// the length of the Grid::nodes array - it should 
  /// be equal to (2*D1 + 1)*(2*D2 + 1)*(2*D3 + 1).
  int total = 0;
  /// This integer property stores the size of the molecular skeleton, 
  /// based on the location of the pharmacophoric nodes that must be included, 
  /// in the positive x direction, so rs1 <= D1 and normally rs1 will 
  /// be significantly less than D1. It allows a variety of Grid methods 
  /// like blank_pharmacophore and fill_interior) to run more quickly 
  /// by reducing the trip count of loops. 
  int rs1 = 0;
  /// This integer property stores the size of the molecular skeleton,
  /// based on the location of the pharmacophoric nodes that must be included,  
  /// in the positive y direction, so rs2 <= D2 and normally rs2 will 
  /// be significantly less than D2. It allows a variety of Grid methods 
  /// like blank_pharmacophore and fill_interior) to run more quickly 
  /// by reducing the trip count of loops. 
  int rs2 = 0;
  /// This integer property stores the size of the molecular skeleton,
  /// based on the location of the pharmacophoric nodes that must be included,  
  /// in the positive z direction, so rs3 <= D3 and normally rs3 will 
  /// be significantly less than D3. It allows a variety of Grid methods 
  /// like blank_pharmacophore and fill_interior) to run more quickly 
  /// by reducing the trip count of loops. 
  int rs3 = 0;
  /// The geometric distance separating the atoms in the 
  /// tetrahedral mesh, measured in Angstroms.  
  double bond_length = 0.0;
  /// This vector contains a list of the atoms that are members of rings and 
  /// is calculated by the Grid::ring_analysis method.
  std::vector<int> ring_info;
  /// This vector stores which nodes are pharmacophoric, i.e. they must be 
  /// included in any molecular skeleton built by this software.
  std::vector<int> pharma_nodes;
  std::vector<State> backup[4];
  /// The main property of this class, a pointer of type Node which 
  /// stores the ensemble of nodes associated with the tetrahedral mesh. 
  Node* nodes;

  /// This method accepts the three mesh coordinates of a node and returns the index to this node in the array Grid::nodes. 
  inline int index1(int,int,int) const;
  /// This method accepts as its arguments the indices of two nodes and calculates the square of the geometric distance separating them. 
  inline double distance(int,int) const;
  /// This method sets the Node::locale and Node::atomic_number properties to zero for each node in the grid. 
  inline void clear();
  /// This method calls the Grid::ring_analysis method and returns its output value.
  inline int ring_count() {return ring_analysis();};
  /// This method accepts as its first three arguments the dimensional indices of a node in the mesh, while the final integer argument represents the state which lies between 1 and 4. The final argument contains the coordinates of this node in the grid. 
  void next_door(int,int,int,int,const double*);
  /// This method is called by the constructor and builds the tetrahedral mesh with the node coordinates and the neighbour relations using the Grid::next_door method. 
  void initialize();
  /// This method is called by the constructor and takes as its unique argument the number of pharmacophoric nodes and then allocates the memory for the Grid::pharma_nodes and the Grid::nodes properties.
  void allocate(int);
  /// This method converts any neighbouring nodes of carbon atoms which are empty to be hydrogen atoms. 
  void add_hydrogens();
  bool initial_deletion(double,int);
  bool path_selection(bool);
  bool secondary_deletion(int,int,int,int);
  bool rationalize(double,int,int);
  /// This method computes the number of rings in the current molecular scaffold and returns this value, as well as writing all of the nodes which are ring members to the Grid::ring_info vector. 
  int ring_analysis();
  void fill_interior();
  void restore(int);
  void save_state(int);
  bool connect_pharmacophores();
  void blank_pharmacophore(double);
  bool create_scaffold();
  /// This method writes the hydrocarbon scaffold in the current grid to an instance of the Molecule class.
  void write_scaffold(Molecule*) const;
 public:
  /// The standard constructor for this class - the first three arguments are the grid dimensions D1, D2 and D3, while the remaining two are the bond_length and number of pharmacophore nodes.
  Grid(int,int,int,double = 1.4,int = 0);
  /// The destructor for this class which frees the memory associated with the Grid::nodes property.
  ~Grid();
  friend class Molecular_Assembler;
};

int Grid::index1(int i,int j,int k) const
{
  int output = (k+D3) + (1+2*D3)*(j+D2) + (1+2*D3)*(1+2*D2)*(i+D1);
  return output;
}

void Grid::clear()
{
  for(int i=0; i<total; ++i) {
    nodes[i].locale = 0;
    nodes[i].atomic_number = 0;
  }
}

double Grid::distance(int in1,int in2) const
{
  double x[3],y[3],output;
  x[0] = nodes[in1].x;
  x[1] = nodes[in1].y;
  x[2] = nodes[in1].z;
  y[0] = nodes[in2].x;
  y[1] = nodes[in2].y;
  y[2] = nodes[in2].z;
  output = (x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]) + (x[2] - y[2])*(x[2] - y[2]);
  return output;
}
#endif
