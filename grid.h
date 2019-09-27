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
  std::vector<int> ring_info;
  /// This vector stores which nodes are pharmacophoric, i.e. they must be 
  /// included in any molecular skeleton built by this software.
  std::vector<int> pharma_nodes;
  std::vector<State> backup[4];
  /// The main property of this class, a pointer of type Node which 
  /// stores the ensemble of nodes associated with the tetrahedral mesh. 
  Node* nodes;

  inline int index1(int,int,int) const;
  inline double distance(int,int) const;
  void next_door(int,int,int,int,const double*);
  void initialize();
  void allocate(int);
  void clear();
  void add_hydrogens();
  bool initial_deletion(double,int);
  bool path_selection(bool);
  bool secondary_deletion(int,int,int,int);
  bool rationalize(double,int,int);
  int ring_count();
  int ring_analysis();
  void fill_interior();
  void restore(int);
  void save_state(int);
  bool connect_pharmacophores();
  void blank_pharmacophore(double);
  bool create_scaffold();
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

double Grid::distance(int in1,int in2) const
{
  double x[3],y[3],output;
  x[0] = nodes[in1].x;
  x[1] = nodes[in1].y;
  x[2] = nodes[in1].z;
  y[0] = nodes[in2].x;
  y[1] = nodes[in2].y;
  y[2] = nodes[in2].z;
  output = (x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]);
  return output;
}
#endif
