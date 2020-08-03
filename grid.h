#include "node.h"
#include "molecule.h"

#ifndef _gridh
#define _gridh

/// A utility class representing a summary of the state of the essential properties of a grid node, used for backing up the grid.
class Summary {
 public:
  /// This property stores a vector containing the index of the node in the Grid::nodes property for each node with a positive locale.
  std::vector<int> index;
  /// This property stores a vector containing the Node::atomic_number property of the grid nodes with a positive locale.
  std::vector<int> atom;
  /// This property stores a vector containing the Node::locale property of the grid nodes with a positive locale.
  std::vector<int> locale;
};

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
  /// This property stores the four different backup states of the grid, which 
  /// is used by the save_state() and restore_state() methods to make a snapshot 
  /// of the grid's state so that the Molecular_Assembler class can work on a 
  /// given hydrocarbon scaffold many times over.
  Summary backup[4];
  /// The main property of this class, a pointer of type Node which 
  /// stores the ensemble of nodes associated with the tetrahedral mesh. 
  Node* nodes;

  /// This method accepts the three mesh coordinates of a node and returns the index to this node in the array Grid::nodes. 
  int index1(int,int,int) const;
  /// This method accepts as its arguments the indices of two nodes and calculates the square of the geometric distance separating them. 
  double distance(int,int) const;
  /// This method sets the Node::locale and Node::atomic_number properties to zero for each node in the grid. 
  void clear();
  /// This method restores the state of the grid properties from one of the four save states in the Grid::backup property, specified by the method's argument.
  void restore_state(int);
  /// This method saves the current state of the grid properties (the Node::locale and Node::atomic_number properties of each node) to the Grid::backup array, with the element specified by the method's argument.
  void save_state(int);
  /// This method calls the Grid::ring_analysis method and returns its output value.
  int ring_count();
  /// This method accepts as its first three arguments the dimensional indices of a node in the mesh, while the final integer argument represents the state which lies between 1 and 4. The final argument contains the coordinates of this node in the grid. 
  void next_door(int,int,int,int,const double*);
  /// This method is called by the constructor and builds the tetrahedral mesh with the node coordinates and the neighbour relations using the Grid::next_door method. 
  void initialize();
  /// This method is called by the constructor and takes as its unique argument the number of pharmacophoric nodes and then allocates the memory for the Grid::pharma_nodes and the Grid::nodes properties.
  void allocate(int);
  /// This method converts any neighbouring nodes of carbon atoms which are empty to be hydrogen atoms. 
  void add_hydrogens();
  /// This method makes a series of deletions of carbon atoms in the initial volume of the molecule via randomly centred "explosions" that eliminate carbon atoms within a given radius, checking that the pharmacophoric connectedness is maintained after each such deletion. The deletions continue until a percentage (first argument) of the original number of carbon atoms is deleted or the maximum number of attempts (second argument) has been made. The method returns true in the first case and false in the second.
  bool initial_deletion(double,int);
  /// This method constructs a path in the grid linking together the pharmacophoric nodes, beginning with either such a node (when the argument is true) or with a random node (when it is false). The method returns true when it is successful in building such a path.
  bool path_selection(bool);
  /// This method carries out a series of secondary deletions of carbon atoms from the molecular volume - the maximum number of attempts is set by the final argument. The criteria for exiting the deletion loop are set by the method's first three arguments, specifying the maximum number of carbon atoms with four carbon neighbours, the maximum number of carbon atoms belonging to four separate rings and the maximum number of rings.
  bool secondary_deletion(int,int,int,int);
  /// This method rationalizes the scaffold by forcing any node that is the neighbour of two or more carbon atoms to become a carbon atom itself. It also eliminates a fixed percentage (the first argument) of randomly selected methyl groups as well as counting the number of rings - if this number is less than the second argument or greater than the third, the method returns false and true otherwise.
  bool rationalize(double,int,int);
  /// This method computes the number of rings in the current molecular scaffold and returns this value, as well as writing all of the nodes which are ring members to the Grid::ring_info vector. 
  int ring_analysis();
  /// This method fills in the interior of the the space (-rs1:rs1,-rs2:rs2,-rs3:rs3) with carbon atoms, if the Node::locale property is equal to two.
  void fill_interior();
  /// This method tests that all of the pharmacophoric nodes are connected together by the molecular scaffold, returning true in that case and false otherwise; carbon nodes which are disconnected from the scaffold are marked as silver for deletion.
  bool connect_pharmacophores();
  /// This method creates a "blank" pharmacophore whose spherical radius is the method's argument, selecting a set of random boundary nodes in number equal to the length of Grid::pharma_nodes.
  void blank_pharmacophore(double);
  /// This method writes the hydrocarbon scaffold in the current grid to an instance of the Molecule class.
  void write_scaffold(Molecule*) const;
 public:
  /// The standard constructor for this class - the first three arguments are the grid dimensions D1, D2 and D3, while the remaining two are the bond_length and number of pharmacophore nodes.
  Grid(int,int,int,double = 1.4,int = 0);
  /// The destructor for this class which frees the memory associated with the Grid::nodes property.
  ~Grid();
  friend class Molecular_Assembler;
};

inline int Grid::index1(int i,int j,int k) const
{
  int output = (k+D3) + (1+2*D3)*(j+D2) + (1+2*D3)*(1+2*D2)*(i+D1);
  return output;
}

inline void Grid::clear()
{
  for(int i=0; i<total; ++i) {
    nodes[i].locale = 0;
    nodes[i].atomic_number = 0;
  }
}

inline int Grid::ring_count()
{
  return ring_analysis();
}

inline void Grid::restore_state(int q)
{
  int i;
  unsigned int n,l;

  // First, initialize the whole grid back to zero...
  clear();

  // Now refill it from the vector created by save_state:
  n = backup[q].index.size();
  for(l=0; l<n; ++l) {
    i = backup[q].index[l];
    nodes[i].locale = backup[q].locale[l];
    nodes[i].atomic_number = backup[q].atom[l];
  }
}

inline void Grid::save_state(int q)
{
  backup[q].index.clear();
  backup[q].locale.clear();
  backup[q].atom.clear();

  for(int i=0; i<total; ++i) {
    if (nodes[i].locale <= 0) continue;
    backup[q].index.push_back(i);
    backup[q].locale.push_back(nodes[i].locale);
    backup[q].atom.push_back(nodes[i].atomic_number);
  }
}

inline double Grid::distance(int in1,int in2) const
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
