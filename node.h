#include "instrumentum.h"

#ifndef _nodeh
#define _nodeh

/// A class representing a single node (i.e. potential atom location) in the tetrahedral mesh. 
class Node {
 private:
  /// This STL vector property stores the three integer co-ordinates of the four neighbours 
  /// of the node and so normally has a length of 12 unless the node is on the boundary of the 
  /// mesh.
  std::vector<int> neighbours;
  /// This integer property is used to store the dimensional orientation of this node during 
  /// the construction of the tetrahedral mesh by the Grid::next_door method.
  int state = 0;
  /// This integer property stores the atomic number of the atom which is located at this node,
  /// based on the value from the periodic table. A value of zero indicates that the node is empty.
  int atomic_number = 0;
  /// This integer property is used to distinguish boundary nodes and pharmacophoric nodes from 
  /// "regular" nodes when building the molecular skeleton. 
  int locale = 0;
  /// This floating point property stores the first of the node's spatial 
  /// coordinates.
  double x = std::numeric_limits<double>::quiet_NaN();
  /// This floating point property stores the second of the node's spatial 
  /// coordinates.
  double y = std::numeric_limits<double>::quiet_NaN(); 
  /// This floating point property stores the third of the node's spatial 
  /// coordinates.
  double z = std::numeric_limits<double>::quiet_NaN();

 public:
  /// The standard default constructor for this class which reserves a length of 12 for the Node::neighbours property.
  Node();
  /// The standard copy constructor which simply copies over the value of the properties of the source instance of this class.
  Node(const Node&);
  /// The standard destructor for the Node class which in this case does nothing.
  ~Node();
  /// The standard assignment operator which simply assigns the value of the properties of the source instance of this class to the target instance.
  Node& operator =(const Node&);
  friend class Grid;
};
#endif

