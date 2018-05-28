#include "global.h"

#ifndef _nodeh
#define _nodeh

class Node {
 public:
  std::vector<int> neighbours;
  int state = 0;
  int atomic_number = 0;
  int atom_index = -1;
  int locale = 0;
  int visited = -1;
  int path_hop = -1;
  double x = std::numeric_limits<double>::quiet_NaN();
  double y = std::numeric_limits<double>::quiet_NaN(); 
  double z = std::numeric_limits<double>::quiet_NaN();

  Node();
  Node(const Node&);
  ~Node();
  Node& operator =(const Node&);
};
#endif

