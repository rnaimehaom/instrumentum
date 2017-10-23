#include "global.h"

#ifndef _nodeh
#define _nodeh

class Node {
 public:
  std::vector<int> neighbours;
  int state;
  double x,y,z;
  int atomic_number;
  int atom_index;
  int locale;
  int visited;
  int path_hop;

  Node();
  Node(const Node&);
  ~Node();
  Node& operator =(const Node&);
};
#endif

