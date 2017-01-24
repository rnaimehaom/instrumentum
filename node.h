#include "global.h"

#ifndef __node_h
#define __node_h

class Node {
 public:
  std::vector<int> neighbours;
  int state;
  double x,y,z;
  int atomic_number;
  int locale;
  int visited;
  int path_hop;

  Node();
  Node(const Node&);
  ~Node();
  Node& operator =(const Node&);
};
#endif

