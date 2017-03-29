#include "node.h"

Node::Node()
{
  neighbours.reserve(12);
  x = std::numeric_limits<double>::quiet_NaN();
  y = std::numeric_limits<double>::quiet_NaN();
  z = std::numeric_limits<double>::quiet_NaN();
  state = 0;
  locale = 0;
  visited = -1;
  path_hop = -1;
  atomic_number = 0;
  atom_index = -1;
}

Node::Node(const Node& source)
{
  neighbours = source.neighbours;
  x = source.x;
  y = source.y;
  z = source.z;
  state = source.state;
  locale = source.locale;
  visited = source.visited;
  path_hop = source.path_hop;
  atomic_number = source.atomic_number;
  atom_index = source.atom_index;
}

Node::~Node()
{

}

Node& Node::operator = (const Node& source)
{
  if (this == &source) return *this;
  neighbours = source.neighbours;
  x = source.x;
  y = source.y;
  z = source.z;
  state = source.state;
  locale = source.locale;
  visited = source.visited;
  path_hop = source.path_hop;
  atomic_number = source.atomic_number;
  atom_index = source.atom_index;
  return *this;
}
