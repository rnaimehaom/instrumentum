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
}

Node::Node(const Node& q)
{
  neighbours = q.neighbours;
  x = q.x;
  y = q.y;
  z = q.z;
  state = q.state;
  locale = q.locale;
  visited = q.visited;
  path_hop = q.path_hop;
  atomic_number = q.atomic_number;
}

Node::~Node()
{

}

Node& Node::operator = (const Node& q)
{
  if (this == &q) return *this;
  neighbours = q.neighbours;
  x = q.x;
  y = q.y;
  z = q.z;
  state = q.state;
  locale = q.locale;
  visited = q.visited;
  path_hop = q.path_hop;
  atomic_number = q.atomic_number;
  return *this;
}
