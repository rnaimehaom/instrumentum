#include "node.h"

Node::Node()
{
  neighbours.reserve(12);
}

Node::Node(const Node& source)
{
  neighbours = source.neighbours;
  x = source.x;
  y = source.y;
  z = source.z;
  state = source.state;
  locale = source.locale;
  atomic_number = source.atomic_number;
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
  atomic_number = source.atomic_number;

  return *this;
}
