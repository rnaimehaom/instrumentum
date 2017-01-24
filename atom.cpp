#include "atom.h"

Atom::Atom()
{
  atomic_number = 0;
  element[0] = '\0';
}

int Atom::get_anumber() const
{
  return atomic_number;
}

void Atom::set_anumber(int n)
{
  atomic_number = n;
  set_element();
}

void Atom::set_element() {
  switch(atomic_number) {
  case 0:
    element[0] = '\0';
    break;
  case 1:
    element[0] = 'H';
    break;
  case 6:
    element[0] = 'C';
    break;
  case 7:
    element[0] = 'N';
    break;
  case 8:
    element[0] = 'O';
    break;
  case 9:
    element[0] = 'F';
    break;
  case 15:
    element[0] = 'P';
    break;
  case 16:
    element[0] = 'S';
    break;
  case 17:
    element[0] = 'C';
    element[1] = 'l';
    break;
  case 35:
    element[0] = 'B';
    element[1] = 'r';
    break;
  case 47:
    element[0] = 'A';
    element[1] = 'r';
    break;
  case 53:
    element[0] = 'I';
  }
}

Atom::Atom(int n)
{
  atomic_number = n;
  this->set_element();
}

