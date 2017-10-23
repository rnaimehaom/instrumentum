#include "global.h"

#ifndef _atomh
#define _atomh

class Atom {
 private:
  int atomic_number;
  char element[2];

  void set_element();

 public:
  Atom(int);
  Atom();
  int get_anumber() const;
  void set_anumber(int);
};
#endif

