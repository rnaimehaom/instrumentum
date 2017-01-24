#include "global.h"

#ifndef __tokenh
#define __tokenh

class Token {
 private:
  int hop_count;
  std::vector<int> visited;

  void one_hop();
  bool exists() const;
 public:
  static int rcount;
  static int nedges;
  static int** edges;
  static int** rinfo;
  
  Token();
  Token(int); 
  Token(const std::vector<int>&); 
  ~Token();
  void circuit();
  void clear();
  void initialize(int);
};
#endif

