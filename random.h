#include "global.h"

#ifndef __randomh
#define __randomh

class Random {
 private:
  // Random number variables
  std::random_device rd;
  std::mt19937* gen;
  std::uniform_real_distribution<>* VRG;

 public:
  Random();
  ~Random();
  void initialize_generator(unsigned long);
  double drandom();
  unsigned int irandom(unsigned int);
  void shuffle(std::vector<int>&);
};
#endif
