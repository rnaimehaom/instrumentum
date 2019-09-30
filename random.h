#include "instrumentum.h"

#ifndef __randomh
#define __randomh

/// A utility class representing a pseudo-random number generator and related methods.
class Random {
 private:
  /// This property is an instance of the 2011 C++ standard 
  /// for generating pseudo-random numbers. 
  std::random_device rd;
  /// This property is the pseudo-random number generator, the well-known 
  /// Mersenne twister with a period length of 19937, part of 2011 C++ 
  /// standard. 
  std::mt19937* gen;
  /// This property is responsible for generating uniform real variates, 
  /// which are the foundation of this class' generation of a wider set 
  /// of random variates. 
  std::uniform_real_distribution<>* VRG;

 public:
  /// The standard constructor for this class, allocating the memory for the Random::gen and Random::VRG properties and initializing the random number seed to the current time.
  Random();
  /// The destructor for this class, freeing the memory associated with the Random::gen and Random::VRG properties.
  ~Random();
  /// This method sets the seed for the pseudo-random number process to the argument, which must be positive. 
  void initialize_generator(long);
  /// This method returns a uniform random variate on the interval [0,1).  
  double drandom();
  /// This method returns a uniform random variate on the interval [0,n) where n is the method's argument, which must be positive.
  int irandom(int);
  /// This method takes its vector argument, which must have a length greater than one, and shuffles the elements using the Fisher-Yates algorithm. 
  void shuffle(std::vector<int>&);
};
#endif
