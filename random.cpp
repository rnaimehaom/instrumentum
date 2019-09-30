#include "random.h"

Random::Random()
{
  gen = new std::mt19937(rd());
  VRG = new std::uniform_real_distribution<>(0.0,1.0);
  gen->seed(long(std::time(NULL)));
}

Random::~Random()
{
  delete gen;
  delete VRG;
}

void Random::initialize_generator(long seed)
{
  if (seed <= 0) throw std::invalid_argument("The random number seed must be positive!");
  gen->seed(seed);
}

int Random::irandom(int nmax)
{
  if (nmax <= 0) throw std::invalid_argument("The argument to Random::irandom must be positive!");
  int output = int(double(nmax)*drandom());
  return output;
}

double Random::drandom()
{
  double out = (*VRG)(*gen);
  return out;
}

void Random::shuffle(std::vector<int>& v)
{
  if (v.size() < 2) throw std::invalid_argument("The vector to be shuffled must have at least two elements!");
  // Fisher-Yates shuffle algorithm
  int i,j,q;
  const int n = (signed) v.size();

  for(i=0; i<n; ++i) {
    v[i] = i;
  }
  for(i=n-1; i>0; --i) {
    j = irandom(1+i);
    q = v[j];
    v[j] = v[i];
    v[i] = q;
  }
}
