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
  assert(seed > 0);
  gen->seed(seed);
}

int Random::irandom(int nmax)
{
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
#ifdef DEBUG
  assert(v.size() > 1);
#endif
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
