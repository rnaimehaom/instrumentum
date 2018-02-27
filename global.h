#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <set>
#include <cassert>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
// To handle the database interactions
#include <sqlite3.h>
#include <boost/filesystem.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _globalh
#define _globalh
void element(int,char*);
bool connected(const std::vector<int>&);
void ring_perception(const std::vector<int>&,std::vector<int>&);

inline void capitalize(std::string& s)
{
  std::transform(s.begin(),s.end(),s.begin(),::toupper);
}

inline int get_index(int x,const std::vector<int>& v)
{
  std::vector<int>::const_iterator it = std::find(v.begin(),v.end(),x);
  if (it == v.end()) return -1;
  int output = it - v.begin();
  return output;
}

inline bool parallel(const double* x,const double* y)
{
  // We take the cross product of these two 3-vectors and see if
  // it's approximately zero
  bool output = false;
  double delta,out[3];
  
  out[0] = x[1]*y[2] - x[2]*y[1];
  out[1] = x[2]*y[0] - x[0]*y[2];
  out[2] = x[0]*y[1] - x[1]*y[0];
  delta = out[0]*out[0] + out[1]*out[1] + out[2]*out[2];

  if (delta < 1e-06) output = true;

  return output;
}
#endif

