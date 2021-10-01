#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <set>
#include <map>
#include <cassert>
#include <chrono>
#include <vector>
#include <string>
#include <random>
#include <mutex>
#include <thread>
#include <algorithm>
#include <unistd.h>
// To handle the database interactions
#include <sqlite3.h>
// To read the XML parameter file
#include <pugixml.hpp>

#ifndef _globalh
#define _globalh
bool file_exists(const std::string&);
bool connected(const std::vector<int>&);
void ring_perception(const std::vector<int>&,std::vector<int>&);

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

