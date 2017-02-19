#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <sqlite3.h>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef __globalh
#define __globalh

void ename(int,char*);
unsigned int irandom(unsigned int);
double rrandom();
void initialize_generator(unsigned int);
bool g_connected(const std::vector<int>&);
int get_index(int,const std::vector<int>&);
void shuffle(std::vector<int>&);
bool parallel(const double*,const double*);
void rperception(const std::vector<int>&,std::vector<int>&);
#endif

