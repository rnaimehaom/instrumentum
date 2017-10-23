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

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _globalh
#define _globalh
const double epsilon = 0.001;

// Random number methods
double drandom();
void initialize_generator(unsigned long);
unsigned int irandom(unsigned int);

void element(int,char*);
bool connected(const std::vector<int>&);
bool file_exists(const std::string&);
void capitalize(std::string&);
int get_index(int,const std::vector<int>&);
void shuffle(std::vector<int>&);
bool parallel(const double*,const double*);
void ring_perception(const std::vector<int>&,std::vector<int>&);
#endif

