#include "global.h"
#include "node.h"
#include "molecule.h"

#ifndef __gridh
#define __gridh

class state {
 public:
  int atom;
  int node;
  int locale;
};

class Grid {
 private:
  int D1,D2,D3;
  int rs1,rs2,rs3;
  Node* nodes;
  std::vector<int> ring_info;
  unsigned int total;
  std::vector<state> backup[4];
  std::vector<int> pnodes;
  double bond_length;

  int index1(int,int,int) const;
  void next_door(int,int,int,int,const double*);
  void initialize();
  void allocate(unsigned int);
  void set_default_values();
  void clear();
  void add_hydrogens();
  bool initial_deletion(double,unsigned int);
  bool path_selection(bool);
  bool secondary_deletion(unsigned int,unsigned int,unsigned int,unsigned int);
  bool rationalize(double,unsigned int,unsigned int);
  unsigned int ring_count();
  unsigned int ring_analysis();
  void fill_interior();
  void restore(int);
  void save_state(int);
  bool connected();
  double distance(int,int) const;

 public:
  Grid(unsigned int);
  Grid(double,unsigned int);
  Grid(int,int,int,unsigned int);
  Grid(int,int,int,double,unsigned int);
  ~Grid();
  void process_pharmacophore(const char*);
  void blank_pharmacophore(double);
  bool create_scaffold();
  void write_scaffold(Molecule*) const;
  friend class Molecular_Assembler;
};
#endif
