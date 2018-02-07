#include "node.h"
#include "molecule.h"

#ifndef _gridh
#define _gridh

typedef struct {
  int atom;
  int node;
  int locale;
} State;

class Grid {
 private:
  int D1,D2,D3;
  int rs1,rs2,rs3;
  Node* nodes;
  std::vector<int> ring_info;
  unsigned int total;
  std::vector<State> backup[4];
  std::vector<int> pnodes;
  double bond_length;

  inline int index1(int,int,int) const;
  inline double distance(int,int) const;
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
  bool connect_pharmacophores();
  void blank_pharmacophore(double);
  bool create_scaffold();
  void write_scaffold(Molecule*) const;
 public:
  Grid(unsigned int);
  Grid(double,unsigned int);
  Grid(int,int,int,unsigned int);
  Grid(int,int,int,double,unsigned int);
  ~Grid();
  friend class Molecular_Assembler;
};

int Grid::index1(int i,int j,int k) const
{
  int output = (k+D3) + (1+2*D3)*(j+D2) + (1+2*D3)*(1+2*D2)*(i+D1);
  return output;
}

double Grid::distance(int in1,int in2) const
{
  double x[3],y[3],output;
  x[0] = nodes[in1].x;
  x[1] = nodes[in1].y;
  x[2] = nodes[in1].z;
  y[0] = nodes[in2].x;
  y[1] = nodes[in2].y;
  y[2] = nodes[in2].z;
  output = (x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]);
  return output;
}
#endif
