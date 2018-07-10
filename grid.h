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
  int D1 = 17;
  int D2 = 17;
  int D3 = 9;
  int total = 23275;
  int rs1 = 0;
  int rs2 = 0;
  int rs3 = 0;
  double bond_length = 1.4;
  std::vector<int> ring_info;
  std::vector<int> pharma_nodes;
  std::vector<State> backup[4];
  Node* nodes;

  inline int index1(int,int,int) const;
  inline double distance(int,int) const;
  void next_door(int,int,int,int,const double*);
  void initialize();
  void allocate(int);
  void clear();
  void add_hydrogens();
  bool initial_deletion(double,int);
  bool path_selection(bool);
  bool secondary_deletion(int,int,int,int);
  bool rationalize(double,int,int);
  int ring_count();
  int ring_analysis();
  void fill_interior();
  void restore(int);
  void save_state(int);
  bool connect_pharmacophores();
  void blank_pharmacophore(double);
  bool create_scaffold();
  void write_scaffold(Molecule*) const;
 public:
  Grid(int);
  Grid(double,int);
  Grid(int,int,int,int);
  Grid(int,int,int,double,int);
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
