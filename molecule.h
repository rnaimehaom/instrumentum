#include "random.h"

#ifndef _moleculeh
#define _moleculeh

class Molecule {
 private:
  int natoms = 0;
  int nrings = 0;
  int p_allocated = 0;
  char ops[7] = {'P','N','S','O','T','F','A'};
  std::string opstring;
  std::vector<int> atom_type;
  std::vector<int> locale;
  std::vector<int> bonds;
  std::vector<int> btype;
  std::vector<int> rbonds;
  std::vector<int> rings;
  std::vector<double> coords;
  std::vector<int>* pieces;

  bool consistent() const;
  bool valence_check() const;
  bool add_oxygen();
  bool add_sulfur();
  bool add_nitrogen();
  bool create_dbond();
  bool create_tbond();
  bool create_penta1();
  bool create_amide();
  void aromatize(std::vector<int>&);
  bool fungrp();
  void get_rings();
  void axial_ring_bonds(std::vector<int>&) const;
  int eliminate_atoms(int*,int);
  int eliminate_atoms(int*,int,std::vector<int>&);
  inline int get_bindex(int,int) const;
  inline int get_rindex(int,int) const;
  inline bool in_ring(int) const;
  inline bool in_aromatic(int) const;
  bool is_aromatic(int) const;
  bool normalize_aromatic_bonds();
  void saturation_check() const;
  void normalize_free_ring(int);
  bool normalize_safe(const std::vector<int>&,bool*);
  void connected_components(int);
  void propagate(std::vector<int>&,int) const;
  void add_atom(int);
  void add_atom(int,const double*,int);
  void add_bond(int,int,int);
  void drop_atom(int);
 public:
  Molecule();
  Molecule(const Molecule&);
  ~Molecule();
  Molecule& operator =(const Molecule&);
  void clear();
  bool decorate(const bool*);
  std::string to_MDLMol() const;
  inline std::string get_opstring() const {return opstring;}; 
  friend std::ostream& operator <<(std::ostream&,const Molecule&);
  friend class Grid;
};

int Molecule::get_rindex(int n,int s) const
{
  for(int i=0; i<6; ++i) {
    if (rings[6*s+i] == n) return i;
  }
  return -1;
}

int Molecule::get_bindex(int n,int s) const
{
  for(int i=0; i<4; ++i) {
    if (bonds[4*s+i] == n) return i;
  }
  return -1;
}

bool Molecule::in_ring(int x) const
{
  // Return one if this atom is contained inside at least one ring
  int i,j;
  for(i=0; i<nrings; ++i) {
    for(j=0; j<6; ++j) {
      if (x == rings[6*i+j]) return true;
    }
  }
  return false;
}

bool Molecule::in_aromatic(int x) const
{
  // Return one if this atom is contained inside at least one aromatic ring
  for(int i=0; i<4; ++i) {
    if (btype[4*x+i] == 4) return true;
  }
  return false;
}
#endif
