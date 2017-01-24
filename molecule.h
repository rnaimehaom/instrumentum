#include "global.h"

#ifndef __moleculeh
#define __moleculeh

class Molecule {
 private:
  std::vector<int> atoms;
  std::vector<int> atom_type;
  std::vector<int> locale;
  std::vector<int> bonds;
  std::vector<int> btype;
  std::vector<int> rbonds;
  std::vector<int>* pieces;
  std::vector<int> rings;
  std::vector<double> coords;
  std::string opstring;
  unsigned int nrings;
  unsigned int p_allocated;
  char ops[7];

  int atom_index(int) const;
  bool add_oxygen();
  bool add_sulfur();
  bool add_nitrogen();
  bool create_dbond();
  bool create_tbond();
  bool create_penta1();
  bool create_amide();
  void aromatize(std::vector<int>&);
  bool fungrp();
  void valence_check() const;
  void get_rings();
  void axial_ring_bonds(std::vector<int>&) const;
  int get_bindex(int,int) const;
  int get_rindex(int,int) const;
  int eliminate_atoms(int*,int);
  int eliminate_atoms(int*,int,std::vector<int>&);
  bool in_ring(int) const;
  bool in_aromatic(int) const;
  bool is_aromatic(int) const;
  bool normalize_aromatic_bonds();
  void resequence(int);
  void saturation_check() const;
  void normalize_free_ring(int);
  bool normalize_safe(const std::vector<int>&,bool*);
  void connected_components(unsigned int);
  void propagate(std::vector<unsigned int>&,unsigned int) const;
  double minimize(unsigned int,double);
  std::string write2string() const;

 public:
  Molecule();
  Molecule(const Molecule&);
  ~Molecule();
  Molecule& operator =(const Molecule&);
  void write2sdf(const char*) const;
  void write2db(unsigned int,pqxx::connection&);
  void add_atom(int,int);
  void add_atom(int,int,const double*,int);
  void add_bond(int,int,int);
  void dump_molecule() const;
  void clear();
  bool decorate(const bool*); 
};
#endif
