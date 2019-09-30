#include "random.h"

#ifndef _moleculeh
#define _moleculeh

/// A class representing a small (approximately 25 atoms) organic molecule, created by "decorating" a hydrocarbon skeleton.  
class Molecule {
 private:
  /// This integer property is the number of atoms in the 
  /// molecule. 
  int natoms = 0;
  /// This integer property is the number of distinct rings in the 
  /// molecule. 
  int nrings = 0;
  /// This integer property is the number of connected components of 
  /// the molecule and is also the length of the Molecule::pieces array. 
  int p_allocated = 0;
  /// This string property is a list of all of the decorative operations
  /// performed successfully on the molecule. 
  std::string opstring = "";
  std::vector<int> atom_type;
  std::vector<int> locale;
  std::vector<int> bonds;
  std::vector<int> btype;
  std::vector<int> rbonds;
  std::vector<int> rings;
  std::vector<double> coords;
  std::vector<int>* pieces;
  static const char ops[7];

  inline int get_bindex(int,int) const;
  inline int get_rindex(int,int) const;
  inline bool in_ring(int) const;
  inline bool in_aromatic(int) const;
  void saturation_check() const;
  bool consistent() const;
  bool valence_check() const;
  bool is_aromatic(int) const;
  void propagate(std::vector<int>&,int) const;
  void axial_ring_bonds(std::vector<int>&) const;
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
  bool eliminate_atoms(int*,int);
  bool eliminate_atoms(int*,int,std::vector<int>&);
  bool normalize_aromatic_bonds();
  void normalize_free_ring(int);
  bool normalize_safe(const std::vector<int>&,bool*);
  void connected_components(int);
  /// This method removes all reference to the atom whose index is the method's argument, reindexing internal arrays to correspond to the fact that the Molecule::natoms property has been decremented. 
  void drop_atom(int);
 public:
  /// The default constructor for this class, which does nothing.
  Molecule();
  /// The copy constructor, which copies over all of the properties of the source instance, deleting and re-allocating the Molecule::pieces property if necessary.
  Molecule(const Molecule&);
  /// The destructor for this class - if Molecule::p_allocted is greater than zero then the memory associated with the Molecule::pieces array is returned. 
  ~Molecule();
  /// The overloaded assignment operator, which copies over all of the properties of the source instance, deleting and re-allocating the Molecule::pieces property if necessary.
  Molecule& operator =(const Molecule&);
  /// This method restores all of the properties of the class to their default values.
  void clear();
  /// This method is the main public one - it takes the molecule, which begins as one consisting exclusively of carbon and hydrogen atoms linked by single bonds, and adds double and triple bonds, aromatic rings, oxygen, nitrogen, sulfur and halogen atoms and so forth. The argument controls which of these decoration operations will be attempted. 
  bool decorate(const bool*);
  /// This method adds an atom to the molecule - the first argument is the atom's type (by atomic number), the second is the atom's geometric coordinates and the final argument is the atom's locale (e.g. whether it is a pharmacophoric atom).
  void add_atom(int,const double*,int);
  /// This method adds a bond, the type of which is the method's final argument, between the two atoms specified by the method's first two arguments.
  void add_bond(int,int,int);
  /// This method writes the content of the molecule to a string corresponding to the MDL MOL format for the molecule, which is how the molecule is stored in an SQLite database. 
  std::string to_MDLMol() const;
  /// This method returns the current value of the Molecule::opstring property.
  inline std::string get_opstring() const {return opstring;}; 
  /// This overloaded ostream operator writes the essential properties (atom type, coordinates, bond table and ring structure) of the instance of this class to the screen.
  friend std::ostream& operator <<(std::ostream&,const Molecule&);
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
