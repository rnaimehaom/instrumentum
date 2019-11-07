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
  /// This property is a vector of the element of each atom in the molecule, 
  /// by atomic number.
  std::vector<int> atom_type;
  /// This property is a vector of the mutability of each atom in the molecule, 
  /// distinguishing pharmacophoric atoms (which cannot be changed) for example. 
  std::vector<int> locale;
  /// This property is a vector storing the bond table of the molecule's atoms, 
  /// indicating which atom is bonded to which others, by the internal index of 
  /// the atom within the molecule. 
  std::vector<int> bonds;
  /// This property is a vector storing the type of the molecule's bonds using the 
  /// notation one for a single shared electron, two for a double bond, three for a 
  /// triple bond and four for an aromatic ring bond.  
  std::vector<int> btype;
  /// This property is a vector storing the molecule bonds which are part of a 
  /// ring. 
  std::vector<int> rbonds;
  /// This property is a vector storing the individuals atoms associated with each 
  /// of the molecule's hexa- and penta-atomic rings. 
  std::vector<int> rings;
  /// This property is a vector storing the geometric coordinates of the molecule's 
  /// atoms, in Angstroms; its length is thus three times the value of Molecule::natoms. 
  std::vector<double> coords;
  /// This array of integer vectors stores the atom indices by connected ring component 
  /// of the molecule.
  std::vector<int>* pieces;
  /// This character array lists the seven different "decorating" operations that 
  /// can be performed on a molecule: O (oxygen substitution), N (nitrogen substitution), 
  /// S (sulfur substitution), T (creation of a triple bond), F (halogen substitution), 
  /// A (amide substitution) and P (conversion of a hexa-atomic to a penta-atomic ring). 
  /// These are the letters of the alphabet which is used to write the value of the 
  /// Molecule::opstring parameter.    
  static const char ops[7];

  /// Given a pair of atoms in the molecule, specified by the method's two arguments, this method returns the index of this bond in the property Molecule::bonds and -1 if no such bond can be found. 
  inline int get_bindex(int,int) const;
  /// Given a pair of atoms in the molecule, specified by the method's two arguments, this method returns the index of the ring in Molecule::rings that contains the two atoms and -1 if no such ring can be found. 
  inline int get_rindex(int,int) const;
  /// This method returns true if the atom whose index is the method's unique argument is a member of a ring, false otherwise. 
  inline bool in_ring(int) const;
  /// This method returns true if the atom whose index is the method's unique argument is a member of an aromatic ring, false otherwise. 
  inline bool in_aromatic(int) const;
  /// This method performs a basic check of the bond valences for the molecule's carbon and hydrogen atoms, which are the main kind of atoms we expect to find in the molecule when this method is called at the beginning of the Molecule::decorate() method; it returns true if all is well and false otherwise.
  bool saturation_check() const;
  /// This method performs a broad check on the consistency of the various Molecule properties in terms of the lengths of the vector properties and the number of atoms, as well as the size of the geometric coordinates and the bond table values. The method returns true if all is well and false otherwise.
  bool consistent() const;
  /// This method performs a check on the valence of the molecule's hydrogen, phosphorus, chlorine and fluorine atoms, returning true if all is well and false otherwise.
  bool valence_check() const;
  /// This method accepts the index of one of the molecule's rings and returns true if this ring is aromatic, false otherwise.
  bool is_aromatic(int) const;
  /// This method paints all of the atoms in the molecule a given colour (the second argument), if their value in the first argument is zero and they belong to a ring which contains an atom of this colour. 
  void propagate(std::vector<int>&,int) const;
  /// This method computes the all axial (as opposed to equatorial) neighbours of the molecule's ring atoms and stores their index in the method's argument. 
  void axial_ring_bonds(std::vector<int>&) const;
  /// This method attempts to convert one of the molecule's atoms to an oxygen atom and returns true if successful, false otherwise. 
  bool add_oxygen();
  /// This method attempts to convert one of the molecule's atoms to a sulfur atom and returns true if successful, false otherwise. 
  bool add_sulfur();
  /// This method attempts to convert one of the molecule's atoms to a nitrogen atom and returns true if successful, false otherwise. 
  bool add_nitrogen();
  /// This method attempts to convert one of the molecule's single bonds to a double bond and returns true if successful, false otherwise. 
  bool create_dbond();
  /// This method attempts to convert one of the molecule's double bonds to a triple bond and returns true if successful, false otherwise. 
  bool create_tbond();
  /// This method attempts to convert one of the molecule's hexa-atom rings to a penta-atom ring and returns true if successful, false otherwise. 
  bool create_penta1();
  /// This method attempts to convert one of the molecule's atoms to an amide group and returns true if successful, false otherwise. 
  bool create_amide();
  /// This method accepts as its argument the vector of axial neighbours of the molecule's ring atoms (computed using Molecule::axial_ring_bonds) and then attempts to aromatize one of the molecule's rings. It returns the number of rings that were aromatized successfully.
  int aromatize(std::vector<int>&);
  /// This method attempts to convert one of the molecule's atoms to a halogen atom and returns true if successful, false otherwise. 
  bool fungrp();
  /// This method identifies all of the rings in the molecule and returns the number of such rings; it also populates the Molecule::rings property.
  int get_rings();
  /// This method method removes from the molecule the atoms stored by index in the method's unique argument. The method iterates in reverse order over the set of atoms and then the Molecule::drop_atom method is called on each member. The method returns true if removing the atoms succeeds and false otherwise. 
  bool eliminate_atoms(std::set<int>&);
  /// This method removes from the molecule the atoms stored by index in the first argument. The method iterates in reverse order over the set of atoms and then the Molecule::drop_atom method is called on each member. The third argument is a vector of axial atoms whose indexing is updated by this method if necessary; the method returns true if removing the atoms succeeds and false otherwise. 
  bool eliminate_atoms(std::set<int>&,std::vector<int>&);
  /// This method converts the molecule's aromatic ring bonds, of type equal to four, and changes them to alternative values of one and two, necessary for being output to an MDL MOL file; it returns true if the bond type conversion is successful.
  bool normalize_aromatic_bonds();
  /// This method converts the aromatic ring bonds of an isolated ring whose index is specified by the method's argument, changing the bond type from four to an alternating pattern of one and two values, needed to satisfy the requirements of the MDL MOL file format.  
  void normalize_free_ring(int);
  /// This method uses a very careful and self-consistent method to convert aromatic ring bonds from being of type four to type one and two in an alternating pattern and returns true if it is successful, false otherwise. The first argument is a vector containing the aromatic rings to convert by ring index and the second argument is a Boolean array of the same length as the first argument. As aromatic rings are normalized, the corresponding array element is set from false to true. 
  bool normalize_safe(const std::vector<int>&,bool*);
  /// This method computes the connected ring components of the molecule, using the propagate() method, where the argument is the number of vertices (i.e. atoms), with the output written to the Molecule::pieces property.
  void connected_components(int);
  /// This method removes all reference to the atom whose index is the method's argument, reindexing internal arrays to correspond to the fact that the Molecule::natoms property has been decremented. It returns false if the atom doesn't exist and true otherwise.
  bool drop_atom(int);
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
  int write(std::ofstream&) const;
  int read(std::ifstream&);
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
