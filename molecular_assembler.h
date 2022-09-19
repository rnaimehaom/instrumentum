#include "grid.h"

#ifndef _masmh
#define _masmh

/// A class representing the director of this program, reading in the parameter file and assembling the molecules that are written to a database.
class Molecular_Assembler {
 private:
  /// This positive integer property is the maximum number of initial deletion attempts;
  /// it is the second argument to the Grid::initial_deletion() method. The value of this
  /// property is read from the XML parameter <CODE>MaximumAttempts</CODE>.
  int max_attempts = 100;
  /// This non-negative integer property is the number of carbon atoms with four carbon
  /// neighbours that will be tolerated in a candidate molecule; it is the first argument 
  /// to the Grid::secondary_deletion() method. The value of this property is read from
  /// the XML parameter <CODE>QuaternaryCarbonAtoms</CODE>. 
  int nc4 = 2;
  /// This non-negative integer property is the maximum number of carbon atoms belonging to four 
  /// separate rings that will be tolerated in a candidate molecule; it is the second argument to
  /// the Grid::secondary_deletion() method. The value of this property is read from the XML
  /// parameter <CODE>FourRingCarbonAtoms</CODE>.
  int nc4rings = 0;
  /// This non-negative integer property is the maximum number of rings that will be tolerated in 
  /// a candidate molecule; it is the third argument to the Grid::secondary_deletion() method. The 
  /// value of this property is read from the XML parameter <CODE>NumberRings</CODE>. 
  int nrings = 4;
  /// This positive integer property is the maximum number of secondary deletion attempts; it is 
  /// the fourth (and final) argument to the Grid::secondary_deletion() method. The value of this 
  /// property is read from the XML parameter <CODE>MaximumSecondary</CODE>.
  int max_secondary = 100;
  /// This non-negative integer property is the minimum number of rings that are acceptable for a 
  /// candidate molecule; it is the second argument of the Grid::rationalize() method. The value 
  /// of this property is read from the XML parameter <CODE>MinimumRings</CODE>.
  int min_rings = 1;
  /// This non-negative integer property is the maximum number of rings that are acceptable for a 
  /// candidate molecule; it is the third argument of the Grid::rationalize() method. The value of
  /// this property must be greater than or equal to that of Molecular_Assembler::min_rings and it
  /// is read from the XML parameter <CODE>MaximumRings</CODE>.
  int max_rings = 6;
  /// This positive integer property sets the number of iterations for the loop in the run() method 
  /// over the call of the Grid::initial_deletion() method. The value of this property is read from the 
  /// XML parameter <CODE>NumberInitial</CODE>.
  int n_initial = 3;
  /// This positive integer property sets the number of iterations for the loop in the run() method 
  /// over the call of the Grid::secondary_deletion() method. The value of this property is read from the 
  /// XML parameter <CODE>NumberSecondary</CODE>.
  int n_secondary = 3;
  /// This positive integer property sets the number of iterations for the loop in the run() method 
  /// over the call of the Grid::path_selection() method. The value of this property is read from the 
  /// XML parameter <CODE>NumberPath</CODE>.
  int n_path = 3;
  /// This positive integer property sets the number of iterations for the loop in the run() method 
  /// over the call of the Grid::rationalize() method. The value of this property is read from the 
  /// XML parameter <CODE>NumberRationalize</CODE>.
  int n_rationalize = 3;
  /// This positive integer property sets the number of iterations for the loop in the run() method 
  /// over the call of the Grid::add_hydrogens() method. The value of this property is read from the 
  /// XML parameter <CODE>NumberDesaturate</CODE>.
  int n_desaturate = 1;
  /// This non-negative integer property is the number of pharmacophoric nodes 
  /// in the grid and is the final argument in the Grid constructor. The value 
  /// of this property is read from the XML parameter <CODE>NumberPharmacophores</CODE>.
  int npharmacophore = 3;
  /// This positive integer property is the row number of this ensemble of 
  /// parameters in the <CODE>Parameter_Set</CODE> table of the SQLite database. 
  int parameter_id = 0;
  /// This positive integer property is the number of C++11 threads that 
  /// will be used to run this software in parallel; its value is read from 
  /// the XML parameter <CODE>NumberThreads</CODE>. 
  int nthread = 1;
  /// The dimension of the grid used to build the hydrocarbon skeleton of the molecule; the value 
  /// of this property sets the size of the \f$x\f$ and \f$y\f$ dimensions, while the \f$z\f$ dimension 
  /// has the value \f$ N/2 + q\f$ where \f$q \equiv N \pmod {N/2}\f$ and \f$N\f$ is the grid size.
  /// The value of this property, which must be positive, is read from the XML parameter <CODE>GridSize</CODE>.  
  int grid_size = 17;
  /// The Unix process ID for this instance of the Molecular_Assembler, which is used to identify 
  /// the binary file in which the output molecules are written.
  int process_id = -1;
  /// This positive integer property is the number of molecules that will be 
  /// created, after which the program exits; its value is read from 
  /// the XML parameter <CODE>NumberMolecules</CODE>. 
  unsigned long n_mols = 50000;

  /// This Boolean property controls whether or not random paths are used 
  /// to keep pharmacophoric nodes connected in the candidate molecule; it 
  /// is the unique argument to the Grid::path_selection() method. The value 
  /// of this property is read from the XML parameter <CODE>PharmacophoreHardening</CODE>.
  bool pharm_hardening = true;
  /// This Boolean property controls whether or not the Molecule::decorate() 
  /// method will attempt to subsitute oxygen atoms in the molecule. The value 
  /// of this property is read from the XML parameter <CODE>SubstituteOxygen</CODE>.
  bool subs_oxygen = true;
  /// This Boolean property controls whether or not the Molecule::decorate() 
  /// method will attempt to subsitute nitrogen atoms in the molecule. The value 
  /// of this property is read from the XML parameter <CODE>SubstituteNitrogen</CODE>.
  bool subs_nitrogen = true;
  /// This Boolean property controls whether or not the Molecule::decorate() 
  /// method will attempt to subsitute sulfur atoms in the molecule. The value 
  /// of this property is read from the XML parameter <CODE>SubstituteSulfur</CODE>.
  bool subs_sulfur = false;
  /// This Boolean property controls whether or not the Molecule::decorate() 
  /// method will attempt to subsitute halogens for terminal hydrogen atoms
  /// in the molecule. The value of this property is read from the XML parameter 
  /// <CODE>SubstituteFunctionalGroups</CODE>.
  bool subs_functional = false;
  /// This Boolean property controls whether or not the Molecule::decorate() 
  /// method will attempt to create amide groups in the molecule. The value 
  /// of this property is read from the XML parameter <CODE>CreateExotic</CODE>.
  bool create_exotic = false;
  /// This Boolean property controls whether or not the Molecule::decorate() 
  /// method will attempt to convert hexa-atomic rings to penta-atomic rings 
  /// in the molecule. The value of this property is read from the XML parameter 
  /// <CODE>CreateFiveMemberRings</CODE>.
  bool create_penta = true;
  /// This Boolean property controls whether or not the Molecule::decorate() 
  /// method will attempt to convert C-C bonds to C=C bonds in the molecule. 
  /// The value of this property is read from the XML parameter <CODE>CreateDoubleBonds</CODE>.
  bool create_double = true;
  /// This Boolean property controls whether or not the Molecule::decorate() 
  /// method will attempt to convert C-C bonds to C≡C bonds in the molecule.  
  /// The value of this property is read from the XML parameter <CODE>CreateTripleBonds</CODE>. 
  bool create_triple = false;
  /// This Boolean property controls whether or not the Molecule::decorate() 
  /// method will attempt to eliminate all of the axial neighbours of ring 
  /// atoms in the molecule, which promotes the creation of aromatic rings.
  /// The value of this property is read from the XML parameter <CODE>StripAxialMethyls</CODE>.
  bool kill_axial = true;
  
  /// This string property is the name of the SQLite database file; its value
  /// is read from the XML parameter <CODE>DatabaseFile</CODE>. 
  std::string database = "";

  /// This floating point property is the percentage of carbon atoms 
  /// that should be removed from the initial grid volume and as such 
  /// must lie between zero and one; it is the first argument to the 
  /// Grid::initial_deletion() method. The value of this property is read 
  /// from the XML parameter <CODE>InitialPercentage</CODE>.
  double percent = 0.5;
  /// This floating point property is the percentage of CH3 groups 
  /// that should be trimmed from the molecular scaffold and as such 
  /// it must lie between zero and one; it is the first argument of 
  /// the Grid::rationalize() method. The value of this property is read 
  /// from the XML parameter <CODE>PercentMethyl</CODE>.
  double percent_methyl = 0.4;
  /// This positive floating point property is the carbon-carbon bond 
  /// length (measured in Angstroms) that is used in assembling the grid, 
  /// where is the inter-node distance; it is the fourth argument to the 
  /// Grid constructor. The value of this property is read from the 
  /// XML parameter <CODE>BondLength</CODE>.
  double bond_length = 1.52;
  /// This positive floating point property is the total volume of the 
  /// grid and thus sets the maximum dimensions of any candidate molecule; 
  /// it is the unique argument of the Grid::blank_pharmacophore() method. 
  /// The value of this property is read from the XML parameter 
  /// <CODE>PharmacophoreRadius</CODE>.
  double pharmacophore_radius = 3.5;

  /// This method builds a string - the method's argument - containing all of the parameter values (i.e. the properties of this class) in a format appropriate for an SQL statement.
  void create_parameter_string(std::string&) const;
  /// This method creates the SQLite database if it doesn't exist, creating the two empty tables Parameter_Set and Compound with the appropriate columns. 
  void create_database() const;
  /// This method uses the Grid and Molecule classes to build molecules, storing them in a binary file in the scratch directory. The first argument is the thread number and the second the program's process ID, needed to name the binary molecule file. 
  void run(int,unsigned long) const;  
 public:
  /// The constructor for this class, which accepts as its unique argument the name of the parameter file that it will parse to obtain the values for its properties.
  Molecular_Assembler(const std::string&);
  /// The principal method for this class, called by the C++ main() program after initializing an instance of this class; it creates the "scratch" directory for storing the binary molecule files, calls the run() method and then reads the binary molecule files to write their content to the SQLite database.
  void assemble();
};
#endif
