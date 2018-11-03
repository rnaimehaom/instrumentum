#include "grid.h"

#ifndef _masmh
#define _masmh

class Molecular_Assembler {
 private:
  int max_attempts = 100;
  int nrings = 4;
  int nc4 = 2;
  int nc4rings = 0;
  int min_rings = 1;
  int max_rings = 6;
  int n_initial = 3;
  int n_secondary = 3;
  int n_path = 3;
  int max_secondary = 100;
  int n_rationalize = 3;
  int n_desaturate = 1;
  int npharmacophore = 3;
  int parameter_id = 0;
  int nthread = 0;

  long n_mols = 50000;
  long seed = 0;    

  bool pharm_hardening = true;
  bool subs_oxygen = true;
  bool subs_nitrogen = true;
  bool subs_sulfur = false;
  bool subs_functional = false;
  bool create_exotic = false;
  bool create_penta = true;
  bool create_double = true;
  bool create_triple = false;
  bool kill_axial = true;
  
  std::string database = "";

  double percent = 0.5;
  double percent_methyl = 0.4;
  double bond_length = 1.52;
  double pharmacophore_radius = 3.5;

  void database_insertion(const std::string&,const std::string&,sqlite3*) const;
  void create_parameter_string(std::string&) const;
  void create_database() const;  
 public:
  Molecular_Assembler(const std::string&);
  void run() const;
};
#endif
