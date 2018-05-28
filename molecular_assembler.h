#include "grid.h"

#ifndef _masmh
#define _masmh

class Molecular_Assembler {
 private:
  unsigned int max_attempts = 100;
  unsigned int nrings = 4;
  unsigned int nc4 = 2;
  unsigned int nc4rings = 0;
  unsigned int min_rings = 1;
  unsigned int max_rings = 6;
  unsigned int n_initial = 3;
  unsigned int n_secondary = 3;
  unsigned int n_path = 3;
  unsigned int max_secondary = 100;
  unsigned int n_rationalize = 3;
  unsigned int n_desaturate = 1;
  unsigned int n_mols = 50000;
  unsigned int npharmacophore = 3;
  unsigned int parameter_id = 0;
    
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
  unsigned long seed = 0;
  unsigned int nthread = 0;

  double percent = 0.5;
  double percent_methyl = 0.4;
  double bond_length = 1.52;
  double pharmacophore_radius = 3.5;

  void database_insertion(const std::string&,const std::string&,sqlite3*) const;
  void create_parameter_string(std::string&) const;
  void create_database() const;  
 public:
  Molecular_Assembler(const char*);
  void run() const;
};
#endif
