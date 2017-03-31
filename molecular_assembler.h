#include "molecule.h"
#include "grid.h"

class Molecular_Assembler {
 private:
  unsigned int max_attempts;
  unsigned int nrings;
  unsigned int nc4;
  unsigned int nc4rings;
  unsigned int min_rings;
  unsigned int max_rings;
  unsigned int n_hardening;
  unsigned int n_initial;
  unsigned int n_secondary;
  unsigned int n_path;
  unsigned int max_secondary;
  unsigned int n_demethylate;
  unsigned int n_desaturate;
  unsigned int n_mols;
  unsigned int npharmacophore;
  unsigned int parameter_id;
    
  bool path_hardening;
  bool subs_oxygen;
  bool subs_nitrogen;
  bool subs_sulfur;
  bool subs_functional;
  bool create_exotic;
  bool create_penta;
  bool create_double;
  bool create_triple;
  bool kill_axial;
  
  std::string database;
  std::string pharmacophore_filename;
  unsigned long seed;

  double percent;
  double percent_methyl;
  double bond_length;
  double pharmacophore_radius;

  void set_default_values();
  void write2disk(const std::string&,const std::string&,sqlite3*) const;
  void create_parameter_string(std::string&) const;
  void create_database() const;  
 public:
  Molecular_Assembler();
  Molecular_Assembler(const char*);
  ~Molecular_Assembler();
  void run() const;
};
