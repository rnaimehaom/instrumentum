#include "molecule.h"
#include "grid.h"

class MASM {
 private:
  enum VTYPE {INT,BOOL,STRING,FLOAT};

  unsigned int snumber;

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
  unsigned int npharm;
    
  bool path_hardening;
  bool subs_oxy;
  bool subs_nit;
  bool subs_sul;
  bool subs_fun;
  bool create_exotic;
  bool create_penta;
  bool create_double;
  bool create_triple;
  bool kill_axial;
  
  std::string filename;
  std::string pharm_fname;

  double percent;
  double percent_methyl;
  double bond_length;
  double pharm_radius;

  void set_default_values();

  std::string create_pstring() const;  
 public:
  MASM();
  MASM(const char*);
  void set_jobid(unsigned int);
  void retrieve_db();
  void run() const;
};
