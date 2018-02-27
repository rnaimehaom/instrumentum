#include "molecular_assembler.h"

Random RND;

Molecular_Assembler::Molecular_Assembler()
{
  set_default_values();
}

void Molecular_Assembler::set_default_values()
{
  percent = 0.5;
  percent_methyl = 0.4;
  bond_length = 1.52;
  pharmacophore_radius = 3.5;

  max_attempts = 100;
  max_secondary = 100;
  npharmacophore = 3;
  nc4 = 2;
  nc4rings = 0;
  nrings = 4;
  min_rings = 1;
  max_rings = 6;
  n_initial = 3;
  n_secondary = 3;
  n_path = 3;
  n_rationalize = 3;
  n_desaturate = 1;
  n_mols = 50000;

  pharm_hardening = true;
  subs_oxygen = true;
  subs_nitrogen = true;
  subs_sulfur = false;
  subs_functional = false;
  create_double = true;
  create_triple = false;
  create_penta = true;
  create_exotic = false;
  kill_axial = true;

  database = "";
  parameter_id = 0;
  seed = 0;
  nthread = 1;
}

Molecular_Assembler::Molecular_Assembler(const char* filename)
{
  // This method reads in the parameters from a file, the name 
  // of which is passed in as an argument to the method
  // The file format should be '#' for a comment, otherwise 
  // varname = value
  unsigned int i,n,bpoint;
  std::string line,name,value,parameter_string;

  set_default_values();

  // Open the file
  std::ifstream s(filename,std::ios_base::in);
  if (!s.is_open()) {
    // File doesn't exist, print an error message and die
    std::cout << "The file " << filename << " cannot be found!" << std::endl;
    std::exit(1);
  }
  // Loop through all lines in the parameter file
  while(std::getline(s,line)) {
    // If it's an empty line, continue
    if (line.empty()) continue;
    // If the line begins with a #, ignore it
    if (line[0] == '#') continue;
    // If there's no equals sign in this line, continue
    if (line.find('=') == std::string::npos) continue;
    n = line.size();
    name = "";
    for(i=0; i<n; ++i) {
      if (line[i] == ' ') continue;
      if (line[i] == '=') {
        bpoint = i;
        break;
      }
      name += line[i];
    }
    value = "";
    for(i=bpoint+1; i<n; ++i) {
      if (line[i] == ' ') continue;
      value += line[i];
    }
    // Now that we have the parameter name, see if it matches
    // any of the known parameters. If so, read in the value and
    // assign it
    if (name == "PharmacophoreRadius") {
      pharmacophore_radius = stod(value);
    }
    else if (name == "InitialPercentage") {
      percent = stod(value);
    }
    else if (name == "PercentMethyl") {
      percent_methyl = stod(value);
    }
    else if (name == "BondLength") {
      bond_length = stod(value);
    }
    else if (name == "DatabaseFile") {
      database = value;
    }
    else if (name == "RandomSeed") {
      seed = stol(value);
    }
    else if (name == "MaximumAttempts") {
      max_attempts = stoi(value);
    }
    else if (name == "NumberRings") {
      nrings = stoi(value);
    }
    else if (name == "NumberC4Atoms") {
      nc4 = stoi(value);
    }
    else if (name == "NumberC4Rings") {
      nc4rings = stoi(value);
    }
    else if (name == "MinimumRings") { 
      min_rings = stoi(value);
    }
    else if (name == "MaximumRings") {
      max_rings = stoi(value);
    }
    else if (name == "NumberInitial") { 
      n_initial = stoi(value);
    }
    else if (name == "NumberSecondary") {
      n_secondary = stoi(value);
    }
    else if (name == "NumberPath") {
      n_path = stoi(value);
    }
    else if (name == "MaximumSecondary") {
      max_secondary = stoi(value);
    }
    else if (name == "NumberRationalize") {
      n_rationalize = stoi(value);
    }
    else if (name == "NumberDesaturate") {
      n_desaturate = stoi(value);
    }
    else if (name == "NumberMolecules") {
      n_mols = stoi(value);
    }
    else if (name == "NumberPharmacophores") {
      npharmacophore = stoi(value);
    }
    else if (name == "CreateFiveMemberRings") {
      capitalize(value);
      create_penta = (value == "YES") ? true : false;
    }
    else if (name == "SubstituteFunctionalGroups") {
      capitalize(value);
      subs_functional = (value == "YES") ? true : false;
    }
    else if (name == "CreateDoubleBonds") { 
      capitalize(value);
      create_double = (value == "YES") ? true : false;
    }
    else if (name == "CreateTripleBonds") {
      capitalize(value);
      create_triple = (value == "YES") ? true : false;
    }
    else if (name == "CreateExotic") {
      capitalize(value);
      create_exotic = (value == "YES") ? true : false;
    }
    else if (name == "PharmacophoreHardening") {
      capitalize(value);
      pharm_hardening = (value == "YES") ? true : false;
    }
    else if (name == "SubstituteOxygen") {
      capitalize(value);
      subs_oxygen = (value == "YES") ? true : false;
    }
    else if (name == "SubstituteNitrogen") {
      capitalize(value);
      subs_nitrogen = (value == "YES") ? true : false;
    }
    else if (name == "SubstituteSulfur") {
      capitalize(value);
      subs_sulfur = (value == "YES") ? true : false;
    }
    else if (name == "StripAxialMethyls") {
      capitalize(value);
      kill_axial = (value == "YES") ? true : false;
    }
  }
  s.close();
  // Sanity checks...
  assert(n_mols > 0);
  assert(n_rationalize > 0);
  assert(n_desaturate > 0);
  assert(n_path > 0);
  assert(n_secondary > 0);
  assert(n_initial > 0);
  assert(min_rings >= 0);
  assert(max_rings >= min_rings);
  assert(max_attempts > 0);
  assert(nc4 >= 0);
  assert(nc4rings >= 0);
  assert(bond_length > std::numeric_limits<double>::epsilon());
  assert(pharmacophore_radius > std::numeric_limits<double>::epsilon());
  assert(percent_methyl > std::numeric_limits<double>::epsilon() && percent_methyl < 1.0);
  assert(percent > std::numeric_limits<double>::epsilon() && percent < 1.0);

  if (seed == 0) seed = (unsigned long) std::time(NULL);

#ifdef _OPENMP
  nthread = atoi(std::getenv("OMP_NUM_THREADS"));
#endif

  // Check if the database exists, if not the tables need to be created!
  if (!boost::filesystem::exists(database)) create_database();

  sqlite3* dbase;
  sqlite3_open(database.c_str(),&dbase);
  create_parameter_string(parameter_string);
  std::string query = "INSERT INTO Parameter_Set(";
  query += "max_attempts,";
  query += "nrings,";
  query += "nc4,";
  query += "nc4rings,";
  query += "min_rings,";
  query += "max_rings,";
  query += "n_initial,";
  query += "n_secondary,";
  query += "n_path,";
  query += "max_secondary,";
  query += "n_rationalize,";
  query += "n_desaturate,";
  query += "npharmacophore,";
  query += "pharm_hardening,";
  query += "subs_oxygen,";
  query += "subs_sulfur,";
  query += "subs_nitrogen,";
  query += "subs_functional,";
  query += "create_exotic,";
  query += "create_penta,";
  query += "create_double,";
  query += "create_triple,";
  query += "kill_axial,";
  query += "percent,";
  query += "percent_methyl,";
  query += "bond_length,";
  query += "pharmacophore_radius,";
  query += "rng_seed,";
  query += "nthread,";
  query += "timestamp) ";
  query +=  "VALUES " + parameter_string + ";";
  sqlite3_exec(dbase,query.c_str(),NULL,NULL,NULL);
  // We need to get the value of the parameter_id from this INSERT statement!
  parameter_id = int(sqlite3_last_insert_rowid(dbase));
  sqlite3_close(dbase);
}

Molecular_Assembler::~Molecular_Assembler()
{

}

void Molecular_Assembler::create_database() const
{
  std::string query; 
  sqlite3* dbase;

  sqlite3_open(database.c_str(),&dbase);

  query = "CREATE TABLE Parameter_Set("; 
  query += "max_attempts INTEGER,";
  query += "nrings INTEGER,";
  query += "nc4 INTEGER,";
  query += "nc4rings INTEGER,";
  query += "min_rings INTEGER,";
  query += "max_rings INTEGER,";
  query += "n_initial INTEGER,";
  query += "n_secondary INTEGER,";
  query += "n_path INTEGER,";
  query += "max_secondary INTEGER,";
  query += "n_rationalize INTEGER,";
  query += "n_desaturate INTEGER,";
  query += "npharmacophore INTEGER,";
  query += "pharm_hardening BOOLEAN,";
  query += "subs_oxygen BOOLEAN,";
  query += "subs_sulfur BOOLEAN,";
  query += "subs_nitrogen BOOLEAN,";
  query += "subs_functional BOOLEAN,";
  query += "create_exotic BOOLEAN,";
  query += "create_penta BOOLEAN,";
  query += "create_double BOOLEAN,";
  query += "create_triple BOOLEAN,";
  query += "kill_axial BOOLEAN,";
  query += "percent REAL,";
  query += "percent_methyl REAL,";
  query += "bond_length REAL,";
  query += "pharmacophore_radius REAL,";
  query += "rng_seed INTEGER,";
  query += "nthread INTEGER,";
  query += "timestamp DATETIME,";
  query += "parameter_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL);";
  sqlite3_exec(dbase,query.c_str(),NULL,NULL,NULL);

  query = "CREATE TABLE Compound(";
  query += "raw_structure TEXT,";
  query += "op_string TEXT,";
  query += "minimized_structure TEXT,";
  query += "root_mean_square REAL,";
  query += "synthetic_feasibility REAL,";
  query += "parameter_id INTEGER,";
  query += "compound_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,";
  query += "FOREIGN KEY(parameter_id) REFERENCES Parameter_Set(parameter_id));";
  sqlite3_exec(dbase,query.c_str(),NULL,NULL,NULL);

  sqlite3_close(dbase);
}

void Molecular_Assembler::create_parameter_string(std::string& output) const
{
  std::ostringstream s; 

  s << "(" << max_attempts << ",";
  s << nrings << ",";
  s << nc4 << ",";
  s << nc4rings << ",";
  s << min_rings << ",";
  s << max_rings << ",";
  s << n_initial << ",";
  s << n_secondary << ",";
  s << n_path << ",";
  s << max_secondary << ",";
  s << n_rationalize << ",";
  s << n_desaturate << ",";
  s << npharmacophore << ",";

  if (pharm_hardening) {
    s << "1,";
  }
  else {
    s << "0,";
  } 
  if (subs_oxygen) {
    s << "1,";
  }
  else {
    s << "0,";
  } 
  if (subs_sulfur) {
    s << "1,";
  }
  else {
    s << "0,";
  } 
  if (subs_nitrogen) {
    s << "1,";
  }
  else {
    s << "0,";
  } 
  if (subs_functional) {
    s << "1,";
  }
  else {
    s << "0,";
  } 
  if (create_exotic) {
    s << "1,";
  }
  else {
    s << "0,";
  } 
  if (create_penta) {
    s << "1,";
  }
  else {
    s << "0,";
  } 
  if (create_double) {
    s << "1,";
  }
  else {
    s << "0,";
  } 
  if (create_triple) {
    s << "1,";
  }
  else {
    s << "0,";
  } 
  if (kill_axial) {
    s << "1,";
  }
  else {
    s << "0,";
  } 
  s << percent << ",";
  s << percent_methyl << ",";
  s << bond_length << ",";
  s << pharmacophore_radius << ",";
  s << seed << ",";
  s << nthread << ",";
  s << "\'DATETIME(\'NOW\')\')";
  output = s.str();
}

void Molecular_Assembler::run() const
{
  sqlite3* dbase;
  unsigned int mol_created = 0;
  bool ornaments[] = {kill_axial,create_penta,create_double,create_triple,create_exotic,subs_oxygen,subs_sulfur,subs_nitrogen,subs_functional};

  sqlite3_open(database.c_str(),&dbase);

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+:mol_created)
  {
#endif
  unsigned int build,i,j,k,l,q;
  unsigned long s = seed;
  bool test;
  std::string mstring,ops;
  Grid* g = new Grid(bond_length,npharmacophore);
  Molecule* m = new Molecule;

#ifdef _OPENMP
  s *= (1 + omp_get_thread_num());
#endif

  RND.initialize_generator(s);

  while(mol_created < n_mols) {
#if VERBOSE
    std::cout << "Putting pharmacophore..." << std::endl;
#endif
    g->blank_pharmacophore(pharmacophore_radius);
    for(i=0; i<n_initial; ++i) {
      g->fill_interior();
      g->save_state(0);
#ifdef VERBOSE
      std::cout << "Initial deletion..." << std::endl;
#endif
      test = g->initial_deletion(percent,max_attempts);
      if (!test) { 
        g->restore(0);
        continue;
      }
      for(j=0; j<n_path; ++j) {
        g->save_state(1);
#ifdef VERBOSE
        std::cout << "Path hardening..." << std::endl;
#endif
        test = g->path_selection(pharm_hardening);
        if (!test) {
          g->restore(1);
          continue;
        }
        for(k=0; k<n_secondary; ++k) {
          g->save_state(2);
#ifdef VERBOSE
          std::cout << "Secondary deletion..." << std::endl;
#endif
          test = g->secondary_deletion(nc4,nc4rings,nrings,max_secondary);
          if (!test) {
            g->restore(2);
            continue;
          }
          for(l=0; l<n_rationalize; ++l) {
            g->save_state(3);
#ifdef VERBOSE
            std::cout << "Rationalizing..." << std::endl;
#endif
            test = g->rationalize(percent_methyl,min_rings,max_rings);
            if (!test) {
              g->restore(3);
              continue;
            }
            g->add_hydrogens();
            build = 0;
            for(q=0; q<n_desaturate; ++q) {
#ifdef VERBOSE
              std::cout << "Writing grid to molecule..." << std::endl;
#endif
              g->write_scaffold(m);
#ifdef VERBOSE
              std::cout << "Decorating molecule..." << std::endl;
#endif
              test = m->decorate(ornaments);
              if (test) {
                build++;
                mstring = m->to_MDLMol();
                ops = m->get_opstring();
#ifdef _OPENMP
#pragma omp critical 
        {
#endif
                database_insertion(ops,mstring,dbase);
#ifdef _OPENMP
        }
#endif
              }
              m->clear();
            }
            mol_created += build;
            g->restore(3);
          }
          g->restore(2);
        }
        g->restore(1);
      }
      g->restore(0);
    }
    g->clear();
  }
  delete g;
  delete m;
#ifdef _OPENMP
  }
#endif
  sqlite3_close(dbase);
}

void Molecular_Assembler::database_insertion(const std::string& opstring,const std::string& molecule,sqlite3* dbase) const
{
  std::string query = "INSERT INTO Compound (parameter_id,op_string,raw_structure) VALUES (";
  std::stringstream sstream;
  sstream << parameter_id;
  query += sstream.str() + ",\'" + opstring + "\',\'" + molecule + "\');";
  sqlite3_exec(dbase,query.c_str(),NULL,NULL,NULL); 
}
