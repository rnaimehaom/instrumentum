#include "masm.h"

MASM::MASM()
{
  set_default_values();
}

void MASM::set_default_values()
{
  snumber = 0;

  percent = 0.5;
  percent_methyl = 0.4;
  bond_length = 1.52;
  pharm_radius = 3.5;

  max_attempts = 100;
  max_secondary = 100;
  npharm = 3;
  nc4 = 2;
  nc4rings = 0;
  nrings = 4;
  min_rings = 1;
  max_rings = 6;
  n_initial = 3;
  n_secondary = 3;
  n_hardening = 3;
  n_path = 3;
  n_demethylate = 3;
  n_desaturate = 1;
  n_mols = 50000;

  path_hardening = true;
  subs_oxy = true;
  create_double = true;
  subs_nit = true;
  subs_sul = false;
  create_triple = false;
  subs_fun = false;
  create_penta = true;
  create_exotic = false;
  kill_axial = true;

  filename = "";
  pharm_fname = ""; 
}

MASM::MASM(const char* fname)
{
  // This method reads in the parameters from a file, the name 
  // of which is passed in as an argument to the method
  // The file format should be '#' for a comment, otherwise 
  // varname = value
  int i,eq_point;
  std::string line,name,value,state_file;

  set_default_values();

  // Open the file
  std::ifstream s(fname,std::ios_base::in);
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
    // Find the position of the equals sign
    eq_point = 0;
    for(i=0; i<(signed) line.length(); ++i) {
      if (line[i] == '=') {
         eq_point = i;
         break;
      }
    }
    // If there's no equals sign in this line, continue
    if (eq_point < 1) continue;
    name = line.substr(0,eq_point-1); trim(name);
    value = line.substr(eq_point+1,line.length()); trim(value);
    // Now that we have the parameter name, see if it matches
    // any of the known parameters. If so, read in the value and
    // assign it
    if (name == "pharm_radius") {
      pharm_radius =  boost::lexical_cast<double>(value);
    }
    else if (name == "percent") {
      percent =  boost::lexical_cast<double>(value);
    }
    else if (name == "percent_methyl") {
      percent_methyl =  boost::lexical_cast<double>(value);
    }
    else if (name == "bond_length") {
      bond_length = boost::lexical_cast<double>(value);
    }
    else if (name == "ofilename") {
      filename = value;
    }
    else if (name == "pfilename") {
      pharm_fname = value;
    }
    else if (name == "max_attempts") {
      max_attempts = boost::lexical_cast<int>(value);
    }
    else if (name == "nrings") {
      nrings =  boost::lexical_cast<int>(value);
    }
    else if (name == "nc4") {
      nc4 =  boost::lexical_cast<int>(value);
    }
    else if (name == "nc4rings") {
      nc4rings =  boost::lexical_cast<int>(value);
    }
    else if (name == "min_rings") { 
      min_rings =  boost::lexical_cast<int>(value);
    }
    else if (name == "max_rings") {
      max_rings =  boost::lexical_cast<int>(value);
    }
    else if (name == "n_initial") { 
      n_initial =  boost::lexical_cast<int>(value);
    }
    else if (name == "n_secondary") {
      n_secondary =  boost::lexical_cast<int>(value);
    }
    else if (name == "n_path") {
      n_path =  boost::lexical_cast<int>(value);
    }
    else if (name == "max_secondary") {
      max_secondary =  boost::lexical_cast<int>(value);
    }
    else if (name == "n_demethylate") {
      n_demethylate =  boost::lexical_cast<int>(value);
    }
    else if (name == "n_desaturate") {
      n_desaturate =  boost::lexical_cast<int>(value);
    }
    else if (name == "n_mols") {
      n_mols =  boost::lexical_cast<int>(value);
    }
    else if (name == "npharm") {
      npharm =  boost::lexical_cast<int>(value);
    }
    else if (name == "n_hardening") {
      n_hardening =  boost::lexical_cast<int>(value);
    }
    else if (name == "create_penta") {
      create_penta = boost::lexical_cast<bool>(value);
    }
    else if (name == "subs_fun") {
      subs_fun = boost::lexical_cast<bool>(value);
    }
    else if (name == "create_double") { 
      create_double = boost::lexical_cast<bool>(value);
    }
    else if (name == "create_triple") {
      create_triple = boost::lexical_cast<bool>(value);
    }
    else if (name == "create_exotic") {
      create_exotic = boost::lexical_cast<bool>(value);
    }
    else if (name == "path_hardening") {
      path_hardening = boost::lexical_cast<bool>(value);
    }
    else if (name == "subs_oxy") {
      subs_oxy = boost::lexical_cast<bool>(value);
    }
    else if (name == "subs_nit") {
      subs_nit = boost::lexical_cast<bool>(value);
    }
    else if (name == "subs_sul") {
      subs_sul = boost::lexical_cast<bool>(value);
    }
    else if (name == "kill_axial") {
      kill_axial = boost::lexical_cast<bool>(value);
    }
  }
  s.close();
#if DBASE
  // Write the parameters to the database and get the PID...
  pqxx::connection c("host=localhost dbname=alchemy user=aurifex password=xabub89");

  pqxx::work w(c);
  std::string ins = create_pstring();
  pqxx::result r = w.exec("INSERT INTO parameters VALUES " + ins + " RETURNING pid;");
  w.commit();
  snumber = boost::lexical_cast<int>(r[0][0].c_str());
#endif
}

std::string MASM::create_pstring() const
{
  std::ostringstream s; 

  s << "(" << max_attempts << ",";
  s << nrings <<",";
  s << nc4 << ",";
  s << nc4rings << ",";
  s << min_rings << ",";
  s << max_rings << ",";
  s << n_hardening << ",";
  s << n_initial << ",";
  s << n_secondary << ",";
  s << n_path << ",";
  s << max_secondary << ",";
  s << n_demethylate << ",";
  s << n_desaturate << ",";
  s << npharm << ",";

  if (path_hardening) {
    s << "TRUE,";
  }
  else {
    s << "FALSE,";
  } 
  if (subs_oxy) {
    s << "TRUE,";
  }
  else {
    s << "FALSE,";
  } 
  if (subs_sul) {
    s << "TRUE,";
  }
  else {
    s << "FALSE,";
  } 
  if (subs_nit) {
    s << "TRUE,";
  }
  else {
    s << "FALSE,";
  } 
  if (subs_fun) {
    s << "TRUE,";
  }
  else {
    s << "FALSE,";
  } 
  if (create_exotic) {
    s << "TRUE,";
  }
  else {
    s << "FALSE,";
  } 
  if (create_penta) {
    s << "TRUE,";
  }
  else {
    s << "FALSE,";
  } 
  if (create_double) {
    s << "TRUE,";
  }
  else {
    s << "FALSE,";
  } 
  if (create_triple) {
    s << "TRUE,";
  }
  else {
    s << "FALSE,";
  } 
  if (kill_axial) {
    s << "TRUE,";
  }
  else {
    s << "FALSE,";
  } 
  s << percent << ",";
  s << percent_methyl << ",";
  s << bond_length << ",";
  s << pharm_radius << ",";
  if (pharm_fname == "") {
    s << "NULL)";
  }
  else {
    s << pharm_fname << ")";
  }
  return s.str();
}

void MASM::set_jobid(unsigned int n)
{
  snumber = n;
}

void MASM::retrieve_db()
{
  // Get the parameters from the database...
  pqxx::connection c("host=localhost dbname=alchemy user=aurifex password=xabub89");
  std::string qstring = "SELECT max_attempts,nrings,nc4,nc4rings,min_rings,max_rings,n_hardening,n_initial,n_secondary,n_path,max_secondary,n_demethylate,n_desaturate,npharm,path_hardening,subs_oxy,subs_sul,subs_nit,subs_fun,create_exotic,create_penta,create_double,create_triple,kill_axial,percent,percent_methyl,bond_length,pharm_radius,pharm_fname FROM parameters WHERE pid=";
  qstring += boost::lexical_cast<std::string>(snumber);
  qstring += ";";
  pqxx::work w(c);
  pqxx::result r = w.exec(qstring.c_str());
  w.commit();

  max_attempts = boost::lexical_cast<int>(r[0][0].c_str());
  nrings = boost::lexical_cast<int>(r[0][1].c_str());
  nc4 = boost::lexical_cast<int>(r[0][2].c_str());
  nc4rings = boost::lexical_cast<int>(r[0][3].c_str());
  min_rings = boost::lexical_cast<int>(r[0][4].c_str());
  max_rings = boost::lexical_cast<int>(r[0][5].c_str());
  n_hardening = boost::lexical_cast<int>(r[0][6].c_str());
  n_initial = boost::lexical_cast<int>(r[0][7].c_str());
  n_secondary = boost::lexical_cast<int>(r[0][8].c_str());
  n_path = boost::lexical_cast<int>(r[0][9].c_str());
  max_secondary = boost::lexical_cast<int>(r[0][10].c_str());
  n_demethylate = boost::lexical_cast<int>(r[0][11].c_str());
  n_desaturate = boost::lexical_cast<int>(r[0][12].c_str());
  npharm = boost::lexical_cast<int>(r[0][13].c_str());

  path_hardening = boost::lexical_cast<bool>(r[0][14].c_str());
  subs_oxy = boost::lexical_cast<bool>(r[0][15].c_str());
  subs_sul = boost::lexical_cast<bool>(r[0][16].c_str());
  subs_nit = boost::lexical_cast<bool>(r[0][17].c_str());
  subs_fun = boost::lexical_cast<bool>(r[0][18].c_str());
  create_exotic = boost::lexical_cast<bool>(r[0][19].c_str());
  create_penta = boost::lexical_cast<bool>(r[0][20].c_str());
  create_double = boost::lexical_cast<bool>(r[0][21].c_str());
  create_triple = boost::lexical_cast<bool>(r[0][22].c_str());
  kill_axial = boost::lexical_cast<bool>(r[0][23].c_str());

  percent = boost::lexical_cast<double>(r[0][24].c_str());
  percent_methyl = boost::lexical_cast<double>(r[0][25].c_str());
  bond_length = boost::lexical_cast<double>(r[0][26].c_str());
  pharm_radius = boost::lexical_cast<double>(r[0][27].c_str());

  pharm_fname = std::string(r[0][28].c_str());
}

void MASM::run() const
{
#ifdef DBASE
  assert(snumber > 0);
#endif

  unsigned int mol_created = 0;
  bool ornaments[] = {kill_axial,create_penta,create_double,create_triple,create_exotic,subs_oxy,subs_sul,subs_nit,subs_fun};

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+:mol_created)
  {
#endif
  unsigned int build,i,j,k,l,q,seed = 1;
  bool test;
  Grid* g = new Grid(bond_length,npharm);
  Molecule* m = new Molecule;

#ifdef DBASE
  pqxx::connection c("host=localhost dbname=alchemy user=aurifex password=xabub89");
#endif

#ifdef _OPENMP
  seed = 1 + omp_get_thread_num();
#endif

  initialize_rand(seed);

  while(mol_created < n_mols) {
#if VERBOSE
    std::cout << "Putting pharmacophore..." << std::endl;
#endif
    g->blank_pharmacophore(pharm_radius);
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
      for(j=0; j<n_hardening; ++j) {
        g->save_state(1);
#ifdef VERBOSE
        std::cout << "Path hardening..." << std::endl;
#endif
        test = g->path_selection(path_hardening);
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
          for(l=0; l<n_demethylate; ++l) {
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
#ifdef DBASE
                m->write2db(snumber,c);
#else
#ifdef _OPENMP
#pragma omp critical
#endif
                m->write2sdf(filename);
#endif
                build++;
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
}
