#include "molecule.h"

int Molecule::write(std::ofstream& s) const
{
  // This method writes the Molecule properties "natoms", "opstring", "bonds", 
  // "btype" and "coords" to a binary diskfile.
  int i,j,n,bcount = 0;
  double x[3];
  Atom_Type a;

  s.write((char*)(&natoms),sizeof(int)); bcount += sizeof(int);
  n = opstring.length();
  s.write((char*)(&n),sizeof(int)); bcount += sizeof(int);
  for(i=0; i<n; ++i) {
    s.write((char*)(&opstring[i]),sizeof(char));   
  }
  bcount += n*sizeof(char);

  for(i=0; i<natoms; ++i) {
    a = atom_type[i];
    s.write((char*)(&a),sizeof(Atom_Type));
  }
  bcount += natoms*sizeof(Atom_Type);

  for(i=0; i<natoms; ++i) {
    n = locale[i];
    s.write((char*)(&n),sizeof(int));
  }
  bcount += natoms*sizeof(int);

  for(i=0; i<natoms; ++i) {
    for(j=0; j<4; ++j) {
      n = bonds[4*i+j];
      s.write((char*)(&n),sizeof(int));
    }
  }
  bcount += 4*natoms*sizeof(int);

  for(i=0; i<natoms; ++i) {
    for(j=0; j<4; ++j) {
      n = btype[4*i+j];
      s.write((char*)(&n),sizeof(int)); 
    }
  }
  bcount += 4*natoms*sizeof(int);

  for(i=0; i<natoms; ++i) {
    x[0] = coords[3*i];
    x[1] = coords[3*i+1];
    x[2] = coords[3*i+2];
    s.write((char*)(&x[0]),3*sizeof(double));
  }
  bcount += 3*natoms*sizeof(double);

  return bcount;
}

int Molecule::read(std::ifstream& s)
{
  int i,j,n,bcount = 0;
  char c;
  double x[3];
  Atom_Type a;

  clear();

  s.read((char*)(&natoms),sizeof(int)); bcount += sizeof(int);
  s.read((char*)(&n),sizeof(int)); bcount += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&c),sizeof(char)); 
    opstring += c; 
  }
  bcount += n*sizeof(char);

  for(i=0; i<natoms; ++i) {
    s.read((char*)(&a),sizeof(Atom_Type));
    atom_type.push_back(a);
  }
  bcount += natoms*sizeof(Atom_Type);

  for(i=0; i<natoms; ++i) {
    s.read((char*)(&n),sizeof(int));
    locale.push_back(n);
  }
  bcount += natoms*sizeof(int);

  for(i=0; i<natoms; ++i) {
    for(j=0; j<4; ++j) {
      s.read((char*)(&n),sizeof(int));
      bonds.push_back(n);
    }
  }
  bcount += 4*natoms*sizeof(int);

  for(i=0; i<natoms; ++i) {
    for(j=0; j<4; ++j) {
      s.read((char*)(&n),sizeof(int));
      btype.push_back(n); 
    }
  }
  bcount += 4*natoms*sizeof(int);

  for(i=0; i<natoms; ++i) {
    s.read((char*)(&x[0]),3*sizeof(double));
    coords.push_back(x[0]);
    coords.push_back(x[1]);
    coords.push_back(x[2]);
  }
  bcount += 3*natoms*sizeof(double);

  return bcount;  
}

std::ostream& operator <<(std::ostream& s,const Molecule& source)
{
  int i,j;
  for(i=0; i<source.natoms; ++i) {
    s << i << "  " << Molecule::atom_names[static_cast<int>(source.atom_type[i])] << "  " << source.locale[i] << std::endl;
    s << "       " << source.coords[3*i] << "  " << source.coords[3*i+1] << "  " << source.coords[3*i+2] << std::endl;
  }
  for(i=0; i<source.natoms; ++i) {
    for(j=0; j<4; ++j) {
      if (source.bonds[4*i+j] != -1) {
        s << i << "  " << source.bonds[4*i+j] << "  " << source.btype[4*i+j] << std::endl;
      }
    }
  }
  for(i=0; i<source.nrings; ++i) {
    for(j=0; j<6; ++j) {
      s << source.rings[6*i+j] << "  ";
    }
    s << std::endl;
  }
  return s;
}

std::string Molecule::to_MDLMol() const
{
  int i,j,bnumber = 0;
  long seconds = long(std::time(nullptr));
  std::string atom;
  std::ostringstream s;

#ifdef DEBUG
  assert(consistent());
#endif

  s << "mol_" << opstring << "_" << seconds << std::endl;
  s << "  MOE2000" << std::endl;
  s << std::endl;
  for(i=0; i<natoms; ++i) {
    for(j=0; j<4; ++j) {
      if (bonds[4*i+j] > i) bnumber++;
    }
  }
  s << std::setw(3) << natoms << std::setw(3) << bnumber << " 0  0  0  0  0  0  0  0   1 V2000" << std::endl;
  for(i=0; i<natoms; ++i) {
    atom = Molecule::atom_names[static_cast<int>(atom_type[i])]; 
    s << std::setw(10) << std::setprecision(4) << std::setiosflags(std::ios::fixed) << coords[3*i] << std::setw(10) << std::setprecision(4) << std::setiosflags(std::ios::fixed) << coords[3*i+1] << std::setw(10) << std::setprecision(4) << std::setiosflags(std::ios::fixed) << coords[3*i+2] << " " << std::setw(3) << std::setiosflags(std::ios::left) << atom << std::resetiosflags(std::ios::left) << " 0  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
  }
  for(i=0; i<natoms; ++i) {
    for(j=0; j<4; ++j) {
      if (bonds[4*i+j] > i) {
        s << std::setw(3) << (i+1) << std::setw(3) << (1+bonds[4*i+j]) << std::setw(3) << btype[4*i+j] << "  0  0  0  0  0" << std::endl;
      }
    }
  }
  s << "M  END" << std::endl;
  s << "$$$$" << std::endl;
  return s.str();
}
