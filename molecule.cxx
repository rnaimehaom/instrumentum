#include "molecule.h"

// Molecular operations...
const char Molecule::ops[7] = {'P','N','S','O','T','F','A'};
// Element symbols...
const std::string Molecule::atom_names[12] = {"","C","H","N","O","F","P","S","Cl","Br","Ag","I"};

Molecule::Molecule(unsigned long s)
{
  if (s == 0) s = (unsigned) std::time(nullptr);
  RND.initialize_generator(s);
}

Molecule::Molecule(const Molecule& source)
{
  if (p_allocated > 0) {
    delete[] pieces;
    p_allocated = 0;
  }
  natoms = source.natoms;
  atom_type = source.atom_type;
  locale = source.locale;
  bonds = source.bonds;
  btype = source.btype;
  rbonds = source.rbonds;
  rings = source.rings;
  coords = source.coords;
  opstring = source.opstring;
  if (source.p_allocated > 0) {
    p_allocated = source.p_allocated;
    pieces = new std::vector<int>[p_allocated];
    for(int i=0; i<p_allocated; ++i) {
      pieces[i] = source.pieces[i];
    }
  }
}

Molecule::~Molecule()
{
  if (p_allocated > 0) delete[] pieces;
}

Molecule& Molecule::operator =(const Molecule& source)
{
  if (this == &source) return *this;

  if (p_allocated > 0) {
    delete[] pieces;
    p_allocated = 0;
  }
  natoms = source.natoms;
  atom_type = source.atom_type;
  locale = source.locale;
  bonds = source.bonds;
  btype = source.btype;
  rbonds = source.rbonds;
  rings = source.rings;
  coords = source.coords;
  opstring = source.opstring;
  if (source.p_allocated > 0) {
    p_allocated = source.p_allocated;
    pieces = new std::vector<int>[p_allocated];
    for(int i=0; i<p_allocated; ++i) {
      pieces[i] = source.pieces[i];
    }
  }

  return *this;
}

void Molecule::clear()
{
  nrings = 0;
  natoms = 0;
  bonds.clear();
  btype.clear();
  atom_type.clear();
  locale.clear();
  rings.clear();
  rbonds.clear();
  coords.clear();
  opstring = "";
  if (p_allocated > 0) delete[] pieces;
  p_allocated = 0;
}

void Molecule::add_atom(const Atom_Type& atomic_number,const double* x,int ltype)
{
#ifdef DEBUG
  assert(atomic_number != Atom_Type::empty);
#endif
  natoms++;
  atom_type.push_back(atomic_number);
  locale.push_back(ltype);
  coords.push_back(x[0]);
  coords.push_back(x[1]);
  coords.push_back(x[2]);
  for(int i=0; i<4; ++i) {
    bonds.push_back(-1);
    btype.push_back(-1);
  }
}

bool Molecule::drop_atom(int n)
{
#ifdef DEBUG
  assert(n >= 0 && n < natoms);
#endif
  int i,j,neg,nafter,fminus,nminus;  
  std::vector<int> nbonds,nbtype,nlocal;
  std::vector<double> ncoords;
  std::vector<Atom_Type> ntype;

  for(i=0; i<n; ++i) {
    nlocal.push_back(locale[i]);
    ntype.push_back(atom_type[i]);
  }
  for(i=1+n; i<natoms; ++i) {
    nlocal.push_back(locale[i]);
    ntype.push_back(atom_type[i]);
  }
  for(i=0; i<4*n; ++i) {
    nbonds.push_back(bonds[i]);
    nbtype.push_back(btype[i]);
  }
  for(i=4*(n+1); i<4*natoms; ++i) {
    nbonds.push_back(bonds[i]);
    nbtype.push_back(btype[i]);
  }
  for(i=0; i<3*n; ++i) {
    ncoords.push_back(coords[i]);
  }
  for(i=3*(n+1); i<3*natoms; ++i) {
    ncoords.push_back(coords[i]);
  }

  // Make the change
  atom_type = ntype;
  locale = nlocal;
  bonds = nbonds;
  btype = nbtype;
  coords = ncoords;
  natoms--;

  // Now the re-indexing of the bond vector...
  for(i=0; i<natoms; ++i) {
    for(j=0; j<4; ++j) {
      if (bonds[4*i+j] == n) {
        bonds[4*i+j] = -1;
        btype[4*i+j] = -1;
      }
      else if (bonds[4*i+j] > n) {
        bonds[4*i+j] -= 1;
      }
    }
  }
  // And the ring list.
  for(i=0; i<nrings; ++i) {
    for(j=0; j<6; ++j) {
      if (rings[6*i+j] == n) {
        rings[6*i+j] = -1;
      }
      else if (rings[6*i+j] > n) {
        rings[6*i+j] -= 1;
      }
    }
  }  

  // Stick the "-1" at the end of bonds and rings
  for(i=0; i<natoms; ++i) {
    do {
      nafter = 0;
      nminus = 0;
      fminus = -1;
      for(j=0; j<4; ++j) {
        if (bonds[4*i+j] == -1 && fminus == -1) {
          fminus = j;
          nminus++;
        }
        else if (bonds[4*i+j] == -1 && fminus >= 0 && nafter == 0) {
          nminus++;
        }
        else if (bonds[4*i+j] >= 0 && fminus >= 0) {
          nafter++;
        }
      }
      if (nafter > 0) {
        for(j=0; j<nafter; ++j) {
          bonds[4*i+fminus+j] = bonds[4*i+fminus+nminus+j];
          btype[4*i+fminus+j] = btype[4*i+fminus+nminus+j];
        }
        for(j=0; j<nminus; ++j) {
          bonds[4*i+fminus+nafter+j] = -1;
          btype[4*i+fminus+nafter+j] = -1;
        }
      }
      else {
        break;
      }
    } while(true);
  }
  for(i=0; i<nrings; ++i) {
    neg = -1;
    for(j=0; j<6; ++j) {
      if (rings[6*i+j] == -1) neg = j;
    }
    if (neg >= 0 && neg < 5) {
      for(j=neg; j<5; ++j) {
        rings[6*i+j] = rings[6*i+j+1];
      }
      rings[6*i+5] = -1;
    }
  }
  return true;
}

void Molecule::add_bond(int atom1,int atom2,int valence)
{
#ifdef DEBUG
  assert(atom1 >= 0 && atom1 < natoms);
  assert(atom2 >= 0 && atom2 < natoms);
#endif
  int i,next;
  bool done;

  done = false;
  for(i=0; i<4; ++i) {
    if (bonds[4*atom1+i] == atom2) {
      done = true;
      break;
    }
    else if (bonds[4*atom1+i] == -1) {
      next = i;
      break;
    }
  }
  if (!done) {
    bonds[4*atom1+next] = atom2;
    btype[4*atom1+next] = valence;
  }

  done = false;
  for(i=0; i<4; ++i) {
    if (bonds[4*atom2+i] == atom1) {
      done = true;
      break;
    }
    else if (bonds[4*atom2+i] == -1) {
      next = i;
      break;
    }
  }
  if (!done) {
    bonds[4*atom2+next] = atom1;
    btype[4*atom2+next] = valence;
  }
}

bool Molecule::eliminate_atoms(std::set<int>& kill) 
{
  if (kill.empty()) return true;
  if (*(kill.begin()) < 0 || *(kill.rbegin()) >= natoms) return false;
  std::set<int>::reverse_iterator rit;

  for(rit=kill.rbegin(); rit!=kill.rend(); ++rit) {
    if (!drop_atom(*rit)) return false;
  }
  return true;
}

bool Molecule::eliminate_atoms(std::set<int>& kill,std::vector<int>& axial) 
{
  if (kill.empty()) return true;
  if (*(kill.begin()) < 0 || *(kill.rbegin()) >= natoms) return false;
  int i,n;
  std::set<int>::reverse_iterator rit;

  for(rit=kill.rbegin(); rit!=kill.rend(); ++rit) {
    n = *rit;
    for(i=0; i<6*nrings; ++i) {
      if (axial[i] < n) continue;
      axial[i] = (axial[i] == n) ? -1 : axial[i] - 1;
    }
    if (!drop_atom(n)) return false;
  }
  return true;
}

bool Molecule::saturation_check() const
{
  int i,j,k,nz,bcount;
  bool found;

  for(i=0; i<natoms; ++i) {
    if (atom_type[i] == Atom_Type::hydrogen) {
      bcount = 0;
      for(j=0; j<natoms; ++j) {
        found = false;
        for(k=0; k<4; ++k) {
          if (bonds[4*j+k] == i) {
            found = true;
            break;
          }
        }
        if (found) bcount++;
      }

      if (bcount != 1) {
#ifdef VERBOSE
        std::cout << "Hydrogen atom has excessive bond count: " << bcount << std::endl;
#endif
        return false;
      }
      nz = 0;
      for(j=0; j<4; ++j) {
        if (bonds[4*i+j] == -1) {
          nz = j;
          break;
        }
      }
      if (nz != 1) {
#ifdef VERBOSE
        std::cout << "Hydrogen atom has valence " << nz << std::endl;
#endif
        return false;
      }
    }
    else if (atom_type[i] == Atom_Type::carbon) {
      bcount = 0;
      for(j=0; j<natoms; ++j) {
        found = false;
        for(k=0; k<4; ++k) {
          if (bonds[4*j+k] == i) {
            found = true;
            break;
          }
        }
        if (found) bcount++;
      }
      if (bcount != 4) {
#ifdef VERBOSE
        std::cout << "Carbon atom has excessive bond count: " << bcount << std::endl;
#endif
        return false;
      }
      nz = 4;
      for(j=0; j<4; ++j) {
        if (bonds[4*i+j] == -1) {
          nz = 1 + j;
          break;
        }
      }
      if (nz != 4) {
#ifdef VERBOSE
        std::cout << "Carbon atom has valence " << nz << std::endl;
#endif
        return false;
      }
    }
  }
  return true;
}

bool Molecule::valence_check() const
{
  int i,j,bcount;

  for(i=0; i<natoms; ++i) {
    bcount = 0;
    for(j=0; j<4; ++j) {
      if (btype[4*i+j] > 0) bcount += btype[4*i+j];
    }
    if (atom_type[i] == Atom_Type::hydrogen && bcount != 1) {
#ifdef VERBOSE
      std::cout << "Hydrogen valence error for " << i << " with " << bcount << std::endl;
#endif
      return false;
    }
    else if (atom_type[i] == Atom_Type::fluorine && bcount != 1) {
#ifdef VERBOSE
      std::cout << "Fluorine valence error for " << i << " with " << bcount << std::endl;
#endif
      return false;
    }
    else if (atom_type[i] == Atom_Type::phosphorus && bcount != 3) {
#ifdef VERBOSE
      std::cout << "Phosphorus valence error for " << i << " with " << bcount << std::endl;
#endif
      return false;
    }
    else if (atom_type[i] == Atom_Type::chlorine && bcount != 1) {
#ifdef VERBOSE
      std::cout << "Chlorine valence error for " << i << " with " << bcount << std::endl;
#endif
      return false;
    }
  }
  return true;
}

bool Molecule::consistent() const
{
  int i,j;
  const unsigned int NA = (unsigned) natoms;

  if (coords.size() != 3*NA) return false;
  if (bonds.size() != 4*NA) return false;
  if (btype.size() != 4*NA) return false;
  if (atom_type.size() != NA) return false;
  if (locale.size() != NA) return false;
  for(i=0; i<natoms; ++i) {
    for(j=0; j<4; ++j) {
      if (bonds[4*i+j] == -1) continue;
      if (bonds[4*i+j] < 0 || bonds[4*i+j] >= natoms) {
#ifdef VERBOSE
        std::cout << "Illegal bond value " << i << "  " << j << "  " << bonds[4*i+j] << "  " << natoms << std::endl;
#endif
        return false;
      }
    }
    for(j=0; j<3; ++j) {
      if (std::abs(coords[3*i+j]) > 1000.0) {
#ifdef VERBOSE
        std::cout << "Illegal coordinate value " << i << "  " << j << "  " << coords[3*i+j] << std::endl;
#endif
        return false;
      }
    }
  }
  return valence_check();
}

bool Molecule::create_tbond()
{
  int i,j,l,k,h1,h2,candidate = 0;
  int temp,a,in1,hydro1[4],hydro2[4];
  bool found;
  double L1[3],L2[3];
  std::set<int> drop;

  // First see if we can find any carbon-carbon double bonds
  for(i=0; i<natoms; ++i) {
    if (atom_type[i] == Atom_Type::carbon) {
      found = false;
      for(j=0; j<4; ++j) {
        temp = bonds[4*i+j];
        if (temp >= 0) {
          if (atom_type[j] == Atom_Type::carbon && btype[4*i+j] == 2) {
            // It's a carbon-carbon double bond, so let's see if each carbon has
            // at least one hydrogen
            candidate = j;
            found = true;
          }
        }
      }
      if (!found) continue;
      h1 = 0;
      a = bonds[4*i+candidate];
      for(j=0; j<4; ++j) {
        temp = bonds[4*i+j];
        if (temp >= 0) {
          if (atom_type[temp] == Atom_Type::hydrogen) {
            hydro1[h1] = temp;
            h1++;
          }
        }
      }
      h2 = 0;
      in1 = -1;
      for(j=0; j<4; ++j) {
        temp = bonds[4*a+j];
        if (temp >= 0) {
          if (temp == i) in1 = j;
          if (atom_type[temp] == Atom_Type::hydrogen) {
            hydro2[h2] = temp;
            h2++;
          }
        }
      }
      if (h1 == 0 || h2 == 0) continue;
      // Now we need to see if these hydrogens are cis or trans
      for(j=0; j<3; ++j) {
        L1[j] = coords[3*i+j] - coords[3*candidate+j];
      }
      for(j=0; j<h1; ++j) {
        for(k=0; k<h2; ++k) {
          for(l=0; l<3; ++l) {
            L2[l] = coords[3*hydro1[j]+l] - coords[3*hydro2[k]+l];
          }
          if (parallel(L1,L2)) {
            // This pair of hydrogens is trans, so we can eliminate them and
            // create a triple bond
            btype[4*i+candidate] = 3;
            btype[4*a+in1] = 3;
            drop.insert(hydro1[j]);
            drop.insert(hydro2[k]);
            eliminate_atoms(drop);
            return true;
          }
        }
      }
    }
  }
  return false;
}

bool Molecule::create_dbond()
{
  int i,temp,q,in1,in2;
  unsigned int j,k,l,m,kount;
  bool done,problem,ring1,ring2;
  double p[4][3],A,B,C,D,delta;
  std::set<int> drop;
  std::vector<int> indices,candidate,h1,h2;

  for(i=0; i<natoms; ++i) {
    indices.push_back(i);
  }
  RND.shuffle(indices);
  // First we enumerate all of the carbon-carbon sp3 bonds
  // in the molecule
  for(i=0; i<natoms; ++i) {
    temp = indices[i];
    if (atom_type[temp] == Atom_Type::carbon) {
      for(j=0; j<4; ++j) {
        q = bonds[4*temp+j];
        if (q >= 0) {
          if (atom_type[q] == Atom_Type::carbon && btype[4*temp+j] == 1) {
            if (temp < q) {
              done = false;
              for(k=0; k<candidate.size(); k+=2) {
                if (candidate[k] == temp && candidate[k+1] == q) {
                  done = true;
                  break;
                }
              }
              if (!done) {
                candidate.push_back(temp);
                candidate.push_back(q);
              }
            }
            else {
              done = false;
              for(k=0; k<candidate.size(); k+=2) {
                if (candidate[k] == q && candidate[k+1] == temp) {
                  done = true;
                  break;
                }
              }
              if (!done) {
                candidate.push_back(q);
                candidate.push_back(temp);
              }
            }
          }
        }
      }
    }
  }
  // Now we go through these C-C sp3 bonds, looking for ones where each carbon in the pair
  // has no exotic bonds or neighbours: aromatic rings, other double/triple bonds etc. If we
  // find such a pair, then we assemble a list of each carbon's hydrogen neighbours for
  // further analysis.
  // The simple way to do this is to simply make sure that each carbon in the pair has two
  // hydrogen neighbours, this will eliminate either one from already participating in a double,
  // aromatic or triple bond, and we are not concerned with one of the carbons being bonded to a
  // heteroatom.
  int nc = (signed) candidate.size();
  for(i=0; i<nc; i+=2) {
    problem = false;
    for(q=0; q<nrings; ++q) {
      ring1 = false;
      ring2 = false;
      for(k=0; k<6; ++k) {
        if (candidate[i] == rings[6*q+k]) ring1 = true;
        if (candidate[i+1] == rings[6*q+k]) ring2 = true;
      }
      if (ring1 && ring2) {
        problem = true;
        break;
      }
    }
    if (problem) continue;
    problem = false;
    h1.clear();
    h2.clear();
    for(j=0; j<4; ++j) {
      temp = bonds[4*candidate[i]+j];
      if (temp >= 0) {
        if (atom_type[temp] == Atom_Type::hydrogen) {
          h1.push_back(temp);
        }
        if (btype[4*candidate[i]+j] != 1) problem = true;
      }
    }
    for(j=0; j<4; ++j) {
      temp = bonds[4*candidate[i+1]+j];
      if (temp >= 0) {
        if (atom_type[temp] == Atom_Type::hydrogen) {
          h2.push_back(temp);
        }
        if (btype[4*candidate[i+1]+j] != 1) problem = true;
      }
    }
    if (h1.size() > 0 && h2.size() > 0 && !problem) {
      // Now, two of these hydrogens, along with the two carbons themselves, need to be
      // coplanar in order to form a double bond. The two remaining hydrogens are the ones
      // which will be deleted.
      for(j=0; j<h1.size(); ++j) {
        kount = 0;
        for(l=0; l<4; ++l) {
          temp = bonds[4*candidate[i]+l];
          if (temp == candidate[i+1] || temp == h1[j]) continue;
          for(m=0; m<3; ++m) {
            p[kount][m] = coords[3*temp+m];
          }
          kount++;
        }
#ifdef DEBUG
        assert(kount == 2);
#endif
        for(k=0; k<h2.size(); ++k) {
          for(l=0; l<4; ++l) {
            temp = bonds[4*candidate[i+1]+l];
            if (temp == candidate[i] || temp == h2[k]) continue;
            for(m=0; m<3; ++m) {
              p[kount][m] = coords[3*temp+m];
            }
            kount++;
          }
#ifdef DEBUG
          assert(kount == 4);
#endif
          A = p[0][1]*(p[1][2]-p[2][2])-p[0][2]*(p[1][1]-p[2][1])+(p[1][1]*p[2][2]-p[2][1]*p[1][2]);
          B = -1.0*(p[0][0]*(p[1][2]-p[2][2])-p[0][2]*(p[1][0]-p[2][0])+(p[1][0]*p[2][2]-p[2][0]*p[1][2]));
          C = p[0][0]*(p[1][1]-p[2][1])-p[0][1]*(p[1][0]-p[2][0])+(p[1][0]*p[2][1]-p[2][0]*p[1][1]);
          D = -1.0*(p[0][0]*(p[1][1]*p[2][2]-p[2][1]*p[1][2])-p[0][1]*(p[1][0]*p[2][2]-p[2][0]*p[1][2])+p[0][2]*(p[1][0]*p[2][1]-p[2][0]*p[1][1]));
          delta = (A*p[3][0]+B*p[3][1]+C*p[3][2]+D)/std::sqrt(A*A+B*B+C*C);
          if (std::abs(delta) < 0.001) {
            // We'll delete h1[j] and h2[k] to create the double bond
            for(l=0; l<4; ++l) {
              if (bonds[4*candidate[i+1]+l] == candidate[i]) in1 = l;
              if (bonds[4*candidate[i]+l] == candidate[i+1]) in2 = l;
            }
            btype[4*candidate[i]+in2] = 2;
            btype[4*candidate[i+1]+in1] = 2;
            drop.insert(h1[j]); drop.insert(h2[k]);
            eliminate_atoms(drop);
            return true;
          }
          else {
            kount = 2;
          }
        }
      }
    }
  }
  return false;
}

bool Molecule::decorate(const bool* ornaments)
{
  // The first thing we will want to do here is to determine all of the rings
  // contained in this molecule, and then if specified, to eliminate all of
  // their axial methyl groups.
  int i,j,k,h,temp,ndouble,nops,alpha,hydrogen[3],opcount = 0;
  bool methyl,test;
  std::set<int> drop;
  std::vector<int> axial;

  bool kill_axial = ornaments[0];
  bool create_penta = ornaments[1];
  bool create_double = ornaments[2];
  bool create_triple = ornaments[3];
  bool create_exotic = ornaments[4];
  bool subs_oxy = ornaments[5];
  bool subs_sul = ornaments[6];
  bool subs_nit = ornaments[7];
  bool subs_fun = ornaments[8];

  assert(saturation_check());
  get_rings();
#ifdef VERBOSE
  std::cout << "There are " << nrings << " rings" << std::endl;
#endif
  if (kill_axial) {
    do {
      methyl = false;
      axial_ring_bonds(axial);     
      for(i=0; i<nrings; ++i) {
        for(j=0; j<6; ++j) {
          temp = axial[6*i+j];
          if (atom_type[temp] == Atom_Type::carbon) {
            // See if it's a methyl
            h = 0;
            for(k=0; k<4; ++k) {
              if (bonds[4*temp+k] >= 0) {
                if (atom_type[bonds[4*temp+k]] == Atom_Type::hydrogen) {
                  hydrogen[h] = bonds[4*temp+k];
                  h++;
                }
              }
            }
            if (h == 3) {
              // It's methyl, get rid of it
              methyl = true;
              break;
            }
          }
        }
        if (methyl) break;
      }
      if (methyl) {
        // Actually eliminate the methyl group here
#ifdef VERBOSE
        std::cout << "Eliminating axial methyl..." << std::endl;
#endif
        atom_type[temp] = Atom_Type::hydrogen;
        drop.clear();
        drop.insert(hydrogen[0]); drop.insert(hydrogen[1]); drop.insert(hydrogen[2]);
        eliminate_atoms(drop);
        //eliminate_atoms(hydrogen,3);
        axial.clear();
      }
      else {
        break;
      }
    } while(true);
  }
  else {
    axial_ring_bonds(axial);
  }
  // Next, we want to see which rings can be aromatized, that is all of the
  // axial neighbours of all of the ring atoms are hydrogen atoms.
  // We will convert all aromatizable rings to benzene, and in later operations
  // these benzene rings can be altered by various changes, subject to the condition
  // that the ring's aromaticity be preserved (e.g. we can add a nitrogen to a
  // benzene, but if we want to add an oxygen, then we must eliminate one of the
  // carbon atoms in the ring)
  aromatize(axial);
  // Now, go create some double bonds...
  if (create_double) {
    ndouble = RND.irandom(4);
    for(i=0; i<ndouble; ++i) {
      test = create_dbond();
#ifdef VERBOSE
      if (test) std::cout << "Double bond created" << std::endl;
#endif
    }
  }
 
  // Next step here involves creating an "op string" of a random sequence of the
  // following: T, O, S, N, P, F, A.
  nops = 3 + RND.irandom(6);
  opstring.clear();
  do {
    alpha = RND.irandom(7);
    if (create_triple == false && alpha == 4) continue;
    if (create_exotic == false && alpha == 6) continue;
    if (subs_fun == false && alpha == 5) continue;
    if (subs_oxy == false && alpha == 3) continue;
    if (subs_sul == false && alpha == 2) continue;
    if (subs_nit == false && alpha == 1) continue;
    if (create_penta == false && alpha == 0) continue;
    opstring += ops[alpha];
    opcount++;
  } while(opcount < nops);

  for(i=0; i<nops; ++i) {
    if (opstring[i] == 'T') {
      test = create_tbond();
    }
    else if (opstring[i] == 'O') {
      test = add_oxygen();
    }
    else if (opstring[i] == 'S') {
      test = add_sulfur();
    }
    else if (opstring[i] == 'N') {
      test = add_nitrogen();
    }
    else if (opstring[i] == 'P') {
      test = create_penta1();
      if (test == 1) break;
    }
    else if (opstring[i] == 'F') {
      test = fungrp();
    }
    else if (opstring[i] == 'A') {
      test = create_amide();
    }
#ifdef DEBUG
    assert(consistent());
#endif
  }
  test = normalize_aromatic_bonds();
  if (!test) {
#ifdef VERBOSE
    std::cout << "Problem in normalizing aromatic ring" << std::endl;
#endif
    return false;
  }
  return true;
}


