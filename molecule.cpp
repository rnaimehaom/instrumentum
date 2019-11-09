#include "molecule.h"

extern Random RND;

extern std::map<int,std::string> element_table;

const char Molecule::ops[7] = {'P','N','S','O','T','F','A'};

Molecule::Molecule()
{

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

void Molecule::add_atom(int atomic_number,const double* x,int ltype)
{
#ifdef DEBUG
  assert(atomic_number > 0);
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
  std::vector<int> nbonds,nbtype,nlocal,ntype;
  std::vector<double> ncoords;

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

int Molecule::write(std::ofstream& s) const
{
  // This method writes the Molecule properties "natoms", "opstring", "bonds", 
  // "btype" and "coords" to a binary diskfile.
  int i,j,n,bcount = 0;
  double x[3];

  s.write((char*)(&natoms),sizeof(int)); bcount += sizeof(int);
  n = opstring.length();
  s.write((char*)(&n),sizeof(int)); bcount += sizeof(int);
  for(i=0; i<n; ++i) {
    s.write((char*)(&opstring[i]),sizeof(char));   
  }
  bcount += n*sizeof(char);

  for(i=0; i<natoms; ++i) {
    n = atom_type[i];
    s.write((char*)(&n),sizeof(int));
  }
  bcount += natoms*sizeof(int);

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

  clear();

  s.read((char*)(&natoms),sizeof(int)); bcount += sizeof(int);
  s.read((char*)(&n),sizeof(int)); bcount += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&c),sizeof(char)); 
    opstring += c; 
  }
  bcount += n*sizeof(char);

  for(i=0; i<natoms; ++i) {
    s.read((char*)(&n),sizeof(int));
    atom_type.push_back(n);
  }
  bcount += natoms*sizeof(int);

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
    s << i << "  " << source.atom_type[i] << "  " << source.locale[i] << std::endl;
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

bool Molecule::saturation_check() const
{
  int i,j,k,nz,bcount;
  bool found;

  for(i=0; i<natoms; ++i) {
    if (atom_type[i] == 1) {
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
    else if (atom_type[i] == 6) {
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
    if (atom_type[i] == 1 && bcount != 1) {
#ifdef VERBOSE
      std::cout << "Hydrogen valence error for " << i << " with " << bcount << std::endl;
#endif
      return false;
    }
    else if (atom_type[i] == 9 && bcount != 1) {
#ifdef VERBOSE
      std::cout << "Fluorine valence error for " << i << " with " << bcount << std::endl;
#endif
      return false;
    }
    else if (atom_type[i] == 15 && bcount != 3) {
#ifdef VERBOSE
      std::cout << "Phosphorus valence error for " << i << " with " << bcount << std::endl;
#endif
      return false;
    }
    else if (atom_type[i] == 17 && bcount != 1) {
#ifdef VERBOSE
      std::cout << "Chlorine valence error for " << i << " with " << bcount << std::endl;
#endif
      return false;
    }
  }
  return true;
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

bool Molecule::add_nitrogen()
{
  // This routine searches for a carbon atom with at least two hydrogen atom
  // bonds, and which is not beside any other heteroatoms or involved in any
  // unsaturated bonds. It will then seek to convert it to an oxygen atom.
  int i,j,k,l,h,c,con,hh,hcount,alpha,nb,naroma,temp1,temp2,tt,in1;
  int carb1,carb2,the_ring,cand,c1,c2,temp,in2;
  int candidate[3],aroma_neighbours[3],hydrogen[3],carbon[3];
  int candidate_rneighbours[3][2],candidate_hneighbours[3][1];
  bool problem,exotic;
  std::set<int> drop;
  std::vector<int> indices;

  for(i=0; i<natoms; ++i) {
    indices.push_back(i);
  }
  RND.shuffle(indices);
  for(i=0; i<natoms; ++i) {
    temp1 = indices[i];
    if (atom_type[temp1] == 6) {
      // Now let's see how many hydrogen atoms this carbon is bonded to
      if (in_ring(temp1)) {
        // See if the ring(s) already contain two heteroatoms
        problem = false;
        for(j=0; j<nrings; ++j) {
          con = 0;
          hcount = 0;
          for(k=0; k<6; ++k) {
            tt = rings[6*j+k];
            if (tt == temp1) con = 1;
            if (tt >= 0) {
              if (atom_type[tt] != 6) hcount++;
            }
          }
          if (con == 1 && hcount >= 2) {
            problem = true;
            break;
          }
        }
        if (problem) continue;
        // Next let's see if the ring is aromatic
        if (in_aromatic(temp1)) {
          // We have here either carbon with two carbon neighbours and a hydrogen, so
          // kill the H and make the C into N, or carbon with *three* carbon neighbours.
          // In the latter case, check to make sure that at least one of these carbons
          // is in a chain or non-aromatic ring, then make the C a N and knock out one
          // of the aromatic ring carbons
          h = 0;
          c = 0;
          exotic = false;
          nb = 0;
          for(j=0; j<4; ++j) {
            temp2 = bonds[4*temp1+j];
            if (temp2 >= 0) {
              nb++;
              if (atom_type[temp2] == 1) {
                hydrogen[h] = temp2;
                h++;
              }
              else if (atom_type[temp2] != 6) {
                exotic = true;
              }
              else if (atom_type[temp2] == 6) {
                carbon[c] = temp2;
                c++;
              }
            }
          }
#ifdef DEBUG
          assert(nb == 3);
#endif
          if (exotic) continue;
          if (h == 0) {
#ifdef DEBUG
            assert(c == 3);
#endif
            naroma = 0;
            for(j=0; j<c; ++j) {
              if (in_aromatic(carbon[j])) {
                aroma_neighbours[naroma] = carbon[j];
                naroma++;
              }
            }
#ifdef DEBUG
            assert(naroma >= 2);
#endif
            if (naroma == 3) continue;
            // Also check that this isn't already a five-membered ring!
            the_ring = -1;
            carb1 = aroma_neighbours[0];
            carb2 = aroma_neighbours[1];
            for(j=0; j<nrings; ++j) {
              if (get_rindex(temp1,j) != -1 && get_rindex(carb1,j) != -1 && get_rindex(carb2,j) != -1) {
                the_ring = j;
                break;
              }
            }
#ifdef DEBUG
            assert(the_ring >= 0);
#endif
            if (rings[6*the_ring+5] == -1) continue;
            // So now we need to find another carbon in this ring that can be
            // sacrificed
            cand = 0;
            for(j=0; j<6; ++j) {
              if (rings[6*the_ring+j] != temp1 && rings[6*the_ring+j] != carb1 && rings[6*the_ring+j] != carb2) {
                if (atom_type[rings[6*the_ring+j]] != 6) continue;
                hh = -1;
                c1 = -1;
                c2 = -1;
                for(k=0; k<4; ++k) {
                  temp = bonds[4*rings[6*the_ring+j]+k];
                  if (temp < 0) continue;
                  if (atom_type[temp] == 1) {
                    hh = temp;
                  }
                  else if (get_rindex(temp,the_ring) != -1) {
                    if (c1 == -1) {
                      c1 = temp;
                    }
                    else if (c1 >= 0) {
                      c2 = temp;
                    }
                  }
                }
                if (hh >= 0) {
                  candidate[cand] = rings[6*the_ring+j];
                  candidate_hneighbours[cand][0] = hh;
                  candidate_rneighbours[cand][0] = c1;
                  candidate_rneighbours[cand][1] = c2;
                  cand++;
                }
              }
            }
            if (cand < 1) continue;
            alpha = RND.irandom(cand);
            atom_type[temp1] = 7;
            drop.clear();
            drop.insert(candidate[alpha]);
            drop.insert(candidate_hneighbours[alpha][0]);
            // Get the two ring neighbours of candidate[alpha] to bond to each other rather than
            // candidate[alpha]
            for(l=0; l<4; ++l) {
              if (candidate[alpha] == bonds[4*candidate_rneighbours[alpha][0]+l]) in1 = l;
              if (candidate[alpha] == bonds[4*candidate_rneighbours[alpha][1]+l]) in2 = l;
            }
            bonds[4*candidate_rneighbours[alpha][0]+in1] = candidate_rneighbours[alpha][1];
            bonds[4*candidate_rneighbours[alpha][1]+in2] = candidate_rneighbours[alpha][0];
            eliminate_atoms(drop);
#ifdef VERBOSE
            std::cout << this << std::endl;
            for(l=0; l<nrings; ++l) {
              for(int m=0; m<6; ++m) {
                std::cout << rings[6*l+m] << "  ";
              }
              std::cout << std::endl;
            }
#endif
          }
          else if (h == 1) {
            // Let's make the conversion...
            atom_type[temp1] = 7;
            drop.clear();
            drop.insert(hydrogen[0]);
            eliminate_atoms(drop);
            return true;
          }
        }
      }
      h = 0;
      exotic = false;
      nb = 0;
      for(j=0; j<4; ++j) {
        temp2 = bonds[4*temp1+j];
        if (temp2 > 0) {
          nb++;
          if (atom_type[temp2] == 1) {
            hydrogen[h] = temp2;
            h++;
          }
          else if (atom_type[temp2] != 6) {
            exotic = true;
          }
        }
        else {
          break;
        }
      }
      if (h >= 1 && !exotic && nb >= 3) {
        // Let's make the conversion...
        atom_type[temp1] = 7;
        drop.clear();
        drop.insert(hydrogen[0]);
        eliminate_atoms(drop);
        return true;
      }
    }
  }
  return false;
}

bool Molecule::add_sulfur()
{
  // This routine searches for a carbon atom with at least two hydrogen atom
  // bonds, and which is not beside any other heteroatoms or involved in any
  // unsaturated bonds. It will then seek to convert it to a sulfur atom.
  int i,j,k,h,con,hcount,nb;
  int hydrogen[3],temp1,temp2,tt;
  bool exotic,problem;
  std::set<int> drop;
  std::vector<int> indices;

  for(i=0; i<natoms; ++i) {
    indices.push_back(i);
  }
  RND.shuffle(indices);
  for(i=0; i<natoms; ++i) {
    temp1 = indices[i];
    if (atom_type[temp1] == 6) {
      // Now let's see how many hydrogen atoms this carbon is bonded to
      if (in_ring(temp1)) {
        // Check to make sure there aren't too many heteroatoms in these ring(s)
        problem = false;
        for(j=0; j<nrings; ++j) {
          con = 0;
          hcount = 0;
          for(k=0; k<6; ++k) {
            tt = rings[6*j+k];
            if (tt == temp1) con = 1;
            if (tt >= 0) {
              if (atom_type[tt] != 6) hcount++;
            }
          }
          if (con == 1 && hcount >= 2) {
            problem = true;
            break;
          }
        }
        if (problem) continue;
      }
      h = 0;
      exotic = false;
      nb = 0;
      for(j=0; j<4; ++j) {
        temp2 = bonds[4*temp1+j];
        if (temp2 > 0) {
          nb++;
          if (atom_type[temp2] == 1) {
            hydrogen[h] = temp2;
            h++;
          }
          else if (atom_type[temp2] != 6) {
            exotic = true;
          }
        }
        else {
          break;
        }
      }
      if (h >= 2 && !exotic && nb == 4) {
        // Let's make the conversion...
        atom_type[temp1] = 16;
        drop.insert(hydrogen[0]); drop.insert(hydrogen[1]);
        eliminate_atoms(drop);
        return true;
      }
    }
  }
  return false;
}

bool Molecule::add_oxygen()
{
  // This routine searches for a carbon atom with at least two hydrogen atom
  // bonds, and which is not beside any other heteroatoms or involved in any
  // unsaturated bonds. It will then seek to convert it to an oxygen atom.
  int i,j,k,h,nb,hcount,con;
  int hydrogen[3],temp1,temp2,tt;
  bool problem,exotic;
  std::set<int> drop;
  std::vector<int> indices;

  for(i=0; i<natoms; ++i) {
    indices.push_back(i);
  }
  RND.shuffle(indices);
  for(i=0; i<natoms; ++i) {
    temp1 = indices[i];
    if (atom_type[temp1] == 6) {
      // Now let's see how many hydrogen atoms this carbon is bonded to
      if (in_ring(temp1)) {
        // Check to make sure there aren't too many heteroatoms in these ring(s)
        problem = false;
        for(j=0; j<nrings; ++j) {
          con = 0;
          hcount = 0;
          for(k=0; k<6; ++k) {
            tt = rings[6*j+k];
            if (tt == temp1) con = 1;
            if (tt >= 0) {
              if (atom_type[tt] != 6) hcount++;
            }
          }
          if (con == 1 && hcount >= 2) {
            problem = true;
            break;
          }
        }
        if (problem) continue;
      }

      h = 0;
      exotic = false;
      nb = 0;
      for(j=0; j<4; ++j) {
        temp2 = bonds[4*temp1+j];
        if (temp2 > 0) {
          nb++;
          if (atom_type[temp2] == 1) {
            hydrogen[h] = temp2;
            h++;
          }
          else if (atom_type[temp2] != 6) {
            exotic = true;
          }
        }
        else {
          break;
        }
      }
      if (h >= 2 && !exotic && nb == 4) {
        // Let's make the conversion...
        atom_type[temp1] = 8;
        drop.insert(hydrogen[0]); drop.insert(hydrogen[1]);
        eliminate_atoms(drop);
        return true;
      }
    }
  }
  return false;
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
    if (atom_type[i] == 6) {
      found = false;
      for(j=0; j<4; ++j) {
        temp = bonds[4*i+j];
        if (temp >= 0) {
          if (atom_type[j] == 6 && btype[4*i+j] == 2) {
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
          if (atom_type[temp] == 1) {
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
          if (atom_type[temp] == 1) {
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

int Molecule::aromatize(std::vector<int>& axial)
{
  int aromatize,in1,temp,rindex,naromatized = 0;
  int i,j,k,l,c,change[4],kount = 0;
  bool good;
  std::set<int> drop;

  do {
    aromatize = -1;
    for(i=kount; i<nrings; ++i) {
      good = true;
      for(j=0; j<6; ++j) {
        if (axial[6*i+j] == -1) {
#ifdef VERBOSE
          std::cout << i << "  Empty" << std::endl;
#endif
          continue;
        }
        else {
#ifdef VERBOSE
          std::cout << i << "  " << atom_type[axial[6*i+j]] << std::endl;
#endif
        }
        if (atom_type[axial[6*i+j]] != 1) {
          good = false;
          break;
        }
      }
      if (good) {
        aromatize = i;
        break;
      }
    }
    if (aromatize == -1) break;
    kount = aromatize + 1;
    // We need to change all the ring bond types to four, and then
    // get rid of all these axial hydrogens
#ifdef VERBOSE
    std::cout << "Aromatizing ring number " << aromatize << std::endl;
#endif
    for(i=0; i<6; ++i) {
      rindex = rings[6*aromatize+i];
      c = 0;
      for(j=0; j<4; ++j) {
        temp = bonds[4*rindex+j];
        if (temp >= 0) {
          for(k=0; k<6; ++k) {
            if (rings[6*aromatize+k] == temp) {
              change[c] = j;
              c++;
              break;
            }
          }
        }
      }
#ifdef DEBUG
      assert(c == 2);
#endif
      for(j=0; j<c; ++j) {
        btype[4*rindex+change[j]] = 4;
        in1 = bonds[4*rindex+change[j]];
        for(k=0; k<4; ++k) {
          if (bonds[4*in1+k] == rindex) {
            l = k;
            break;
          }
        }
        btype[4*in1+l] = 4;
      }
    }
    drop.clear();
    for(l=0; l<6; ++l) {
      if (axial[6*aromatize+l] != -1) drop.insert(axial[6*aromatize+l]);
    }
    assert(eliminate_atoms(drop,axial));
    naromatized++;
  } while(true);

  return naromatized;
}

int Molecule::get_rings()
{
  int i,j,k,in1,in2;
  unsigned int l;
  bool found;
  std::vector<int> ring_vertices,ring_edges,redges,ring_data,wcopy;

  wcopy = bonds;

  for(i=0; i<natoms; ++i) {
    for(j=0; j<4; ++j) {
      if (wcopy[4*i+j] == -1) break;
      k = wcopy[4*i+j];
      wcopy[4*i+j] = -1;
      for(l=0; l<4; ++l) {
        if (wcopy[4*k+l] == i) {
          wcopy[4*k+l] = -1;
          break;
        }
      }
      if (connected(wcopy)) {
        found = false;
        for(l=0; l<ring_vertices.size(); ++l) {
          if (ring_vertices[l] == i) {
            found = true;
            break;
          }
        }
        if (!found) ring_vertices.push_back(i);
        found = false;
        for(l=0; l<ring_vertices.size(); ++l) {
          if (ring_vertices[l] == k) {
            found = true;
            break;
          }
        }
        if (!found) ring_vertices.push_back(k);
        found = false;
        for(l=0; l<ring_edges.size(); l+=2) {
          in1 = ring_edges[l];
          in2 = ring_edges[l+1];
          if ((in1 == i && in2 == k) || (in1 == k && in2 == i)) {
            found = true;
            break;
          }
        }
        if (!found) {
          ring_edges.push_back(i);
          ring_edges.push_back(k);
        }
      }
      wcopy = bonds;
    }
  }

  if (ring_vertices.empty()) {
    nrings = 0;
    return 0;
  }

  for(l=0; l<4*ring_vertices.size(); ++l) {
    redges.push_back(-1);
  }
  for(l=0; l<ring_edges.size(); l+=2) {
    in1 = get_index(ring_edges[l],ring_vertices);
    in2 = get_index(ring_edges[l+1],ring_vertices);
    for(j=0; j<4; ++j) {
      if (redges[4*in1+j] >= 0) continue;
      redges[4*in1+j] = in2;
      break;
    }
    for(j=0; j<4; ++j) {
      if (redges[4*in2+j] >= 0) continue;
      redges[4*in2+j] = in1;
      break;
    }
  }
  ring_perception(redges,ring_data);
  rings.clear();
  for(l=0; l<ring_data.size(); ++l) {
    rings.push_back(ring_vertices[ring_data[l]]);
  }
  nrings = rings.size()/6;
  return nrings;
}

bool Molecule::fungrp()
{
  // For now, this routine searches for terminal hydrogens that can be converted
  // to a halogen, essentially chlorine or fluorine.
  int i,alpha,temp1,temp2;
  std::vector<int> indices;
  double p1[3],p2[3],xc[3];

  for(i=0; i<natoms; ++i) {
    indices.push_back(i);
  }
  RND.shuffle(indices);
  for(i=0; i<natoms; ++i) {
    temp1 = indices[i];
    if (atom_type[temp1] == 1) {
      // Now see if this hydrogen is bonded to a carbon
      temp2 = bonds[4*temp1];
      if (atom_type[temp2] == 6) {
        // Do the conversion...
        alpha = RND.irandom(3);
        if (alpha == 0) {
          atom_type[temp1] = 9;
        }
        else if (alpha == 1) {
          atom_type[temp1] = 17;
        }
        else {
          // Phosphorus - so we need to add two hydrogens that will bond to this
          // P atom, since it has a valence of three, like nitrogen.
          // Something broken here with the bonds I think!!!
          atom_type[temp1] = 15;

          p1[0] = coords[3*temp2];
          p1[1] = coords[3*temp2+1];
          p1[2] = coords[3*temp2+2];
          p2[0] = coords[3*temp1] - coords[3*temp2];
          p2[1] = coords[3*temp1+1] - coords[3*temp2+1];
          p2[2] = coords[3*temp1+2] - coords[3*temp2+2];

          xc[0] = p1[0] + 2.2*p2[0];
          xc[1] = p1[1] + 2.2*p2[1] + 0.3;
          xc[2] = p1[2] + 2.2*p2[2];
          add_atom(1,xc,0);
          add_bond(temp1,natoms-1,1);

          xc[0] = p1[0] + 2.2*p2[0];
          xc[1] = p1[1] + 2.2*p2[1] - 0.3;
          xc[2] = p1[2] + 2.2*p2[2];
          add_atom(1,xc,0);
          add_bond(temp1,natoms-1,1);
        }
        return true;
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
    if (atom_type[temp] == 6) {
      for(j=0; j<4; ++j) {
        q = bonds[4*temp+j];
        if (q >= 0) {
          if (atom_type[q] == 6 && btype[4*temp+j] == 1) {
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
        if (atom_type[temp] == 1) {
          h1.push_back(temp);
        }
        if (btype[4*candidate[i]+j] != 1) problem = true;
      }
    }
    for(j=0; j<4; ++j) {
      temp = bonds[4*candidate[i+1]+j];
      if (temp >= 0) {
        if (atom_type[temp] == 1) {
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

bool Molecule::is_aromatic(int rnumber) const
{
  // Returns true if this ring is aromatic, false otherwise
  int i,j,k;
  int temp,atom;
  bool found;
  for(i=0; i<6; ++i) {
    atom = rings[6*rnumber+i];
    if (atom < 0) continue;
    for(j=0; j<4; ++j) {
      temp = bonds[4*atom+j];
      if (temp < 0) continue;
      found = false;
      for(k=0; k<6; ++k) {
        if (temp == rings[6*rnumber+k]) {
          found = true;
          break;
        }
      }
      if (found) {
        if (btype[4*atom+j] != 4) return false;
      }
    }
  }
  return true;
}

bool Molecule::create_penta1()
{
  int i,j,k,cand,nitro,n_atoms[3],bcount1,bcount2,n,alpha,alpha1,alpha2,h;
  int candidates[6],candidate_hneighbours[6][2],candidate_rneighbours[6][2];
  int c1,c2,temp,in1,in2,the_nitro,hydro[2];
  double p2[3],rcentre[3],xc[3];
  std::set<int> drop;

  for(i=0; i<nrings; ++i) {
    if (rings[6*i+5] == -1) continue;
    if (is_aromatic(i)) {
      // Check to see if this ring contains zero, one or two nitrogen atoms
      nitro = 0;
      for(j=0; j<6; ++j) {
        if (atom_type[rings[6*i+j]] == 7) nitro++;
      }
      if (nitro == 0) {
        // Benzene: convert to furan or thiophene
        cand = 0;
        // So, find two disposable CH groups, one to delete and the other to convert
        // to oxygen
        for(j=0; j<6; ++j) {
          h = 0;
          c1 = -1;
          c2 = -1;
          for(k=0; k<4; ++k) {
            temp = bonds[4*rings[6*i+j]+k];
            if (temp >= 0) {
              if (atom_type[temp] == 1) {
                h = temp;
              }
              else if (get_rindex(temp,i) != -1) {
                if (c1 == -1) {
                  c1 = temp;
                }
                else if (c1 >= 0) {
                  c2 = temp;
                }
              }
            }
          }
          if (h > 0) {
            candidates[cand] = rings[6*i+j];
            candidate_hneighbours[cand][0] = h;
            candidate_rneighbours[cand][0] = c1;
            candidate_rneighbours[cand][1] = c2;
            cand++;
          }
        }
        if (cand < 2) continue;
        // Choose two at random, and do the necessary changes
        alpha1 = RND.irandom(cand);
        do {
          alpha2 = RND.irandom(cand);
          if (alpha2 != alpha1) break;
        } while(true);
        drop.insert(candidate_hneighbours[alpha1][0]);
        drop.insert(candidates[alpha2]);
        drop.insert(candidate_hneighbours[alpha2][0]);
        alpha = RND.irandom(2);
        if (alpha == 0) {
          atom_type[candidates[alpha1]] = 8;
        }
        else {
          atom_type[candidates[alpha1]] = 16;
        }
        in1 = get_bindex(candidates[alpha2],candidate_rneighbours[alpha2][0]);
        in2 = get_bindex(candidates[alpha2],candidate_rneighbours[alpha2][1]);
        bonds[4*candidate_rneighbours[alpha2][0]+in1] = candidate_rneighbours[alpha2][1];
        bonds[4*candidate_rneighbours[alpha2][1]+in2] = candidate_rneighbours[alpha2][0];
        eliminate_atoms(drop);
        return true;
      }
      else if (nitro == 1) {
        // Delete a carbon, and add a hydrogen to this nitrogen
        cand = 0;
        the_nitro = -1;
        for(j=0; j<6; ++j) {
          if (atom_type[rings[6*i+j]] == 7) the_nitro = rings[6*i+j];
          h = 0;
          c1 = -1;
          c2 = -1;
          for(k=0; k<4; ++k) {
            temp = bonds[4*rings[6*i+j]+k];
            if (temp >= 0) {
              if (atom_type[temp] == 1) {
                h = temp;
              }
              else if (get_rindex(temp,i) != -1) {
                if (c1 == -1) {
                  c1 = temp;
                }
                else if (c1 >= 0) {
                  c2 = temp;
                }
              }
            }
          }
          if (h > 0) {
            candidates[cand] = rings[6*i+j];
            candidate_hneighbours[cand][0] = h;
            candidate_rneighbours[cand][0] = c1;
            candidate_rneighbours[cand][1] = c2;
            cand++;
          }
        }
        bcount1 = 0;
        for(k=0; k<4; ++k) {
          if (bonds[4*the_nitro+k] >= 0) {
            bcount1++;
          }
        }
#ifdef DEBUG
        assert(bcount1 == 2);
#endif
        if (cand < 1) continue;
        // Choose two at random, and do the necessary changes
        alpha1 = RND.irandom(cand);
        drop.insert(candidate_hneighbours[alpha1][0]);
        drop.insert(candidates[alpha1]);
        // Let's compute the ring centre
        rcentre[0] = 0.0;
        rcentre[1] = 0.0;
        rcentre[2] = 0.0;
        for(j=0; j<6; ++j) {
          if (rings[6*i+j] >= 0) {
            for(k=0; k<3; ++k) {
              rcentre[k] += coords[3*rings[6*i+j]+k];
            }
          }
        }
        p2[0] = coords[3*the_nitro] - rcentre[0];
        p2[1] = coords[3*the_nitro+1] - rcentre[1];
        p2[2] = coords[3*the_nitro+2] - rcentre[2];
        xc[0] = rcentre[0] + 2.2*p2[0];
        xc[1] = rcentre[1] + 2.2*p2[1];
        xc[2] = rcentre[2] + 2.2*p2[2];
        // Now add a hydrogen to this lone nitrogen atom
        add_atom(1,xc,0);
        add_bond(the_nitro,natoms-1,1);
        in1 = get_bindex(candidates[alpha1],candidate_rneighbours[alpha1][0]);
        in2 = get_bindex(candidates[alpha1],candidate_rneighbours[alpha1][1]);
        bonds[4*candidate_rneighbours[alpha1][0]+in1] = candidate_rneighbours[alpha1][1];
        bonds[4*candidate_rneighbours[alpha1][1]+in2] = candidate_rneighbours[alpha1][0];
        eliminate_atoms(drop);
        return true;
      }
      else {
#ifdef DEBUG
        assert(nitro <= 2);
#endif        
        // Delete a carbon and add a hydrogen to one of the nitrogens
        cand = 0;
        n = 0;
        for(j=0; j<6; ++j) {
          if (atom_type[rings[6*i+j]] == 7) {
            n_atoms[n] = rings[6*i+j];
            n++;
          }
          h = 0;
          c1 = -1;
          c2 = -1;
          for(k=0; k<4; ++k) {
            temp = bonds[4*rings[6*i+j]+k];
            if (temp >= 0) {
              if (atom_type[temp] == 1) {
                h = temp;
              }
              else if (get_rindex(temp,i) != -1) {
                if (c1 == -1) {
                  c1 = temp;
                }
                else if (c1 >= 0) {
                  c2 = temp;
                }
              }
            }
          }
          if (h > 0) {
            candidates[cand] = rings[6*i+j];
            candidate_hneighbours[cand][0] = h;
            candidate_rneighbours[cand][0] = c1;
            candidate_rneighbours[cand][1] = c2;
            cand++;
          }
        }
        bcount1 = 0;
        bcount2 = 0;
        for(k=0; k<4; ++k) {
          if (bonds[4*n_atoms[0]+k] >= 0) {
            bcount1++;
          }
        }
#ifdef DEBUG
        assert(bcount1 == 2);
#endif
        for(k=0; k<4; ++k) {
          if (bonds[4*n_atoms[1]+k] >= 0) {
            bcount2++;
          }
        }
#ifdef DEBUG
        assert(bcount2 == 2);
#endif
        if (cand < 1) continue;
        // Choose two at random, and do the necessary changes
        alpha1 = RND.irandom(cand);
        alpha2 = RND.irandom(2);
        drop.insert(candidate_hneighbours[alpha1][0]);
        drop.insert(candidates[alpha1]);

        rcentre[0] = 0.0;
        rcentre[1] = 0.0;
        rcentre[2] = 0.0;
        for(j=0; j<6; ++j) {
          if (rings[6*i+j] >= 0) {
            for(k=0; k<3; ++k) {
              rcentre[k] += coords[3*rings[6*i+j]+k];
            }
          }
        }
        p2[0] = coords[3*n_atoms[alpha2]] - rcentre[0];
        p2[1] = coords[3*n_atoms[alpha2]+1] - rcentre[1];
        p2[2] = coords[3*n_atoms[alpha2]+2] - rcentre[2];
        xc[0] = rcentre[0] + 2.2*p2[0];
        xc[1] = rcentre[1] + 2.2*p2[1];
        xc[2] = rcentre[2] + 2.2*p2[2];
        // Now add a hydrogen to this lone nitrogen atom
        add_atom(1,xc,0);
        add_bond(n_atoms[alpha2],natoms-1,1);

        in1 = get_bindex(candidates[alpha1],candidate_rneighbours[alpha1][0]);
        in2 = get_bindex(candidates[alpha1],candidate_rneighbours[alpha1][1]);
        bonds[4*candidate_rneighbours[alpha1][0]+in1] = candidate_rneighbours[alpha1][1];
        bonds[4*candidate_rneighbours[alpha1][1]+in2] = candidate_rneighbours[alpha1][0];
        eliminate_atoms(drop);
        return true;
      }
    }
    else {
      // Just find a carbon with two hydrogen neighbours and kill it
      cand = 0;
      for(j=0; j<6; ++j) {
        if (atom_type[rings[6*i+j]] != 6) continue;
        h = 0;
        c1 = 0;
        c2 = 0;
        for(k=0; k<4; ++k) {
          temp = bonds[4*rings[6*i+j]+k];
          if (temp >= 0) {
            if (atom_type[temp] == 1) {
              hydro[h] = temp;
              h++;
            }
            else if (get_rindex(temp,i) != -1) {
              if (c1 == 0) {
                c1 = temp;
              }
              else {
                c2 = temp;
              }
            }
          }
        }
        if (h == 2) {
          candidates[cand] = rings[6*i+j];
          candidate_hneighbours[cand][0] = hydro[0];
          candidate_hneighbours[cand][1] = hydro[1];
          candidate_rneighbours[cand][0] = c1;
          candidate_rneighbours[cand][1] = c2;
          cand++;
        }
      }
      if (cand < 1) continue;
      // Choose two at random, and do the necessary changes
      alpha1 = RND.irandom(cand);
      drop.insert(candidate_hneighbours[alpha1][0]);
      drop.insert(candidates[alpha1]);
      drop.insert(candidate_hneighbours[alpha1][1]);
      in1 = get_bindex(candidates[alpha1],candidate_rneighbours[alpha1][0]);
      in2 = get_bindex(candidates[alpha1],candidate_rneighbours[alpha1][1]);
      bonds[4*candidate_rneighbours[alpha1][0]+in1] = candidate_rneighbours[alpha1][1];
      bonds[4*candidate_rneighbours[alpha1][1]+in2] = candidate_rneighbours[alpha1][0];
      eliminate_atoms(drop);
      return true;
    }
  }
  return false;
}

bool Molecule::create_amide()
{
  // First try to find a C-C double bond such that each carbon has at least one other
  // carbon neighbour, and one of them also has a hydrogen neighbour, something like
  //  H - C = C - C
  //      |   |
  //      C   C
  // We will transform this into one of
  // C - N - C = O
  //     |   |           AMIDE
  //     C   C
  // or
  // C - N - C = S
  //     |   |           SULFONAMIDE
  //     C   C
  // or finally
  // C - O - C = O
  //         |           ESTER
  //         C
  int i,j,a,in1,h1,c1,h2,c2,alpha,candidate,hydro1[4],hydro2[4];
  std::set<int> drop;

  // First see if we can find any carbon-carbon double bonds
  for(i=0; i<natoms; ++i) {
    if (atom_type[i] == 6) {
      candidate = -1;
      for(j=0; j<4; ++j) {
        if (bonds[4*i+j] < 0) continue;
        if (atom_type[j] == 6 && btype[4*i+j] == 2) {
          // It's a carbon-carbon double bond, so let's see if each carbon has
          // at least one hydrogen
          candidate = j;
        }
      }
      if (candidate == -1) continue;
      h1 = 0;
      c1 = 0;
      a = bonds[4*i+candidate];
      for(j=0; j<4; ++j) {
        if (bonds[4*i+j] < 0 || bonds[4*i+j] == a) continue;
        if (atom_type[bonds[4*i+j]] == 1) {
          hydro1[h1] = bonds[4*i+j];
          h1++;
        }
        else if (atom_type[bonds[4*i+j]] == 6) {
          c1++;
        }
      }
      h2 = 0;
      c2 = 0;
      for(j=0; j<4; ++j) {
        if (bonds[4*a+j] < 0 || bonds[4*a+j] == i) continue;
        if (atom_type[bonds[4*a+j]] == 1) {
          hydro2[h2] = bonds[4*a+j];
          h2++;
        }
        else if (atom_type[bonds[4*a+j]] == 6) {
          c2++;
        }
      }
      if (c1 > 0 && c2 > 0) {
        // For an amide, the hydrogen becomes an oxygen, the other carbon becomes a nitrogen, 
        // and what was a carbon-carbon double bond becomes a nitrogen-carbon single bond.
        // So, first find out who has a hydrogen
        if ((h1 + h2) == 1) {
          alpha = RND.irandom(2);
          if (alpha == 0) {
            if (h1 > 0) {
              atom_type[a] = 7;
              atom_type[hydro1[0]] = 8;
              in1 = get_bindex(hydro1[0],i);
              btype[4*i+in1] = 2;
              btype[4*hydro1[0]] = 2;
              in1 = get_bindex(a,i);
              btype[4*i+in1] = 1;
              in1 = get_bindex(i,a);
              btype[4*a+in1] = 1;
            }
            else {
              atom_type[i] = 7;
              atom_type[hydro2[0]] = 8;
              in1 = get_bindex(hydro2[0],a);
              btype[4*a+in1] = 2;
              btype[4*hydro2[0]+0] = 2;
              in1 = get_bindex(a,i);
              btype[4*i+in1] = 1;
              in1 = get_bindex(i,a);
              btype[4*a+in1] = 1;
            }
            return true;
          }
          else {
            // Sulfonamide
            if (h1 > 0) {
              atom_type[a] = 7;
              atom_type[hydro1[0]] = 16;
              in1 = get_bindex(hydro1[0],i);
              btype[4*i+in1] = 2;
              btype[4*hydro1[0]] = 2;
              in1 = get_bindex(a,i);
              btype[4*i+in1] = 1;
              in1 = get_bindex(i,a);
              btype[4*a+in1] = 1;
            }
            else {
              atom_type[i] = 7;
              atom_type[hydro2[0]] = 16;
              in1 = get_bindex(hydro2[0],a);
              btype[4*a+in1] = 2;
              btype[4*hydro2[0]] = 2;
              in1 = get_bindex(a,i);
              btype[4*i+in1] = 1;
              in1 = get_bindex(i,a);
              btype[4*a+in1] = 1;
            }
            return true;
          }
        }
        else if ((h1 + h2) == 2) {
          alpha = RND.irandom(3);
          if (alpha == 0) {
            atom_type[a] = 7;
            atom_type[hydro1[0]] = 16;
            in1 = get_bindex(hydro1[0],i);
            btype[4*i+in1] = 2;
            btype[4*hydro1[0]] = 2;
            in1 = get_bindex(a,i);
            btype[4*i+in1] = 1;
            in1 = get_bindex(i,a);
            btype[4*a+in1] = 1;
          }
          else if (alpha == 1) {
            atom_type[a] = 7;
            atom_type[hydro1[0]] = 8;
            in1 = get_bindex(hydro1[0],i);
            btype[4*i+in1] = 2;
            btype[4*hydro1[0]] = 2;
            in1 = get_bindex(a,i);
            btype[4*i+in1] = 1;
            in1 = get_bindex(i,a);
            btype[4*a+in1] = 1;
          }
          else {
            // An ester
#ifdef DEBUG
            assert(h1 > 0 && h2 > 0);
#endif
            atom_type[a] = 8;
            atom_type[hydro1[0]] = 8;
            in1 = get_bindex(hydro1[0],i);
            btype[4*i+in1] = 2;
            btype[4*hydro1[0]] = 2;
            in1 = get_bindex(a,i);
            btype[4*i+in1] = 1;
            in1 = get_bindex(i,a);
            btype[4*a+in1] = 1;
            drop.insert(hydro2[0]);
            eliminate_atoms(drop);
          }
          return true;
        }
      }
    }
  }
  return false;
}

void Molecule::axial_ring_bonds(std::vector<int>& axial) const
{
  // This method fills the array axial with the atom indices of the axial
  // neighbours for the molecule's ring atoms
  int i,j,k,l,n1,n2,n,a1,a2,temp,rneighbour;
  int rindex,rindex1,rindex2;
  int neighbours1[4],neighbours2[4],neighbours[4];
  bool inring;
  double l1[3],l2[3];

  axial.clear();

  for(l=0; l<nrings; ++l) {
    rindex1 = rings[6*l];
    rneighbour = -1;
    for(i=0; i<4; ++i) {
      temp = bonds[4*rindex1+i];
      if (temp >= 0) {
        for(j=1; j<6; ++j) {
          if (rings[6*l+j] == temp) {
            rneighbour = j;
            break;
          }
        }
      }
      if (rneighbour > 0) break;
    }
    rindex2 = rings[6*l+rneighbour];
    // So now we have a pair of neighbouring ring atoms, rindex1 and rindex2
    // We now need to grab each atom's two non-ring neighbours
    n1 = 0;
    for(i=0; i<4; ++i) {
      temp = bonds[4*rindex1+i];
      if (temp >= 0) {
        inring = false;
        for(j=0; j<6; ++j) {
          if (rings[6*l+j] == temp) inring = true;
        }
        if (!inring) {
          neighbours1[n1] = temp;
          n1++;
        }
      }
    }
    n2 = 0;
    for(i=0; i<4; ++i) {
      temp = bonds[4*rindex2+i];
      if (temp >= 0) {
        inring = false;
        for(j=0; j<6; ++j) {
          if (rings[6*l+j] == temp) inring = true;
        }
        if (!inring) {
          neighbours2[n2] = temp;
          n2++;
        }
      }
    }
#ifdef DEBUG
    assert(n1 == 2 && n2 == 2);
#endif
    // Now we simply need to compare these lines to see which is parallel
    a1 = -1;
    a2 = -1;
    for(i=0; i<2; ++i) {
      for(j=0; j<3; ++j) {
        l1[j] = coords[3*rindex1+j] - coords[3*neighbours1[i]+j];
      }
      for(j=0; j<2; ++j) {
        for(k=0; k<3; ++k) {
          l2[k] = coords[3*rindex2+k] - coords[3*neighbours2[j]+k];
        }
        if (parallel(l1,l2)) {
          a1 = neighbours1[i];
          a2 = neighbours2[j];
          break;
        }
      }
      if (a1 >= 0) break;
    }
    // Why does it not find a parallel pair among these vertices?
    axial.push_back(a1);
    for(i=0; i<3; ++i) {
      l1[i] = coords[3*rindex1+i] - coords[3*a1+i];
    }
    for(i=1; i<6; ++i) {
      rindex = rings[6*l+i];
      if (rindex == rindex2) {
        axial.push_back(a2);
        continue;
      }
      // Grab this ring atom's two non-ring neighbours and see which is parallel to our
      // known axial bond, l1
      n = 0;
      for(j=0; j<4; ++j) {
        temp = bonds[4*rindex+j];
        if (temp >= 0) {
          inring = false;
          for(k=0; k<6; ++k) {
            if (temp == rings[6*l+k]) inring = true;
          }
          if (!inring) {
            neighbours[n] = temp;
            n++;
          }
        }
      }
#ifdef DEBUG
      assert(n == 2);
#endif
      for(j=0; j<2; ++j) {
        for(k=0; k<3; ++k) {
          l2[k] = coords[3*rindex+k] - coords[3*neighbours[j]+k];
        }
        if (parallel(l1,l2)) {
          axial.push_back(neighbours[j]);
          break;
        }
      }
    }
  }
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
          if (atom_type[temp] == 6) {
            // See if it's a methyl
            h = 0;
            for(k=0; k<4; ++k) {
              if (bonds[4*temp+k] >= 0) {
                if (atom_type[bonds[4*temp+k]] == 1) {
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
        atom_type[temp] = 1;
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

void Molecule::normalize_free_ring(int the_aring)
{
  int test,kount,temp,delta,v1,v2,v11,v22,a,a2,in1,done[5];
  int oxy,sul,nit,oxygen,sulfur,nitrogen,rn,rneighbour[2];
  bool found;
  int i,j;
  if (rings[6*the_aring+5] == -1) {
    // We need to see what kind of heteroatom pattern is here:
    // O, S, N, O+N, S+N, N+N
    // The oxygen and sulfur bonds are set to one, or the N with
    // an extra-ring bond
    oxygen = -1;
    sulfur = -1;
    v11 = -1;
    v22 = -1;
    oxy = 0;
    sul = 0;
    nit = 0;
    rn = 0;
    for(i=0; i<5; ++i) {
      if (atom_type[rings[6*the_aring+i]] == 8) oxy++;
      if (atom_type[rings[6*the_aring+i]] == 16) sul++;
      if (atom_type[rings[6*the_aring+i]] == 7) nit++;
    }
#ifdef DEBUG
    assert(oxy < 2 && sul < 2);
#endif
    if (oxy == 1) {
      for(i=0; i<5; ++i) {
        if (atom_type[rings[6*the_aring+i]] == 8) {
          oxygen = rings[6*the_aring+i];
          break;
        }
      }
      for(i=0; i<4; ++i) {
        temp = bonds[4*oxygen+i];
        if (temp < 0) continue;
        if (get_rindex(temp,the_aring) != -1) {
          rneighbour[rn] = i;
          rn++;
        }
      }
#ifdef DEBUG
      assert(rn == 2);
#endif
      btype[4*oxygen+rneighbour[0]] = 1;
      btype[4*oxygen+rneighbour[1]] = 1;
      v1 = bonds[4*oxygen+rneighbour[0]];
#ifdef DEBUG
      assert(get_bindex(oxygen,v1) != -1);
#endif
      btype[4*v1+get_bindex(oxygen,v1)] = 1;
      v2 = bonds[4*oxygen+rneighbour[1]];
      btype[4*v2+get_bindex(oxygen,v2)] = 1;
      // Now the two non-oxygen ring neighbours of v1 and v2
      // will have bond type two:
      for(i=0; i<4; ++i) {
        temp = bonds[4*v1+i];
        if (temp < 0) continue;
        if (temp == oxygen) continue;
        if (get_rindex(temp,the_aring) != -1) v11 = temp;
      }
      btype[4*v1+get_bindex(v11,v1)] = 2;
      btype[4*v11+get_bindex(v1,v11)] = 2;
      for(i=0; i<4; ++i) {
        temp = bonds[4*v2+i];
        if (temp < 0) continue;
        if (temp == oxygen) continue;
        if (get_rindex(temp,the_aring) != -1) v22 = temp;
      }
      btype[4*v2+get_bindex(v22,v2)] = 2;
      btype[4*v22+get_bindex(v2,v22)] = 2;
      // Finally the bond between v11 and v22 needs to be set
      // to be of type one:
      btype[4*v11+get_bindex(v22,v11)] = 1;
      btype[4*v22+get_bindex(v11,v22)] = 1;
    }
    else if (sul == 1) {
      for(i=0; i<5; ++i) {
        if (atom_type[rings[6*the_aring+i]] == 16) {
          sulfur = rings[6*the_aring+i];
          break;
        }
      }
      for(i=0; i<4; ++i) {
        temp = bonds[4*sulfur+i];
        if (temp < 0) continue;
        if (get_rindex(temp,the_aring) != -1) {
          rneighbour[rn] = i;
          rn++;
        }
      }
#ifdef DEBUG
      assert(rn == 2);
#endif
      btype[4*sulfur+rneighbour[0]] = 1;
      btype[4*sulfur+rneighbour[1]] = 1;
      v1 = bonds[4*sulfur+rneighbour[0]];
      btype[4*v1+get_bindex(sulfur,v1)] = 1;
      v2 = bonds[4*sulfur+rneighbour[1]];
      btype[4*v2+get_bindex(sulfur,v2)] = 1;
      // Now the two non-sulfur ring neighbours of v1 and v2
      // will have bond type two:
      for(i=0; i<4; ++i) {
        temp = bonds[4*v1+i];
        if (temp < 0) continue;
        if (temp == sulfur) continue;
        if (get_rindex(temp,the_aring) != -1) v11 = temp;
      }
      btype[4*v1+get_bindex(v11,v1)] = 2;
      btype[4*v11+get_bindex(v1,v11)] = 2;
      for(i=0; i<4; ++i) {
        temp = bonds[4*v2+i];
        if (temp < 0) continue;
        if (temp == sulfur) continue;
        if (get_rindex(temp,the_aring) != -1) v22 = temp;
      }
      btype[4*v2+get_bindex(v22,v2)] = 2;
      btype[4*v22+get_bindex(v2,v22)] = 2;
      // Finally the bond between v11 and v22 needs to be set
      // to be of type one:
      btype[4*v11+get_bindex(v22,v11)] = 1;
      btype[4*v22+get_bindex(v11,v22)] = 1;
    }
    else {
      // Need to find the nitrogen with three bonds
      nitrogen = -1;
      for(i=0; i<5; ++i) {
        if (atom_type[rings[6*the_aring+i]] == 7) {
          test = rings[6*the_aring+i];
          kount = 0;
          for(j=0; j<4; ++j) {
            temp = bonds[4*test+j];
            if (temp >= 0) kount++;
          }
          if (kount == 3) {
            nitrogen = test;
            break;
          }
#ifdef VERBOSE
          std::cout << "Nitrogen has " << kount << " bonds" << std::endl;
#endif
        }
      }
#ifdef DEBUG
      assert(nitrogen != -1);
#endif
      for(i=0; i<4; ++i) {
        temp = bonds[4*nitrogen+i];
        if (temp < 0) continue;
        if (get_rindex(temp,the_aring) != -1) {
          rneighbour[rn] = i;
          rn++;
        }
      }
#ifdef DEBUG
      assert(rn == 2);
#endif
      btype[4*nitrogen+rneighbour[0]] = 1;
      btype[4*nitrogen+rneighbour[1]] = 1;
      v1 = bonds[4*nitrogen+rneighbour[0]];
      btype[4*v1+get_bindex(nitrogen,v1)] = 1;
      v2 = bonds[4*nitrogen+rneighbour[1]];
      btype[4*v2+get_bindex(nitrogen,v2)] = 1;
      // Now the two non-sulfur ring neighbours of v1 and v2
      // will have bond type two:
      for(i=0; i<4; ++i) {
        temp = bonds[4*v1+i];
        if (temp < 0) continue;
        if (temp == nitrogen) continue;
        if (get_rindex(temp,the_aring) != -1) v11 = temp;
      }
      btype[4*v1+get_bindex(v11,v1)] = 2;
      btype[4*v11+get_bindex(v1,v11)] = 2;
      for(i=0; i<4; ++i) {
        temp = bonds[4*v2+i];
        if (temp < 0) continue;
        if (temp == nitrogen) continue;
        if (get_rindex(temp,the_aring) != -1) v22 = temp;
      }
      btype[4*v2+get_bindex(v22,v2)] = 2;
      btype[4*v22+get_bindex(v2,v22)] = 2;
      // Finally the bond between v11 and v22 needs to be set
      // to be of type one:
      btype[4*v11+get_bindex(v22,v11)] = 1;
      btype[4*v22+get_bindex(v11,v22)] = 1;
    }
  }
  else {
    a = rings[6*the_aring];
    done[0] = a;
    rn = 0;
    for(j=0; j<4; ++j) {
      temp = bonds[4*a+j];
      if (temp < 0) continue;
      if (get_rindex(temp,the_aring) != -1) {
        rneighbour[rn] = temp;
        rn++;
      }
    }
    btype[4*a+get_bindex(rneighbour[0],a)] = 1;
    btype[4*a+get_bindex(rneighbour[1],a)] = 2;
    btype[4*rneighbour[0]+get_bindex(a,rneighbour[0])] = 1;
    btype[4*rneighbour[1]+get_bindex(a,rneighbour[1])] = 2;
    done[1] = rneighbour[0];
    done[2] = rneighbour[1];
    for(j=0; j<4; ++j) {
      temp = bonds[4*rneighbour[0]+j];
      if (temp < 0) continue;
      if (temp == a) continue;
      delta = get_rindex(temp,the_aring);
      if (delta != -1) {
        in1 = rings[6*the_aring+delta];
        done[3] = in1;
        btype[4*rneighbour[0]+j] = 2;
        btype[4*in1+get_bindex(rneighbour[0],in1)] = 2;
      }
    }
    for(j=0; j<4; ++j) {
      temp = bonds[4*rneighbour[1]+j];
      if (temp < 0) continue;
      if (temp == a) continue;
      delta = get_rindex(temp,the_aring);
      if (delta != -1) {
        in1 = rings[6*the_aring+delta];
        done[4] = in1;
        btype[4*rneighbour[1]+j] = 1;
        btype[4*in1+get_bindex(rneighbour[1],in1)] = 1;
      }
    }
    // Get the one ring node we haven't done yet
    for(i=0; i<6; ++i) {
      in1 = rings[6*the_aring+i];
      found = false;
      for(j=0; j<5; ++j) {
        if (in1 == done[j]) {
          found = true;
          break;
        }
      }
      if (!found) {
        a2 = rings[6*the_aring+i];
        in1 = get_bindex(done[3],a2);
        btype[4*a2+in1] = 1;
        btype[4*done[3]+get_bindex(a2,done[3])] = 1;
        in1 = get_bindex(done[4],a2);
        btype[4*a2+in1] = 2;
        btype[4*done[4]+get_bindex(a2,done[4])] = 2;
        break;
      }
    }
  }
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

bool Molecule::normalize_safe(const std::vector<int>& aromatic,bool* done)
{
  // This routine will go searching for aromatic bonds to convert from 4 to
  // 1 or 2, but in a very careful manner to make sure that no contradictions
  // arise - if so, return false.
  int i,j,k,l,atype,winner,test1,test2,ratom,atom,atom2,btype1,btype2;
  int valence,nringp,ring_partners[2];
  bool found,change;

  do {
    winner = -1;
    for(i=0; i<(signed) aromatic.size(); ++i) {
      if (done[i]) continue;
      // Now see if this ring has any bonds which have been converted to
      // one or two
      for(j=0; j<6; ++j) {
        ratom = rings[6*aromatic[i]+j];
        if (ratom == -1) continue;
        test1 = -1;
        test2 = -1;
        for(k=0; k<4; ++k) {
          if (bonds[4*ratom+k] < 0) continue;
          found = false;
          for(l=0; l<6; ++l) {
            if (bonds[4*ratom+k] == rings[6*aromatic[i]+l]) {
              found = true;
              break;
            }
          }
          if (found) {
            if (test1 < 0) {
              test1 = btype[4*ratom+k];
            }
            else {
              test2 = btype[4*ratom+k];
            }
          }
        }
        if (test1 == 1 || test1 == 2) {
          winner = i;
          break;
        }
        if (test2 == 1 || test2 == 2) {
          winner = i;
          break;
        }
      }
      if (winner >= 0) break;
    }
    if (winner >= 0) {
      done[winner] = true;
      // Let's find out which bonds have already been done...
      do {
        change = false;
        for(j=0; j<6; ++j) {
          atom = rings[6*aromatic[winner]+j];
          if (atom == -1) continue;
          nringp = 0;
          valence = 0;
          for(k=0; k<4; ++k) {
            if (bonds[4*atom+k] < 0) continue;
            found = false;
            for(l=0; l<6; ++l) {
              if (bonds[4*atom+k] == rings[6*aromatic[winner]+l]) {
                found = true;
                break;
              }
            }
            if (found) {
              ring_partners[nringp] = k;
              nringp++;
            }
            else {
              if (btype[4*atom+k] != 4) {
                valence += btype[4*atom+k];
              }
              else {
                valence += 1;
              }
            }
          }
#ifdef VERBOSE
          std::cout << "For atom " << atom << " type = " << atom_type[atom] << " we have nringp " << nringp << " and valence " << valence << std::endl;
#endif
#ifdef DEBUG
          assert(nringp == 2);
          assert(valence < 3);
#endif
          atype = atom_type[atom];
          btype1 = btype[4*atom+ring_partners[0]];
          btype2 = btype[4*atom+ring_partners[1]];
          if (btype1 == 4) {
            if (btype2 == 1) {
              change = true;
              if (valence == 1 && atype == 6) {
                btype[4*atom+ring_partners[0]] = 2;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 2;
              }
              else if (valence == 1 && atype == 7) {
                btype[4*atom+ring_partners[0]] = 1;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              else if (valence == 2) {
                btype[4*atom+ring_partners[0]] = 1;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              else if (valence == 0 && atype == 7) {
                btype[4*atom+ring_partners[0]] = 2;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 2;
              }
              else if (valence == 0 && (atype == 8 || atype == 16)) {
                btype[4*atom+ring_partners[0]] = 1;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              break;
            }
            else if (btype2 == 2) {
              change = true;
              if (valence == 1) {
                btype[4*atom+ring_partners[0]] = 1;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              else if (valence == 2) {
                return false;
              }
              else if (valence == 0) {
                btype[4*atom+ring_partners[0]] = 1;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              break;
            }
          }
          if (btype2 == 4) {
            if (btype1 == 1) {
              change = true;
              if (valence == 1 && atype == 6) {
                btype[4*atom+ring_partners[1]] = 2;
                atom2 = bonds[4*atom+ring_partners[1]];
                btype[4*atom2+get_bindex(atom,atom2)] = 2;
              }
              else if (valence == 1 && atype == 7) {
                btype[4*atom+ring_partners[1]] = 1;
                atom2 = bonds[4*atom+ring_partners[1]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              else if (valence == 2) {
                btype[4*atom+ring_partners[1]] = 1;
                atom2 = bonds[4*atom+ring_partners[1]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              else if (valence == 0 && atype == 7) {
                btype[4*atom+ring_partners[1]] = 2;
                atom2 = bonds[4*atom+ring_partners[1]];
                btype[4*atom2+get_bindex(atom,atom2)] = 2;
              }
              else if (valence == 0 && (atype == 8 || atype == 16)) {
                btype[4*atom+ring_partners[1]] = 1;
                atom2 = bonds[4*atom+ring_partners[1]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              break;
            }
            else if (btype1 == 2) {
              change = true;
              if (valence == 1) {
                btype[4*atom+ring_partners[0]] = 1;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              else if (valence == 2) {
                return false;
              }
              else if (valence == 0) {
                btype[4*atom+ring_partners[0]] = 1;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              break;
            }
          }
        }
        if (!change) break;
      } while(true);
      // Check the valence of each of the ring atoms to see that it is
      // appropriate
      for(j=0; j<6; ++j) {
        atom = rings[6*aromatic[winner]+j];
        if (atom == -1) continue;
        valence = 0;
        atype = atom_type[atom];
        for(k=0; k<4; ++k) {
          if (bonds[4*atom+k] >= 0) valence += btype[4*atom+k];
        }
        if (atype == 6 && valence != 4) return false;
        if (atype == 7 && valence != 3) return false;
        if (atype == 8 && valence != 2) return false;
        if (atype == 16 && valence != 2) return false;
      }
    }
    else {
      break;
    }
  } while(true);
  return true;
}

void Molecule::connected_components(int vertices)
{
  std::vector<int> visited;
  int i,j,nzero,nc = 1;
  bool first = true;

  for(i=0; i<vertices; ++i) {
    visited.push_back(0);
  }
  nzero = 0;
  do {
    visited[nzero] = nc;
    propagate(visited,nc);
    first = true;
    for(j=0; j<vertices; ++j) {
      if (visited[j] == 0 && first) {
        first = false;
        nzero = j;
      }
    }
    if (first) break;
    nc++;
  } while(true);
  pieces = new std::vector<int>[nc];
  for(i=0; i<vertices; ++i) {
    pieces[visited[i]-1].push_back(i);
  }
  p_allocated = nc;
}

void Molecule::propagate(std::vector<int>& visited,int nc) const
{
  unsigned int i,j;
  std::vector<unsigned int> convert;

  do {
    convert.clear();
    for(i=0; i<visited.size(); ++i) {
      if (visited[i] != 0) continue;
      for(j=0; j<4; ++j) {
        if (rbonds[4*i+j] < 0) continue;
        if (visited[rbonds[4*i+j]] == nc) {
          convert.push_back(i);
          break;
        }
      }
    }
    if (convert.empty()) break;
    for(i=0; i<convert.size(); ++i) {
      visited[convert[i]] = nc;
    }
  } while(true);
}

bool Molecule::normalize_aromatic_bonds()
{
  // This method's job is to convert all these aromatic ring bonds
  // of type "4" to proper "1" or "2" bonds for the SD file output
  // The strategy here is as follows: first look for
  // systems of fused rings. If there are none (i.e. no two
  // rings have at least one edge in common), then we may
  // apply the usual strategy to aromatic rings. For each
  // fused ring cluster, check to see if there are any
  // aromatic rings contained within it, if no, pass, if
  // yes, see if any of the aromatic are five-membered.
  // First some simple escape routes:
  int i,j,k,a,temp,delta,need_to_do,kount,the_ring,fiver,sum,leave = 0,the_aring = -1;
  unsigned int l;
  bool found,test,done[6];
  bool* ring_cluster;
  std::vector<int> aromatic,ratom;

  for(i=0; i<nrings; ++i) {
    if (is_aromatic(i)) {
      leave++;
      the_aring = i;
    }
  }
#ifdef VERBOSE
  std::cout << "Found " << leave << " aromatic rings" << std::endl;
#endif
  if (leave == 0) return true;

  if (leave == 1) {
    // This will be trivial to fix up
    normalize_free_ring(the_aring);
    return true;
  }
  // First need to look for ring clusters: we will remove all
  // non-ring atoms from the molecule, and see what the connected
  // components are.
  for(i=0; i<natoms; ++i) {
    for(j=0; j<nrings; ++j) {
      found = false;
      for(k=0; k<6; ++k) {
        if (rings[6*j+k] == i) {
          found = true;
          break;
        }
      }
      if (found) {
        ratom.push_back(i);
        break;
      }
    }
  }
#ifdef VERBOSE
  std::cout << "Number of ring atoms is " << ratom.size() << std::endl;
#endif
  for(i=0; i<4*natoms; ++i) {
    rbonds.push_back(-1);
  }

  for(l=0; l<ratom.size(); ++l) {
    a = ratom[l];
    kount = 0;
    for(j=0; j<4; ++j) {
      temp = bonds[4*a+j];
      if (temp >= 0) {
        delta = get_index(temp,ratom);
        if (delta != -1) {
          rbonds[4*l+kount] = delta;
          kount++;
        }
      }
    }
  }

  connected_components(ratom.size());

  if (p_allocated == nrings) {
    // This is also quite simple - there are no fused rings at
    // all
    for(i=0; i<nrings; ++i) {
      if (is_aromatic(i)) {
        // This ring needs to have its bonds changed to one and
        // two
        normalize_free_ring(i);
      }
    }
    return true;
  }
  for(i=0; i<p_allocated; ++i) {
    for(l=0; l<pieces[i].size(); ++l) {
      pieces[i][l] = ratom[pieces[i][l]];
    }
  }
  ring_cluster = new bool[p_allocated*nrings];
  for(i=0; i<p_allocated; ++i) {
    for(j=0; j<nrings; ++j) {
      ring_cluster[nrings*i+j] = false;
      found = false;
      for(l=0; l<pieces[i].size(); ++l) {
        if (rings[6*j] == pieces[i][l]) {
          found = true;
          break;
        }
      }
      if (found) {
        ring_cluster[nrings*i+j] = true;
      }
    }
  }
  for(i=0; i<p_allocated; ++i) {
    // We'll go through each ring cluster to see firstly if it
    // a) has more than one aromatic ring and
    // b) one of these is five-membered
    sum = 0;
    for(j=0; j<nrings; ++j) {
      if (ring_cluster[nrings*i+j]) sum++;
    }
#ifdef DEBUG
    assert(sum != 0);
#endif
    if (sum > 1) {
      for(j=0; j<nrings; ++j) {
        if (ring_cluster[nrings*i+j]) {
          if (is_aromatic(j)) aromatic.push_back(j);
        }
      }
      if (aromatic.empty()) {
        continue;
      }
      else if (aromatic.size() == 1) {
        normalize_free_ring(aromatic[0]);
      }
      else {
        // The hard case... see if any of the aromatic rings are
        // five-membered
        for(l=0; l<aromatic.size(); ++l) {
          done[l] = false;
        }
        fiver = -1;
        for(l=0; l<aromatic.size(); ++l) {
          if (rings[6*aromatic[l]+5] == -1) {
            fiver = l;
            break;
          }
        }
        if (fiver == -1) {
          normalize_free_ring(aromatic[0]);
          done[0] = true;
          do {
            test = normalize_safe(aromatic,done);
            if (!test) {
              delete[] ring_cluster;
              return false;
            }
            else {
              need_to_do = -1;
              for(l=0; l<aromatic.size(); ++l) {
                if (!done[l]) {
                  need_to_do = l;
                  break;
                }
              }
              if (need_to_do == -1) {
                break;
              }
              else {
                normalize_free_ring(aromatic[need_to_do]);
                done[need_to_do] = true;
              }
            }
          } while(true);
        }
        else {
          // Take the five-membered ring and do it
          normalize_free_ring(aromatic[fiver]);
          done[fiver] = true;
          // Now go through the others and do them one by one, and if we reach a
          // contradiction, just return zero
          do {
            test = normalize_safe(aromatic,done);
            if (!test) {
              delete[] ring_cluster;
              return false;
            }
            else {
              need_to_do = -1;
              for(l=0; l<aromatic.size(); ++l) {
                if (!done[l]) {
                  need_to_do = l;
                  break;
                }
              }
              if (need_to_do == -1) {
                break;
              }
              else {
                normalize_free_ring(aromatic[need_to_do]);
                done[need_to_do] = true;
              }
            }
          } while(true);
        }
      }
    }
    else {
      the_ring = -1;
      for(j=0; j<nrings; ++j) {
        if (ring_cluster[nrings*i+j] == 1) the_ring = j;
      }
#ifdef DEBUG
      assert(the_ring != -1);
#endif
      if (is_aromatic(the_ring)) normalize_free_ring(the_ring);
    }
  }
  delete[] ring_cluster;
  return true;
}

std::string Molecule::to_MDLMol() const
{
  int i,j,bnumber = 0;
  long seconds = long(std::time(nullptr));
  std::string atom;
  std::ostringstream s;
  std::map<int,std::string>::const_iterator it;

  assert(consistent());

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
    it = element_table.find(atom_type[i]);
#ifdef DEBUG
    assert(it != element_table.end());
#endif
    atom = it->second; 
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
