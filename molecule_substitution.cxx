#include "molecule.h"

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
    if (atom_type[temp1] == Atom_Type::carbon) {
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
              if (atom_type[tt] != Atom_Type::carbon) hcount++;
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
              if (atom_type[temp2] == Atom_Type::hydrogen) {
                hydrogen[h] = temp2;
                h++;
              }
              else if (atom_type[temp2] != Atom_Type::carbon) {
                exotic = true;
              }
              else if (atom_type[temp2] == Atom_Type::carbon) {
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
                if (atom_type[rings[6*the_ring+j]] != Atom_Type::carbon) continue;
                hh = -1;
                c1 = -1;
                c2 = -1;
                for(k=0; k<4; ++k) {
                  temp = bonds[4*rings[6*the_ring+j]+k];
                  if (temp < 0) continue;
                  if (atom_type[temp] == Atom_Type::hydrogen) {
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
            atom_type[temp1] = Atom_Type::nitrogen;
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
            atom_type[temp1] = Atom_Type::nitrogen;
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
          if (atom_type[temp2] == Atom_Type::hydrogen) {
            hydrogen[h] = temp2;
            h++;
          }
          else if (atom_type[temp2] != Atom_Type::carbon) {
            exotic = true;
          }
        }
        else {
          break;
        }
      }
      if (h >= 1 && !exotic && nb >= 3) {
        // Let's make the conversion...
        atom_type[temp1] = Atom_Type::nitrogen;
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
    if (atom_type[temp1] == Atom_Type::carbon) {
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
              if (atom_type[tt] != Atom_Type::carbon) hcount++;
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
          if (atom_type[temp2] == Atom_Type::hydrogen) {
            hydrogen[h] = temp2;
            h++;
          }
          else if (atom_type[temp2] != Atom_Type::carbon) {
            exotic = true;
          }
        }
        else {
          break;
        }
      }
      if (h >= 2 && !exotic && nb == 4) {
        // Let's make the conversion...
        atom_type[temp1] = Atom_Type::sulfur;
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
    if (atom_type[temp1] == Atom_Type::carbon) {
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
              if (atom_type[tt] != Atom_Type::carbon) hcount++;
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
          if (atom_type[temp2] == Atom_Type::hydrogen) {
            hydrogen[h] = temp2;
            h++;
          }
          else if (atom_type[temp2] != Atom_Type::carbon) {
            exotic = true;
          }
        }
        else {
          break;
        }
      }
      if (h >= 2 && !exotic && nb == 4) {
        // Let's make the conversion...
        atom_type[temp1] = Atom_Type::oxygen;
        drop.insert(hydrogen[0]); drop.insert(hydrogen[1]);
        eliminate_atoms(drop);
        return true;
      }
    }
  }
  return false;
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
    if (atom_type[temp1] == Atom_Type::hydrogen) {
      // Now see if this hydrogen is bonded to a carbon
      temp2 = bonds[4*temp1];
      if (atom_type[temp2] == Atom_Type::carbon) {
        // Do the conversion...
        alpha = RND.irandom(3);
        if (alpha == 0) {
          atom_type[temp1] = Atom_Type::fluorine;
        }
        else if (alpha == 1) {
          atom_type[temp1] = Atom_Type::chlorine;
        }
        else {
          // Phosphorus - so we need to add two hydrogens that will bond to this
          // P atom, since it has a valence of three, like nitrogen.
          // Something broken here with the bonds I think!!!
          atom_type[temp1] = Atom_Type::phosphorus;

          p1[0] = coords[3*temp2];
          p1[1] = coords[3*temp2+1];
          p1[2] = coords[3*temp2+2];
          p2[0] = coords[3*temp1] - coords[3*temp2];
          p2[1] = coords[3*temp1+1] - coords[3*temp2+1];
          p2[2] = coords[3*temp1+2] - coords[3*temp2+2];

          xc[0] = p1[0] + 2.2*p2[0];
          xc[1] = p1[1] + 2.2*p2[1] + 0.3;
          xc[2] = p1[2] + 2.2*p2[2];
          add_atom(Atom_Type::hydrogen,xc,0);
          add_bond(temp1,natoms-1,1);

          xc[0] = p1[0] + 2.2*p2[0];
          xc[1] = p1[1] + 2.2*p2[1] - 0.3;
          xc[2] = p1[2] + 2.2*p2[2];
          add_atom(Atom_Type::hydrogen,xc,0);
          add_bond(temp1,natoms-1,1);
        }
        return true;
      }
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
    if (atom_type[i] == Atom_Type::carbon) {
      candidate = -1;
      for(j=0; j<4; ++j) {
        if (bonds[4*i+j] < 0) continue;
        if (atom_type[j] == Atom_Type::carbon && btype[4*i+j] == 2) {
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
        if (atom_type[bonds[4*i+j]] == Atom_Type::hydrogen) {
          hydro1[h1] = bonds[4*i+j];
          h1++;
        }
        else if (atom_type[bonds[4*i+j]] == Atom_Type::carbon) {
          c1++;
        }
      }
      h2 = 0;
      c2 = 0;
      for(j=0; j<4; ++j) {
        if (bonds[4*a+j] < 0 || bonds[4*a+j] == i) continue;
        if (atom_type[bonds[4*a+j]] == Atom_Type::hydrogen) {
          hydro2[h2] = bonds[4*a+j];
          h2++;
        }
        else if (atom_type[bonds[4*a+j]] == Atom_Type::carbon) {
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
              atom_type[a] = Atom_Type::nitrogen;
              atom_type[hydro1[0]] = Atom_Type::oxygen;
              in1 = get_bindex(hydro1[0],i);
              btype[4*i+in1] = 2;
              btype[4*hydro1[0]] = 2;
              in1 = get_bindex(a,i);
              btype[4*i+in1] = 1;
              in1 = get_bindex(i,a);
              btype[4*a+in1] = 1;
            }
            else {
              atom_type[i] = Atom_Type::nitrogen;
              atom_type[hydro2[0]] = Atom_Type::oxygen;
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
              atom_type[a] = Atom_Type::nitrogen;
              atom_type[hydro1[0]] = Atom_Type::sulfur;
              in1 = get_bindex(hydro1[0],i);
              btype[4*i+in1] = 2;
              btype[4*hydro1[0]] = 2;
              in1 = get_bindex(a,i);
              btype[4*i+in1] = 1;
              in1 = get_bindex(i,a);
              btype[4*a+in1] = 1;
            }
            else {
              atom_type[i] = Atom_Type::nitrogen;
              atom_type[hydro2[0]] = Atom_Type::sulfur;
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
            atom_type[a] = Atom_Type::nitrogen;
            atom_type[hydro1[0]] = Atom_Type::sulfur;
            in1 = get_bindex(hydro1[0],i);
            btype[4*i+in1] = 2;
            btype[4*hydro1[0]] = 2;
            in1 = get_bindex(a,i);
            btype[4*i+in1] = 1;
            in1 = get_bindex(i,a);
            btype[4*a+in1] = 1;
          }
          else if (alpha == 1) {
            atom_type[a] = Atom_Type::nitrogen;
            atom_type[hydro1[0]] = Atom_Type::oxygen;
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
            atom_type[a] = Atom_Type::oxygen;
            atom_type[hydro1[0]] = Atom_Type::oxygen;
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
