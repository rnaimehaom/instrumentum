#include "molecule.h"

bool Molecule::normalize_safe(const std::vector<int>& aromatic,bool* done)
{
  // This routine will go searching for aromatic bonds to convert from 4 to
  // 1 or 2, but in a very careful manner to make sure that no contradictions
  // arise - if so, return false.
  int i,j,k,l,winner,test1,test2,ratom,atom,atom2,btype1,btype2;
  int valence,nringp,ring_partners[2];
  bool found,change;
  Atom_Type atype;

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
          std::cout << "For atom " << atom << " type = " << Molecule::atom_names[static_cast<int>(atom_type[atom])] << " we have nringp " << nringp << " and valence " << valence << std::endl;
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
              if (valence == 1 && atype == Atom_Type::carbon) {
                btype[4*atom+ring_partners[0]] = 2;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 2;
              }
              else if (valence == 1 && atype == Atom_Type::nitrogen) {
                btype[4*atom+ring_partners[0]] = 1;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              else if (valence == 2) {
                btype[4*atom+ring_partners[0]] = 1;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              else if (valence == 0 && atype == Atom_Type::nitrogen) {
                btype[4*atom+ring_partners[0]] = 2;
                atom2 = bonds[4*atom+ring_partners[0]];
                btype[4*atom2+get_bindex(atom,atom2)] = 2;
              }
              else if (valence == 0 && (atype == Atom_Type::oxygen || atype == Atom_Type::sulfur)) {
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
              if (valence == 1 && atype == Atom_Type::carbon) {
                btype[4*atom+ring_partners[1]] = 2;
                atom2 = bonds[4*atom+ring_partners[1]];
                btype[4*atom2+get_bindex(atom,atom2)] = 2;
              }
              else if (valence == 1 && atype == Atom_Type::nitrogen) {
                btype[4*atom+ring_partners[1]] = 1;
                atom2 = bonds[4*atom+ring_partners[1]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              else if (valence == 2) {
                btype[4*atom+ring_partners[1]] = 1;
                atom2 = bonds[4*atom+ring_partners[1]];
                btype[4*atom2+get_bindex(atom,atom2)] = 1;
              }
              else if (valence == 0 && atype == Atom_Type::nitrogen) {
                btype[4*atom+ring_partners[1]] = 2;
                atom2 = bonds[4*atom+ring_partners[1]];
                btype[4*atom2+get_bindex(atom,atom2)] = 2;
              }
              else if (valence == 0 && (atype == Atom_Type::oxygen || atype == Atom_Type::sulfur)) {
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
        if (atype == Atom_Type::carbon && valence != 4) return false;
        if (atype == Atom_Type::nitrogen && valence != 3) return false;
        if (atype == Atom_Type::oxygen && valence != 2) return false;
        if (atype == Atom_Type::sulfur && valence != 2) return false;
      }
    }
    else {
      break;
    }
  } while(true);
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
      if (atom_type[rings[6*the_aring+i]] == Atom_Type::oxygen) oxy++;
      if (atom_type[rings[6*the_aring+i]] == Atom_Type::sulfur) sul++;
      if (atom_type[rings[6*the_aring+i]] == Atom_Type::nitrogen) nit++;
    }
#ifdef DEBUG
    assert(oxy < 2 && sul < 2);
#endif
    if (oxy == 1) {
      for(i=0; i<5; ++i) {
        if (atom_type[rings[6*the_aring+i]] == Atom_Type::oxygen) {
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
        if (atom_type[rings[6*the_aring+i]] == Atom_Type::sulfur) {
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
        if (atom_type[rings[6*the_aring+i]] == Atom_Type::nitrogen) {
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
        if (atom_type[rings[6*i+j]] == Atom_Type::nitrogen) nitro++;
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
              if (atom_type[temp] == Atom_Type::hydrogen) {
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
          atom_type[candidates[alpha1]] = Atom_Type::oxygen;
        }
        else {
          atom_type[candidates[alpha1]] = Atom_Type::sulfur;
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
          if (atom_type[rings[6*i+j]] == Atom_Type::nitrogen) the_nitro = rings[6*i+j];
          h = 0;
          c1 = -1;
          c2 = -1;
          for(k=0; k<4; ++k) {
            temp = bonds[4*rings[6*i+j]+k];
            if (temp >= 0) {
              if (atom_type[temp] == Atom_Type::hydrogen) {
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
        add_atom(Atom_Type::hydrogen,xc,0);
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
          if (atom_type[rings[6*i+j]] == Atom_Type::nitrogen) {
            n_atoms[n] = rings[6*i+j];
            n++;
          }
          h = 0;
          c1 = -1;
          c2 = -1;
          for(k=0; k<4; ++k) {
            temp = bonds[4*rings[6*i+j]+k];
            if (temp >= 0) {
              if (atom_type[temp] == Atom_Type::hydrogen) {
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
        add_atom(Atom_Type::hydrogen,xc,0);
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
        if (atom_type[rings[6*i+j]] != Atom_Type::carbon) continue;
        h = 0;
        c1 = 0;
        c2 = 0;
        for(k=0; k<4; ++k) {
          temp = bonds[4*rings[6*i+j]+k];
          if (temp >= 0) {
            if (atom_type[temp] == Atom_Type::hydrogen) {
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
          std::cout << i << "  " << Molecule::atom_names[static_cast<int>(atom_type[axial[6*i+j]])] << std::endl;
#endif
        }
        if (atom_type[axial[6*i+j]] != Atom_Type::hydrogen) {
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
