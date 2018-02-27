#include "global.h"

void element(int atomic_number,char* element)
{
  switch(atomic_number) {
  case 0:
    element[0] = '\0';
    break;
  case 1:
    element[0] = 'H';
    element[1] = '\0';
    return;
  case 6:
    element[0] = 'C';
    element[1] = '\0';
    return;
  case 7:
    element[0] = 'N';
    element[1] = '\0';
    return;
  case 8:
    element[0] = 'O';
    element[1] = '\0';
    return;
  case 9:
    element[0] = 'F';
    element[1] = '\0';
    return;
  case 15:
    element[0] = 'P';
    element[1] = '\0';
    return;
  case 16:
    element[0] = 'S';
    element[1] = '\0';
    return;
  case 17:
    element[0] = 'C';
    element[1] = 'l';
    element[2] = '\0';
    return;
  case 35:
    element[0] = 'B';
    element[1] = 'r';
    element[2] = '\0';
    return;
  case 47:
    element[0] = 'A';
    element[1] = 'r';
    element[2] = '\0';
    return;
  case 53:
    element[0] = 'I';
    element[1] = '\0';
    return;
  }
}

bool connected(const std::vector<int>& bonds)
{
  int i,j;
  std::set<int> current,next;
  std::set<int>::const_iterator it;
  const int n = bonds.size()/4;
  bool visited[n];

  for(i=0; i<n; ++i) {
    visited[i] = false;
  }
  visited[0] = true;
  current.insert(0);

  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      if (visited[i]) {
        for(j=0; j<4; ++j) {
          if (bonds[4*i+j] < 0) continue;
          if (!visited[bonds[4*i+j]]) next.insert(bonds[4*i+j]);
        }
      }
    }
    if (next.empty()) break;
    for(it=next.begin(); it!=next.end(); ++it) {
      visited[*it] = true;
    }
    current = next;
    next.clear();
  } while(true);

  for(i=0; i<n; ++i) {
    if (!visited[i]) return false;
  }

  return true;
}

void ring_perception(const std::vector<int>& rbonds,std::vector<int>& rings)
{
  unsigned int i,j,k,l,m;
  int in1,last,first,current;
  bool found,deja = false;
  std::vector<int> vlist,nvlist,nvertex,ring,leftover;

  for(i=0; i<rbonds.size()/4; ++i) {
    vlist.clear();
    leftover.clear();
    nvlist.clear();
    nvlist.push_back(i);
    for(j=0; j<5; ++j) {
      current = nvlist[j];
      nvertex.clear();
      for(k=0; k<4; ++k) {
        in1 = rbonds[4*current+k];
        if (in1 == -1) break;
        found = false;
        for(l=0; l<nvlist.size(); ++l) {
          if (in1 == nvlist[l]) {
            found = true;
            break;
          }
        }
        if (!found) nvertex.push_back(in1);
      }
      if (nvertex.size() == 1) {
        // The simple case with only one option
        nvlist.push_back(nvertex[0]);
      }
      else {
        // Some kind of fusion or junction vertex where
        // several rings come together
        for(k=0; k<nvertex.size()-1; ++k) {
          for(l=0; l<nvlist.size(); ++l) {
            leftover.push_back(nvlist[l]);
          }
          leftover.push_back(nvertex[k]);
          for(l=0; l<6-(1+nvlist.size()); ++l) {
            leftover.push_back(-1);
          }
        }
        nvlist.push_back(nvertex[nvertex.size()-1]);
      }
    }
    for(j=0; j<6; ++j) {
      vlist.push_back(nvlist[j]);
    }
    for(j=0; j<leftover.size()/6; ++j) {
      nvlist.clear();
      for(k=0; k<6; ++k) {
        in1 = 6*j + k;
        if (leftover[in1] >= 0) nvlist.push_back(leftover[in1]);
      }
      for(k=nvlist.size()-1; k<5; ++k) {
        current = nvlist[k];
        nvertex.clear();
        for(l=0; l<4; ++l) {
          in1 = rbonds[4*current+l];
          if (in1 == -1) break;
          found = false;
          for(m=0; m<nvlist.size(); ++m) {
            if (in1 == nvlist[m]) {
              found = true;
              break;
            }
          }
          if (!found) nvertex.push_back(in1);
        }
        if (nvertex.size() == 1) {
          // The simple case with only one option
          nvlist.push_back(nvertex[0]);
        }
        else {
          // Some kind of fusion or junction vertex where
          // several rings come together
          for(l=0; l<nvertex.size()-1; ++l) {
            for(m=0; m<nvlist.size(); ++m) {
              leftover.push_back(nvlist[m]);
            }
            leftover.push_back(nvertex[l]);
            for(m=0; m<6-(1+nvlist.size()); ++m) {
              leftover.push_back(-1);
            }
          }
          nvlist.push_back(nvertex[nvertex.size()-1]);
        }
      }
      for(k=0; k<6; ++k) {
        vlist.push_back(nvlist[k]);
      }
    }
    for(j=0; j<vlist.size()/6; ++j) {
      last = vlist[6*j+5];
      first = vlist[6*j];
      found = false;
      for(k=0; k<4; ++k) {
        in1 = rbonds[4*last+k];
        if (in1 == -1) continue;
        if (in1 == first) {
          found = true;
          break;
        }
      }
      if (found) {
        // It's a ring
        ring.clear();
        for(k=0; k<6; ++k) {
          ring.push_back(vlist[6*j+k]);
        }
        std::sort(ring.begin(),ring.end());
        for(k=0; k<rings.size()/6; ++k) {
          deja = true;
          for(l=0; l<6; ++l) {
            if (rings[6*k+l] != ring[l]) {
              deja = false;
              break;
            }
          }
          if (deja) break;
        }
        if (!deja) {
          for(k=0; k<6; ++k) {
            rings.push_back(ring[k]);
          }
        }
      }
    }
  }
}


