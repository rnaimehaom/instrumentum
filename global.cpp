#include "global.h"

typedef boost::mt19937 base_generator_type;
base_generator_type generator(42u);
boost::uniform_real<> uni_dist(0,1);
boost::variate_generator<base_generator_type&,boost::uniform_real<> > brandom(generator,uni_dist);

void initialize_generator(unsigned int seed)
{
  unsigned int s = seed*((unsigned int) std::time(NULL));
  generator.seed(s);
}

unsigned int irandom(unsigned int nmax)
{
  unsigned int output = (unsigned) int(double(nmax)*rrandom());
  return output;
}

double rrandom()
{
  double output = double(brandom());
  return output;
}

bool parallel(const double* x,const double* y)
{
  // We take the cross product of these two 3-vectors and see if
  // it's approximately zero
  double delta,out[3];
  out[0] = x[1]*y[2]-x[2]*y[1];
  out[1] = x[2]*y[0]-x[0]*y[2];
  out[2] = x[0]*y[1]-x[1]*y[0];
  delta = std::sqrt(out[0]*out[0] + out[1]*out[1] + out[2]*out[2]);
  if (delta < 0.001) return true;
  return false;
}

void shuffle(std::vector<int>& v)
{
  int test;
  bool good;
  unsigned int i,kount = 0;

  for(i=0; i<v.size(); ++i) {
    v[i] = i;
  }
  do {
    test = irandom(v.size());
    good = true;
    for(i=0; i<kount; ++i) {
      if (v[i] == test) {
        good = false;
        break;
      }
    }
    if (good) {
      v[kount] = test;
      kount++;
      if (kount == v.size()) break;
    }
  } while(true);
}

int get_index(int element,const std::vector<int>& v)
{
  unsigned int i;
  for(i=0; i<v.size(); ++i) {
    if (v[i] == element) return i;
  }
  return -1;
}

bool g_connected(const std::vector<int>& bonds)
{
  unsigned int i,j,n = bonds.size()/4;
  int in1;
  bool output = true;
  bool* visited = new bool[n];
  std::vector<int> convert;

  for(i=0; i<n; ++i) {
    visited[i] = false;
  }

  visited[0] = true;
  do {
    convert.clear();
    for(i=0; i<n; ++i) {
      if (visited[i]) {
        for(j=0; j<4; ++j) {
          in1 = bonds[4*i+j];
          if (in1 >= 0) {
            if (!visited[in1]) convert.push_back(in1);
          }
        }
      }
    }
    if (convert.empty()) break;
    for(i=0; i<convert.size(); ++i) {
      visited[convert[i]] = true;
    }
  } while(true);

  for(i=0; i<n; ++i) {
    if (!visited[i]) {
      output = false;
      break;
    }
  }
  delete[] visited;
  return output;
}

void ename(int anumber,char* element)
{
  switch(anumber) {
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

void rperception(const std::vector<int>& rbonds,std::vector<int>& rings)
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


