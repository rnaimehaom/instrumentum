#include "grid.h"

int Grid::index1(int i,int j,int k) const
{
  int output = (k+D3) + (1+2*D3)*(j+D2) + (1+2*D3)*(1+2*D2)*(i+D1);
  return output;
}

void Grid::next_door(int m,int n,int k,int i,const double* deltas)
{
  int in1,in2,x_new,y_new,z_new,n_state;
  double x_c,y_c,z_c,dx,dy1,dy2,dz1,dz2;

  dx = deltas[0];
  dy1 = deltas[1];
  dy2 = deltas[2];
  dz1 = deltas[3];
  dz2 = deltas[4];

  in1 = index1(m,n,k);
  x_c = nodes[in1].x;
  y_c = nodes[in1].y;
  z_c = nodes[in1].z;
  n_state = nodes[in1].state;
  x_new = 0;
  y_new = 0;
  z_new = 0;

  if (i == 1) {
    if (n_state == 1) {
      x_new = m - 1;
      x_c = x_c - dx;
      y_new = n;
      y_c = y_c + dy1;
      z_new = k;
      z_c = z_c - dz1;
      n_state = 2;
    }
    else if (n_state == 2) {
      x_new = m + 1;
      x_c = x_c + dx;
      y_new = n;
      y_c = y_c - dy1;
      z_new = k;
      z_c = z_c + dz1;
      n_state = 1;
    }
    else if (n_state == 3) {
      x_new = m - 1;
      x_c = x_c - dx;
      y_new = n;
      y_c = y_c - dy1;
      z_new = k;
      z_c = z_c + dz1;
      n_state = 4;
    }
    else if (n_state == 4) {
      x_new = m + 1;
      x_c = x_c + dx;
      y_new = n;
      y_c = y_c + dy1;
      z_new = k;
      z_c = z_c - dz1;
      n_state = 3;
    }
  }
  else if (i == 2) {
    if (n_state == 1) {
      x_new = m + 1;
      x_c = x_c + dx;
      y_new = n;
      y_c = y_c + dy1;
      z_new = k;
      z_c = z_c - dz1;
      n_state = 2;
    }
    else if (n_state == 2) {
      x_new = m - 1;
      x_c = x_c - dx;
      y_new = n;
      y_c = y_c - dy1;
      z_new = k;
      z_c = z_c + dz1;
      n_state = 1;
    }
    else if (n_state == 3) {
      x_new = m + 1;
      x_c = x_c + dx;
      y_new = n;
      y_c = y_c - dy1;
      z_new = k;
      z_c = z_c + dz1;
      n_state = 4;
    }
    else if (n_state == 4) {
      x_new = m - 1;
      x_c = x_c - dx;
      y_new = n;
      y_c = y_c + dy1;
      z_new = k;
      z_c = z_c - dz1;
      n_state = 3;
    }
  }
  else if (i == 3) {
    if (n_state == 1) {
      x_new = m;
      y_new = n - 1;
      y_c = y_c - dy2;
      z_new = k;
      z_c = z_c - dz1;
      n_state = 2;
    }
    else if (n_state == 2) {
      x_new = m;
      y_new = n + 1;
      y_c = y_c + dy2;
      z_new = k;
      z_c = z_c + dz1;
      n_state = 1;
    }
    else if (n_state == 3) {
      x_new = m;
      y_new = n + 1;
      y_c = y_c + dy2;
      z_new = k;
      z_c = z_c + dz1;
      n_state = 4;
    }
    else if (n_state == 4) {
      x_new = m;
      y_new = n - 1;
      y_c = y_c - dy2;
      z_new = k;
      z_c = z_c - dz1;
      n_state = 3;
    }
  }
  else if (i == 4) {
    if (n_state == 1) {
      x_new = m;
      y_new = n;
      z_new = k + 1;
      z_c = z_c + dz2;
      n_state = 3;
    }
    else if (n_state == 2) {
      x_new = m;
      y_new = n;
      z_new = k - 1;
      z_c = z_c - dz2;
      n_state = 4;
    }
    else if (n_state == 3) {
      x_new = m;
      y_new = n;
      z_new = k - 1;
      z_c = z_c - dz2;
      n_state = 1;
    }
    else if (n_state == 4) {
      x_new = m;
      y_new = n;
      z_new = k + 1;
      z_c = z_c + dz2;
      n_state = 2;
    }
  }

  in2 = index1(x_new,y_new,z_new);
  if (in2 >= 0 && in2 < (signed) total) {
    nodes[in2].x = x_c;
    nodes[in2].y = y_c;
    nodes[in2].z = z_c;
    nodes[in2].state = n_state;
    nodes[in1].neighbours.push_back(x_new);
    nodes[in1].neighbours.push_back(y_new);
    nodes[in1].neighbours.push_back(z_new);
  }
}

Grid::~Grid()
{
  delete[] nodes;
}

Grid::Grid(unsigned int np)
{
  set_default_values();

  allocate(np);
  initialize();
}

Grid::Grid(double blength,unsigned int np)
{
  set_default_values();

  bond_length = blength;

  allocate(np);
  initialize();
}

Grid::Grid(int i,int j,int k,unsigned int np)
{
  set_default_values();

  D1 = i;
  D2 = j;
  D3 = k;

  allocate(np);
  initialize();
}

Grid::Grid(int i,int j,int k,double blength,unsigned int np)
{
  set_default_values();

  D1 = i;
  D2 = j;
  D3 = k;
  bond_length = blength;

  allocate(np);
  initialize();
}

void Grid::set_default_values()
{
  D1 = 17;
  D2 = 17;
  D3 = 9;
  rs1 = 0;
  rs2 = 0;
  rs3 = 0;
  bond_length = 1.4;
}

void Grid::allocate(unsigned int np)
{
  unsigned int i;
  for(i=0; i<np; ++i) {
    pnodes.push_back(-1);
  }
  total = (unsigned) (2*D3+1)*(2*D2+1)*(2*D1+1);
  nodes = new Node[total];
}

void Grid::initialize()
{
  int i,j,k,m,n,in1,in2,state;
  double deltas[5];
  std::vector<int> v;
  unsigned int l;
  const double dx = 0.81635*bond_length;
  const double dy1 = 0.4713*bond_length;
  const double dy2 = 0.9426*bond_length;
  const double dz1 = 0.3338*bond_length;
  const double dz2 = bond_length;

  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        nodes[in1].neighbours.reserve(12);
        nodes[in1].locale = 0;
        nodes[in1].state = 0;
        nodes[in1].atomic_number = 0;
        nodes[in1].x = 0.0;
        nodes[in1].y = 0.0;
        nodes[in1].z = 0.0;
      }
    }
  }

  state = 1;

  k = 0;
  for(i=-D1+1; i<D1; i=i+2) {
    for(j=-D2+1; j<D2; j=j+2) {
      in1 = index1(i,j,k);
      nodes[in1].x = dx*double(i);
      nodes[in1].y = (dy1 + dy2)*double(j);
      nodes[in1].z = (dz1 + dz2)*double(k);
      nodes[in1].state = state;
    }
  }
  for(i=-D1; i<=D1; i=i+2) {
    for(j=-D2; j<=D2; j=j+2) {
      in1 = index1(i,j,k);
      nodes[in1].x = dx*double(i);
      nodes[in1].y = (dy1 + dy2)*double(j);
      nodes[in1].z = (dz1 + dz2)*double(k);
      nodes[in1].state = state;
    }
  }

  deltas[0] = dx;
  deltas[1] = dy1;
  deltas[2] = dy2;
  deltas[3] = dz1;
  deltas[4] = dz2;

  for(n=-D2+4; n<D2-3; n=n+2) {
    for(m=-D1; m<D1+1; m=m+2) {
      for(i=1; i<5; ++i) {
        next_door(m,n,k,i,deltas);
      }
    }
  }
  for(n=-D2+3; n<D2-2; n=n+2) {
    for(m=-D1+1; m<D1; m=m+2) {
      for(i=1; i<5; ++i) {
        next_door(m,n,k,i,deltas);
      }
    }
  }
  for(n=-D2+3; n<D2-2; n=n+2) {
    for(m=-D1; m<D1+1; m=m+2) {
      for(i=1; i<5; ++i) {
        next_door(m,n,k,i,deltas);
      }
    }
  }
  for(n=-D2+4; n<D2-3; n=n+2) {
    for(m=-D1+1; m<D1; m=m+2) {
      for(i=1; i<5; ++i) {
        next_door(m,n,k,i,deltas);
      }
    }
  }

  for(k=-1; k>-(D3+1); k--) {
    if (abs(k)%2 == 1) {
      for(n=-D2+3; n<D2-2; n=n+2) {
        for(m=-D1; m<1+D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+4; n<D2-3; n=n+2) {
        for(m=-D1+1; m<D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+4; n<D2-3; n=n+2) {
        for(m=-D1; m<1+D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+3; n<D2-2; n=n+2) {
        for(m=-D1+1; m<D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
    }
    else {
      for(n=-D2+4; n<D2-3; n=n+2) {
        for(m=-D1; m<1+D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+3; n<D2-2; n=n+2) {
        for(m=-D1+1; m<D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+3; n<D2-2; n=n+2) {
        for(m=-D1; m<1+D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+4; n<D2-3; n=n+2) {
        for(m=-D1+1; m<D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
    }
  }
  for(k=1; k<D3+1; ++k) {
    if ((k%2) == 1) {
      for(n=-D2+4; n<D2-3; n=n+2) {
        for(m=-D1; m<1+D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+3; n<D2-2; n=n+2) {
        for(m=-D1+1; m<D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+3; n<D2-2; n=n+2) {
        for(m=-D1; m<1+D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+4; n<D2-3; n=n+2) {
        for(m=-D1+1; m<D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
    }
    else {
      for(n=-D2+3; n<D2-2; n=n+2) {
        for(m=-D1; m<1+D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+4; n<D2-3; n=n+2) {
        for(m=-D1+1; m<D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+4; n<D2-3; n=n+2) {
        for(m=-D1; m<1+D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
      for(n=-D2+3; n<D2-2; n=n+2) {
        for(m=-D1+1; m<D1; m=m+2) {
          for(i=1; i<5; ++i) {
            next_door(m,n,k,i,deltas);
          }
        }
      }
    }
  }

  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        v = nodes[in1].neighbours;
        nodes[in1].neighbours.clear();
        for(l=0; l<v.size(); l+=3) {
          in2 = index1(v[l],v[l+1],v[l+2]);
          nodes[in1].neighbours.push_back(in2);
        }
      }
    }
  }
}

void Grid::process_pharmacophore(const char* filename)
{


}

void Grid::blank_pharmacophore(double radius)
{
  double delta,x[3];
  int i,j,k,in1,in2,temp;
  unsigned int test,m,l,kount = 0;
  bool bad,done;

  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        x[0] = nodes[in1].x;
        x[1] = nodes[in1].y;
        x[2] = nodes[in1].z;
        delta = std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        if (delta <= radius) nodes[in1].locale = 2;
      }
    }
  }
  // Now that the interior nodes have been marked, we can
  // mark frontier nodes, namely locale=0 nodes that are
  // neighbours of a locale=2 node.
  std::vector<int> frontier_nodes;
  frontier_nodes.reserve(100);
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].locale == 0) {
          done = false;
          for(l=0; l<nodes[in1].neighbours.size(); ++l) {
            in2 = nodes[in1].neighbours[l];
            if (nodes[in2].locale == 2) {
              nodes[in1].locale = 1;
              frontier_nodes.push_back(in1);
              done = true;
              break;
            }
            if (done) break;
          }
        }
      }
    }
  }
  // Now choose two or three of these frontier nodes to be
  // "must have" nodes, to aid in constructing connected
  // molecules
  while(kount < pnodes.size()) {
    test = irandom(frontier_nodes.size());
    bad = false;
    for(m=0; m<kount; ++m) {
      for(l=0; l<nodes[pnodes[m]].neighbours.size(); ++l) {
        in2 = nodes[pnodes[m]].neighbours[l];
        if (in2 == frontier_nodes[test]) {
          bad = true;
          break;
        }
      }
      if (bad) break;
    }
    if (!bad) {
      pnodes[kount] =  frontier_nodes[test];
      nodes[frontier_nodes[test]].locale = 5;
      kount++;
    }
  }
  // Lastly, set the values for rs1, rs2 and rs3
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        temp = nodes[in1].locale;
        if (temp == 1 || temp == 2) {
          if (std::abs(i) > rs1) rs1 = std::abs(i);
          if (std::abs(j) > rs2) rs2 = std::abs(j);
          if (std::abs(k) > rs3) rs3 = std::abs(k);
        }
      }
    }
  }
}

void Grid::fill_interior()
{
  int i,j,k,in1,kount = 0;

  for(i=-rs1; i<=rs1; ++i) {
    for(j=-rs2; j<=rs2; ++j) {
      for(k=-rs3; k<=rs3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].locale == 2) {
          nodes[in1].atomic_number = 6;
          kount++;
        }
      }
    }
  }
#ifdef VERBOSE
  std::cout << "Created " << kount << " carbon nodes" << std::endl;
#endif
}

double Grid::distance(int in1,int in2) const
{
  double x[3],y[3],output;
  x[0] = nodes[in1].x;
  x[1] = nodes[in1].y;
  x[2] = nodes[in1].z;
  y[0] = nodes[in2].x;
  y[1] = nodes[in2].y;
  y[2] = nodes[in2].z;
  output = std::sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]));
  return output;
}

bool Grid::connected()
{
  // To do this, we're first going to paint all of the carbon nodes,
  // and the pharmacophore nodes, as unvisited, and then choose one
  // at random. We then jump outwards to all of its carbon neighbours,
  // and so on, marking each in turn as visited. When we have marched
  // as much as possible, we check to see if any of the pharmacophore
  // nodes are unvisited. If yes, return 0, otherwise, mark up all of
  // unvisited carbon nodes as silver and return 1.
  int i,j,k,in1,in2;
  unsigned int l;
  std::vector<int> convert;

  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number == 6 || nodes[in1].locale == 5) {
          nodes[in1].visited = 0;
        }
        else {
          nodes[in1].visited = -1;
        }
      }
    }
  }

  // Now choose one of the pharmacophoric nodes to begin the process...
  nodes[pnodes[0]].visited = 1;
  do {
    convert.clear();
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].visited == 0) {
            // See if any of my neighbours have been visited
            for(l=0; l<nodes[in1].neighbours.size(); ++l) {
              in2 = nodes[in1].neighbours[l];
              if (nodes[in2].visited == 1) {
                convert.push_back(in1);
                break;
              }
            }
          }
        }
      }
    }
    if (convert.empty()) break;
    for(l=0; l<convert.size(); ++l) {
      nodes[convert[l]].visited = 1;
    }
  } while(true);
  // Last step is to see if all of the pnodes have visited set equal to one,
  // and if so, to set any unvisited carbons to also be silver
  for(l=0; l<pnodes.size(); ++l) {
    if (nodes[pnodes[l]].visited != 1) return false;
  }
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number == 6 && nodes[in1].visited == 0) nodes[in1].atomic_number = 47;
      }
    }
  }
  return true;
}

bool Grid::initial_deletion(double percent,unsigned int max_attempts)
{
  // This function will make a series of random deletions in the
  // initial volume of carbon atoms, checking to preserve the
  // connectivity at each stage, until either "percent" percentage
  // of the original number of heavy atoms remain, or the number
  // "max_attempts" of deletions has been exhausted
  int i,j,k,alpha1,alpha2,alpha3,target,in1;
  unsigned int nheavy,kount,nkill,iterations = 0;
  double radius,delta,n_debut;

  nheavy = 0;
  for(i=-rs1; i<=rs1; ++i) {
    for(j=-rs2; j<=rs2; ++j) {
      for(k=-rs3; k<=rs3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].locale == 2) nheavy++;
      }
    }
  }
  n_debut = double(nheavy);

#ifdef VERBOSE
  std::cout << "Beginning deletions with " << n_debut << " carbon atoms" << std::endl;
#endif

  if (!connected()) {
#ifdef VERBOSE
    std::cout << "Bad initial connectivity..." << std::endl;
#endif
    return false;
  }
  do {
    alpha1 = -rs1 + irandom(2*rs1);
    alpha2 = -rs2 + irandom(2*rs2);
    alpha3 = -rs3 + irandom(2*rs3);
    target = index1(alpha1,alpha2,alpha3);
    radius = 1.0 + rrandom();
    nkill = 0;
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].locale != 2 || nodes[in1].atomic_number == 0) continue;
          delta = distance(in1,target);
          if (delta <= radius) {
            nodes[in1].atomic_number = 47;
            nkill++;
          }
        }
      }
    }
    // Now let us see if eliminating these pseudo-silver atoms will cause
    // the molecule to become disconnected...
    if (nkill == 0) continue;
#ifdef VERBOSE
    std::cout << "Checking connectivity..." << std::endl;
#endif
    if (connected()) {
      // It's alright to get rid of the silver...
      kount = 0;
      for(i=-rs1; i<=rs1; ++i) {
        for(j=-rs2; j<=rs2; ++j) {
          for(k=-rs3; k<=rs3; ++k) {
            in1 = index1(i,j,k);
            if (nodes[in1].atomic_number == 47) {
              nodes[in1].atomic_number = 0;
              kount++;
            }
          }
        }
      }
      nheavy -= kount;
#ifdef VERBOSE
      std::cout << "Good deletion, there are now " << nheavy << " carbon atoms" << std::endl;
#endif
      if ((double(nheavy)/n_debut) <= percent) return true;
    }
    else {
      // It's not alright - change the pseudo-silver back to carbon
      for(i=-rs1; i<=rs1; ++i) {
        for(j=-rs2; j<=rs2; ++j) {
          for(k=-rs3; k<=rs3; ++k) {
            in1 = index1(i,j,k);
            if (nodes[in1].atomic_number == 47) nodes[in1].atomic_number = 6;
          }
        }
      }
    }
    iterations++;
  } while(iterations < max_attempts);
  return false;
}

bool Grid::rationalize(double percent_methyl,unsigned int rings_min,unsigned int rings_max)
{
  // This function eliminates flagpole interactions among
  // hydrogen atoms by forcing any node which is a neighbour
  // of two or more carbon atoms to become carbon as well
  int i,j,k,in1,in2;
  unsigned int alpha,nrings,l,ncarbon,nkill,kount = 0;
  std::vector<int> methyl,convert;

  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number > 0 && nodes[in1].locale == 2) {
          // See if exactly one neighbour is carbon
          ncarbon = 0;
          for(l=0; l<nodes[in1].neighbours.size(); ++l) {
            in2 = nodes[in1].neighbours[l];
            if (nodes[in2].atomic_number == 6) ncarbon++;
          }
          if (ncarbon == 1) methyl.push_back(in1);
        }
      }
    }
  }
  nkill = (unsigned) int(percent_methyl*double(methyl.size()));

  if (nkill > 0) {
    do {
      alpha = irandom(methyl.size());
      in1 = methyl[alpha];
      if (nodes[in1].atomic_number > 0) {
        nodes[in1].atomic_number = 0;
        kount++;
      }
      if (kount == nkill) break;
    } while(true);
  }

  do {
    convert.clear();
    for(i=-D1; i<=D1; ++i) {
      for(j=-D2; j<=D2; ++j) {
        for(k=-D3; k<=D3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].atomic_number == 0) {
            // See if more than one neighbour is carbon
            ncarbon = 0;
            for(l=0; l<nodes[in1].neighbours.size(); ++l) {
              in2 = nodes[in1].neighbours[l];
              if (nodes[in2].atomic_number == 6) ncarbon++;
            }
            if (ncarbon > 1) convert.push_back(in1);
          }
        }
      }
    }
    if (convert.empty()) break;
#ifdef VERBOSE
    std::cout << "Converting " << convert.size() << " nodes to carbon" << std::endl;
#endif
    for(l=0; l<convert.size(); ++l) {
      nodes[convert[l]].atomic_number = 6;
    }
  } while(true);
  // Now the ring counting...
#ifdef DEBUG
  assert(connected());
#endif
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number == 47) nodes[in1].atomic_number = 0;
      }
    }
  }
  nrings = ring_count();
#ifdef VERBOSE
  std::cout << "Found " << nrings << " rings" << std::endl;
#endif
  if (nrings < rings_min || nrings > rings_max) return false;
  return true;
}

void Grid::add_hydrogens()
{
  // Any empty nodes that are neighbours of a carbon atom
  // need to be converted to hydrogen
  int i,j,k,in1,in2;
  unsigned int l;
  std::vector<int> convert;
  bool carbon;

  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number == 0) {
          // See if at least one neighbour is carbon
          carbon = false;
          for(l=0; l<nodes[in1].neighbours.size(); ++l) {
            in2 = nodes[in1].neighbours[l];
            if (nodes[in2].atomic_number == 6) {
              carbon = true;
              break;
            }
          }
          if (carbon) convert.push_back(in1);
        }
      }
    }
  }
#ifdef VERBOSE
  std::cout << "Converting " << convert.size() << " nodes to H" << std::endl;
#endif
  for(l=0; l<convert.size(); ++l) {
    nodes[convert[l]].atomic_number = 1;
  }
}

unsigned int Grid::ring_count()
{
  unsigned int nring = ring_analysis();
  return nring;
}

unsigned int Grid::ring_analysis()
{
  // The first item of business in this routine is to find all
  // of the edges (and so vertices) that are involved in rings.
  // We will do this by removing a bond from the molecule, and
  // seeing if it remains connected.
  int i,j,k,in1,in2;
  unsigned int l,nbonds;
  bool found;
  std::vector<int> vertices,ring_vertices,ring_edges,rings,redges,bonds,wcopy;

  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number > 0 || nodes[in1].locale == 5) vertices.push_back(in1);
      }
    }
  }

  const int nvertices = (signed) vertices.size();
  for(i=0; i<nvertices; ++i) {
    for(j=0; j<4; ++j) {
      bonds.push_back(-1);
    }
  }

  for(i=0; i<nvertices; ++i) {
    nbonds = 0;
    for(l=0; l<nodes[vertices[i]].neighbours.size(); ++l) {
      in2 = nodes[vertices[i]].neighbours[l];
      if (nodes[in2].atomic_number > 0 || nodes[in2].locale == 5) {
        bonds[4*i+nbonds] = get_index(in2,vertices);
        nbonds++;
      }
    }
    if (nbonds == 0) {
#ifdef VERBOSE
      std::cout << "Problem: isolated atom at " << vertices[i] << " with locale " << nodes[vertices[i]].locale << '\n';
#endif
    }
  }
  wcopy = bonds;

  // A sanity check...
  if (!g_connected(bonds)) {
#ifdef VERBOSE
    std::cout << "Error raised..." << std::endl;
#endif
    std::exit(1);
  }

  // We will methodically go through this bond table, and eliminate
  // a bond at each occasion, checking to see if the resulting molecule
  // is still connected
  for(i=0; i<nvertices; ++i) {
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
      if (g_connected(wcopy)) {
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

  if (ring_vertices.size() == 0) return 0;

  for(l=0; l<4*ring_vertices.size(); ++l) {
    redges.push_back(-1);
  }
  for(l=0; l<ring_edges.size(); l+=2) {
    in1 = get_index(ring_edges[l],ring_vertices);
    in2 = get_index(ring_edges[l+1],ring_vertices);
#ifdef DEBUG
    assert(in1 >= 0 && in1 < (signed) (4*ring_vertices.size()));
    assert(in2 >= 0 && in1 < (signed) (4*ring_vertices.size()));
#endif
    for(j=0; j<4; ++j) {
      if (redges[4*in1+j] == -1) {
        redges[4*in1+j] = in2;
        break;
      }
    }
    for(j=0; j<4; ++j) {
      if (redges[4*in2+j] == -1) {
        redges[4*in2+j] = in1;
        break;
      }
    }
  }
  rperception(redges,rings);
  ring_info.clear();
  for(l=0; l<rings.size(); ++l) {
    ring_info.push_back(vertices[ring_vertices[rings[l]]]);
  }
  return ring_info.size()/6;
}

bool Grid::secondary_deletion(unsigned int nc4,unsigned int nc4rings,unsigned int nrings,unsigned int attempts)
{
  int i,j,k,alpha1,alpha2,alpha3,target,in1,in2,temp;
  unsigned int m,l,q,nkill,cring,ring_count,ring_member,ncarbon,silver,nbonds,iterations = 0;
  std::vector<int> vertices,c4;
  double radius,delta;

  do {
    alpha1 = -rs1 + irandom(2*rs1);
    alpha2 = -rs2 + irandom(2*rs2);
    alpha3 = -rs3 + irandom(2*rs3);
    target = index1(alpha1,alpha2,alpha3);
    radius = 2.0 + 2.5*rrandom();
    nkill = 0;
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].locale > 2 || nodes[in1].atomic_number == 0) continue;
          delta = distance(in1,target);
          if (delta <= radius) {
            nodes[in1].atomic_number = 0;
            nkill++;
          }
        }
      }
    }
    if (nkill > 0) {
      // Let's make sure it worked smoothly...
      silver = 0;
      vertices.clear();
      for(i=-D1; i<=D1; ++i) {
        for(j=-D2; j<=D2; ++j) {
          for(k=-D3; k<=D3; ++k) {
            in1 = index1(i,j,k);
            if (nodes[in1].atomic_number > 0 || nodes[in1].locale == 5) vertices.push_back(in1);
          }
        }
      }
      for(m=0; m<vertices.size(); ++m) {
        nbonds = 0;
        for(l=0; l<nodes[vertices[m]].neighbours.size(); ++l) {
          in2 = nodes[vertices[m]].neighbours[l];
          if (nodes[in2].atomic_number > 0 || nodes[in2].locale == 5) nbonds++;
        }
        if (nbonds == 0) {
          silver++;
          nodes[vertices[m]].atomic_number = 0;
        }
      }
#ifdef DEBUG
      assert(connected());
#endif
      silver = 0;
      for(i=-rs1; i<=rs1; ++i) {
        for(j=-rs2; j<=rs2; ++j) {
          for(k=-rs3; k<=rs3; ++k) {
            in1 = index1(i,j,k);
            if (nodes[in1].atomic_number == 47) {
              nodes[in1].atomic_number = 0;
              silver++;
            }
          }
        }
      }
      // First, see how many carbons with four carbon neighbours we have
      c4.clear();
      for(i=-rs1; i<=rs1; ++i) {
        for(j=-rs2; j<=rs2; ++j) {
          for(k=-rs3; k<=rs3; ++k) {
            in1 = index1(i,j,k);
            if (nodes[in1].atomic_number > 0) {
              // See if more than neighbour is carbon
              ncarbon = 0;
              for(l=0; l<nodes[in1].neighbours.size(); ++l) {
                in2 = nodes[in1].neighbours[l];
                if (nodes[in2].atomic_number == 6) ncarbon++;
              }
              if (ncarbon == 4) c4.push_back(in1);
            }
          }
        }
      }
      // Next, see how many of these four carbon neighbours are in four separate
      // rings
      cring = 0;
      ring_count = ring_analysis();
      for(l=0; l<c4.size(); ++l) {
        temp = c4[l];
        ring_member = 0;
        for(m=0; m<ring_count; ++m) {
          for(q=0; q<6; ++q) {
            if (ring_info[6*m+q] == temp) {
              ring_member++;
              break;
            }
          }
        }
        if (ring_member == 4) cring++;
      }
      if (c4.size() <= nc4 && cring <= nc4rings && ring_count <= nrings) return true;
    }
    iterations++;
  } while(iterations < attempts);

  return false;
}

bool Grid::path_selection(bool random)
{
  int i,j,k,inode,temp,pcount,in1,in2,kount,current_pt,next_node;
  unsigned int l,m,nbonds;
  bool found,done;
  std::vector<int> candidate,pathback,convert,winner,vertices,bonds;

  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        nodes[in1].path_hop = -1;
      }
    }
  }

  if (random) {
    inode = pnodes[irandom(pnodes.size())];
  }
  else {
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].atomic_number > 0) candidate.push_back(in1);
        }
      }
    }
    inode = candidate[irandom(candidate.size())];
  }
  nodes[inode].path_hop = 0;

  kount = 0;
  do {
    kount++;
    convert.clear();
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].path_hop < 0) {
            if (nodes[in1].atomic_number > 0 || nodes[in1].locale == 5) {
              for(l=0; l<nodes[in1].neighbours.size(); ++l) {
                in2 = nodes[in1].neighbours[l];
                if (nodes[in2].path_hop >= 0) {
                  found = false;
                  for(m=0; m<convert.size(); ++m) {
                    if (convert[m] == in1) {
                      found = true;
                      break;
                    }
                  }
                  if (!found) {
                    convert.push_back(in1);
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
    if (convert.size() == 0) break;
    for(i=0; i<(signed) convert.size(); ++i) {
      nodes[convert[i]].path_hop = kount;
    }
  } while(true);

  winner.reserve(4);
  for(m=0; m<pnodes.size(); ++m) {
    if (pnodes[m] == inode) continue;
    current_pt = pnodes[m];
    done = false;
    do {
      pcount = nodes[current_pt].path_hop;
      winner.clear();
      for(l=0; l<nodes[current_pt].neighbours.size(); ++l) {
        in1 = nodes[current_pt].neighbours[l];
        temp = nodes[in1].path_hop;
        if (temp >= 0) {
          if (temp < pcount) {
            winner.clear();
            winner.push_back(in1);
          }
          else if (temp == pcount) {
            winner.push_back(in1);
          }
        }
      }
#ifdef DEBUG
      assert(winner.size() > 0);
#endif
      next_node = winner[irandom(winner.size())];
      if (nodes[next_node].locale == 6) {
        for(l=0; l<pathback.size(); ++l) {
          if (nodes[pathback[l]].locale == 2) nodes[pathback[l]].locale = 6;
        }
        done = true;
      }
      else if (next_node == inode) {
        for(l=0; l<pathback.size(); ++l) {
          if (nodes[pathback[l]].locale == 2) nodes[pathback[l]].locale = 6;
        }
        done = true;
      }
      else {
        found = false;
        for(l=0; l<pathback.size(); ++l) {
          if (pathback[l] == next_node) {
            found = true;
            break;
          }
        }
        if (!found) {
          pathback.push_back(next_node);
        }
        current_pt = next_node;
      }
    } while(!done);
  }

  // A final check here to make sure that the molecule, if it consisted of only
  // the paths selected, really is connected:
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        in2 = nodes[in1].locale;
        if (in2 == 5 || in2 == 6) vertices.push_back(in1);
      }
    }
  }

  for(l=0; l<4*vertices.size(); ++l) {
    bonds.push_back(-1);
  }
  for(m=0; m<vertices.size(); ++m) {
    nbonds = 0;
    for(l=0; l<nodes[vertices[m]].neighbours.size(); ++l) {
      in1 = nodes[vertices[m]].neighbours[l];
      in2 = nodes[in1].locale;
      if (in2 == 5 || in2 == 6) {
        bonds[4*m+nbonds] = get_index(in1,vertices);
        nbonds++;
      }
    }
    if (nbonds == 0) {
#ifdef VERBOSE
      std::cout << "Problem: isolated atom at " << vertices[i] << " with locale " << nodes[vertices[i]].locale << std::endl;
#endif
      std::exit(1);
    }
  }
  // A sanity check...
  if (!g_connected(bonds)) {
#ifdef VERBOSE
    std::cout << "Error raised in path routine..." << std::endl;
#endif
    std::exit(1);
  }
  return true;
}

void Grid::clear()
{
  int i,j,k,in1;
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        nodes[in1].locale = 0;
        nodes[in1].atomic_number = 0;
      }
    }
  }
}

void Grid::restore(int q)
{
  int i,j,k,in1;
  unsigned int l;
  state s;

  // First, initialize the whole grid back to zero...
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        nodes[in1].locale = 0;
        nodes[in1].atomic_number = 0;
      }
    }
  }
  // Now refill it from the vector created by save_state:
  for(l=0; l<backup[q].size(); ++l) {
    s = backup[q][l];
    in1 = s.node;
    nodes[in1].locale = s.locale;
    nodes[in1].atomic_number = s.atom;
  }
}

void Grid::save_state(int q)
{
  int i,j,k,in1;
  state s;

  backup[q].clear();
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].locale > 0) {
          // Save this node's info
          s.node = in1;
          s.locale = nodes[in1].locale;
          s.atom = nodes[in1].atomic_number;
          backup[q].push_back(s);
        }
      }
    }
  }
}

bool Grid::create_scaffold()
{
  bool test;

  fill_interior();

  test = initial_deletion(0.5,100);
  if (!test) return false;

  test = path_selection(1);
  if (!test) return false;

  test = secondary_deletion(2,0,4,100);
  if (!test) return false;

  test = rationalize(0.4,1,6);
  if (!test) return false;

  add_hydrogens();

  return true;
}

void Grid::write_scaffold(Molecule* output) const
{
  int i,j,k,cc,in1,in2,na = 0;
  unsigned int l;
  double x[3];
  std::vector<int> carbon;

  // First add the atoms...
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number > 0) {
          x[0] = nodes[in1].x;
          x[1] = nodes[in1].y;
          x[2] = nodes[in1].z;
          output->add_atom(nodes[in1].atomic_number,x,nodes[in1].locale);
          nodes[in1].atom_index = na; na++;
        }
      }
    }
  }
  // Now do the bonds...
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number <= 0) continue;
        if (nodes[in1].atomic_number == 1) {
          // Hydrogen - among your four neighbours, which contains
          // a heavy atom?
          cc = 0;
          for(l=0; l<nodes[in1].neighbours.size(); ++l) {
            in2 = nodes[in1].neighbours[l];
            if (nodes[in2].atomic_number == 6) {
              cc = in2;
              break;
            }
          }
          output->add_bond(nodes[in1].atom_index,nodes[cc].atom_index,1);
        }
        else {
          // Heavy atom, probably carbon, should be bonded to all
          // of its non-empty neighbours
          carbon.clear();
          for(l=0; l<nodes[in1].neighbours.size(); ++l) {
            in2 = nodes[in1].neighbours[l];
            if (nodes[in2].atomic_number > 0) carbon.push_back(in2);
          }
          for(l=0; l<carbon.size(); ++l) {
            output->add_bond(nodes[in1].atom_index,nodes[carbon[l]].atom_index,1);
          }
        }
      }
    }
  }
}






