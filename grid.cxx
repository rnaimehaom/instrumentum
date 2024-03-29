#include "grid.h"

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
  if (total > 0) delete[] nodes;
  delete RND;
}

Grid::Grid(int i,int j,int k,unsigned long s,double lambda,int np)
{
  assert(i >= 1 && j >= 1 && k >= 1);
  assert(lambda > std::numeric_limits<double>::epsilon());
  D1 = i;
  D2 = j;
  D3 = k;
  bond_length = lambda;

  allocate(np);
  RND->initialize_generator(s);
  initialize();
}

void Grid::allocate(int np)
{
  for(int i=0; i<np; ++i) {
    pharma_nodes.push_back(-1);
  }
  total = (2*D3 + 1)*(2*D2 + 1)*(2*D1 + 1);
  nodes = new Node[total];
  RND = new Random;
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
        nodes[in1].atomic_number = Atom_Type::empty;
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

void Grid::blank_pharmacophore(double radius)
{
  int i,j,k,in1,candidate;
  unsigned int l,kt = 0;
  bool found,done;
  double delta;
  std::vector<int> vx,frontier_nodes;
  const double r2 = radius*radius;

  for(i=0; i<total; ++i) {
    delta = nodes[i].x*nodes[i].x + nodes[i].y*nodes[i].y + nodes[i].z*nodes[i].z;
    if (delta < r2) nodes[i].locale = 2;
  }

  // Now that the interior nodes have been marked, we can
  // mark frontier nodes, namely locale=0 nodes that are
  // neighbours of a locale=2 node.
  for(i=0; i<total; ++i) {
    if (nodes[i].locale != 0) continue;
    done = false;
    for(l=0; l<nodes[i].neighbours.size(); ++l) {
      in1 = nodes[i].neighbours[l];
      if (nodes[in1].locale == 2) {
        nodes[i].locale = 1;
        frontier_nodes.push_back(i);
        done = true;
        break;
      }
      if (done) break;
    }    
  }

  // Now choose two or three of these frontier nodes to be
  // "must have" nodes, to aid in constructing connected
  // molecules
  while(kt < pharma_nodes.size()) {
    in1 = RND->irandom(frontier_nodes.size());
    candidate = frontier_nodes[in1];
    found = false;
    for(l=0; l<kt; ++l) {
      vx = nodes[pharma_nodes[l]].neighbours;
      if (std::count(vx.begin(),vx.end(),candidate) > 0) {
        found = true;
        break;
      }
    }
    if (!found) {
      pharma_nodes[kt] = candidate;
      nodes[candidate].locale = 5;
      kt++;
    }
  }
  // Lastly, set the values for rs1, rs2 and rs3
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].locale == 1 || nodes[in1].locale == 2) {
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
  int i,j,k,in1,ccount = 0;

  for(i=-rs1; i<=rs1; ++i) {
    for(j=-rs2; j<=rs2; ++j) {
      for(k=-rs3; k<=rs3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].locale == 2) {
          nodes[in1].atomic_number = Atom_Type::carbon;
          ccount++;
        }
      }
    }
  }
#ifdef VERBOSE
  std::cout << "Created " << ccount << " carbon nodes" << std::endl;
#endif
}

bool Grid::connect_pharmacophores()
{
  // To do this, we're first going to paint all of the carbon nodes,
  // and the pharmacophore nodes, as unvisited, and then choose one
  // at random. We then jump outwards to all of its carbon neighbours,
  // and so on, marking each in turn as visited. When we have marched
  // as much as possible, we check to see if any of the pharmacophore
  // nodes are unvisited. If yes, return 0, otherwise, mark up all of
  // unvisited carbon nodes as silver and return 1.
  int i,in1,in2;
  unsigned int l;
  std::vector<int> visited;
  std::set<int> convert,current;
  std::set<int>::const_iterator it;

  for(i=0; i<total; ++i) {
    visited.push_back(-1);
  }

  for(i=0; i<total; ++i) {
    if (nodes[i].atomic_number == Atom_Type::carbon || nodes[i].locale == 5) visited[i] = 0;
  }

  // Now choose one of the pharmacophoric nodes to begin the process...
  visited[pharma_nodes[0]] = 1;
  current.insert(pharma_nodes[0]);
  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      in1 = *it;
      for(l=0; l<nodes[in1].neighbours.size(); ++l) {
        in2 = nodes[in1].neighbours[l];
        if (visited[in2] == 0) convert.insert(in2);
      }
    }
    if (convert.empty()) break;
    for(it=convert.begin(); it!=convert.end(); ++it) {
      visited[*it] = 1;
    }
    current = convert;
    convert.clear();
  } while(true);
  // Last step is to see if all of the pharma_nodes have visited set equal to one,
  // and if so, to set any unvisited carbons to also be silver
  for(l=0; l<pharma_nodes.size(); ++l) {
    if (visited[pharma_nodes[l]] != 1) return false;
  }
  for(i=0; i<total; ++i) {
    if (nodes[i].atomic_number == Atom_Type::carbon && visited[i] == 0) nodes[i].atomic_number = Atom_Type::silver;
  }

  return true;
}

bool Grid::initial_deletion(double percent,int max_attempts)
{
  // This function will make a series of random deletions in the
  // initial volume of carbon atoms, checking to preserve the
  // connectivity at each stage, until either "percent" percentage
  // of the original number of heavy atoms remain, or the number
  // "max_attempts" of deletions has been exhausted.
  int i,j,k,alpha1,alpha2,alpha3,target,in1,kount,nkill,nheavy = 0,iterations = 0;
  double r2,radius,delta,n_debut;

  assert(max_attempts > 0);

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

  if (!connect_pharmacophores()) {
#ifdef VERBOSE
    std::cout << "Bad initial connectivity..." << std::endl;
#endif
    return false;
  }
  do {
    alpha1 = -rs1 + RND->irandom(2*rs1);
    alpha2 = -rs2 + RND->irandom(2*rs2);
    alpha3 = -rs3 + RND->irandom(2*rs3);
    target = index1(alpha1,alpha2,alpha3);
    radius = 1.0 + RND->drandom();
    r2 = radius*radius;
    nkill = 0;
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].locale != 2 || nodes[in1].atomic_number == Atom_Type::empty) continue;
          delta = distance(in1,target);
          if (delta <= r2) {
            nodes[in1].atomic_number = Atom_Type::silver;
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
    if (connect_pharmacophores()) {
      // It's alright to get rid of the silver...
      kount = 0;
      for(i=-rs1; i<=rs1; ++i) {
        for(j=-rs2; j<=rs2; ++j) {
          for(k=-rs3; k<=rs3; ++k) {
            in1 = index1(i,j,k);
            if (nodes[in1].atomic_number == Atom_Type::silver) {
              nodes[in1].atomic_number = Atom_Type::empty;
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
            if (nodes[in1].atomic_number == Atom_Type::silver) nodes[in1].atomic_number = Atom_Type::carbon;
          }
        }
      }
    }
    iterations++;
  } while(iterations < max_attempts);

  return false;
}

bool Grid::rationalize(double percent_methyl,int rings_min,int rings_max)
{
  // This function eliminates flagpole interactions among
  // hydrogen atoms by forcing any node which is a neighbour
  // of two or more carbon atoms to become carbon as well.
  int i,j,k,in1,in2,alpha,nrings,ncarbon,nkill,kount = 0;
  unsigned int l;
  std::vector<int> methyl,convert;

  assert(rings_min >= 0 && rings_max >= 0);

  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number != Atom_Type::empty && nodes[in1].locale == 2) {
          // See if exactly one neighbour is carbon
          ncarbon = 0;
          for(l=0; l<nodes[in1].neighbours.size(); ++l) {
            in2 = nodes[in1].neighbours[l];
            if (nodes[in2].atomic_number == Atom_Type::carbon) ncarbon++;
          }
          if (ncarbon == 1) methyl.push_back(in1);
        }
      }
    }
  }
  nkill = int(percent_methyl*double(methyl.size()));

  if (nkill > 0) {
    do {
      alpha = RND->irandom(methyl.size());
      in1 = methyl[alpha];
      if (nodes[in1].atomic_number != Atom_Type::empty) {
        nodes[in1].atomic_number = Atom_Type::empty;
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
          if (nodes[in1].atomic_number == Atom_Type::empty) {
            // See if more than one neighbour is carbon
            ncarbon = 0;
            for(l=0; l<nodes[in1].neighbours.size(); ++l) {
              in2 = nodes[in1].neighbours[l];
              if (nodes[in2].atomic_number == Atom_Type::carbon) ncarbon++;
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
      nodes[convert[l]].atomic_number = Atom_Type::carbon;
    }
  } while(true);
#ifdef DEBUG
  assert(connect_pharmacophores());
#endif
  // Now the ring counting...
  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number == Atom_Type::silver) nodes[in1].atomic_number = Atom_Type::empty;
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
        if (nodes[in1].atomic_number == Atom_Type::empty) {
          // See if at least one neighbour is carbon
          carbon = false;
          for(l=0; l<nodes[in1].neighbours.size(); ++l) {
            in2 = nodes[in1].neighbours[l];
            if (nodes[in2].atomic_number == Atom_Type::carbon) {
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
    nodes[convert[l]].atomic_number = Atom_Type::hydrogen;
  }
}

int Grid::ring_analysis()
{
  // The first item of business in this routine is to find all
  // of the edges (and so vertices) that are involved in rings.
  // We will do this by removing a bond from the molecule, and
  // seeing if it remains connected.
  int i,j,k,in1,in2,nbonds;
  unsigned int l;
  bool found;
  std::vector<int> vertices,ring_vertices,ring_edges,rings,redges,bonds,wcopy;

  for(i=-D1; i<=D1; ++i) {
    for(j=-D2; j<=D2; ++j) {
      for(k=-D3; k<=D3; ++k) {
        in1 = index1(i,j,k);
        if (nodes[in1].atomic_number != Atom_Type::empty || nodes[in1].locale == 5) vertices.push_back(in1);
      }
    }
  }

  const int nvertices = (signed) vertices.size();
  for(i=0; i<4*nvertices; ++i) {
    bonds.push_back(-1);
  }

  for(i=0; i<nvertices; ++i) {
    nbonds = 0;
    for(l=0; l<nodes[vertices[i]].neighbours.size(); ++l) {
      in2 = nodes[vertices[i]].neighbours[l];
      if (nodes[in2].atomic_number != Atom_Type::empty || nodes[in2].locale == 5) {
        bonds[4*i+nbonds] = get_index(in2,vertices);
        nbonds++;
      }
    }
  }
  wcopy = bonds;

  // A sanity check...
  if (!connected(bonds)) return -1;

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

  if (ring_vertices.size() == 0) return 0;

  for(l=0; l<4*ring_vertices.size(); ++l) {
    redges.push_back(-1);
  }
#ifdef DEBUG
  const int nv = (signed) ring_vertices.size();
#endif  
  for(l=0; l<ring_edges.size(); l+=2) {
    in1 = get_index(ring_edges[l],ring_vertices);
    in2 = get_index(ring_edges[l+1],ring_vertices);
#ifdef DEBUG
    assert(in1 >= 0 && in1 < 4*nv);
    assert(in2 >= 0 && in2 < 4*nv);
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
  ring_perception(redges,rings);
  ring_info.clear();
  for(l=0; l<rings.size(); ++l) {
    ring_info.push_back(vertices[ring_vertices[rings[l]]]);
  }
  return (signed) ring_info.size()/6;
}

bool Grid::secondary_deletion(int nc4,int nc4rings,int nrings,int attempts)
{
  int i,j,k,q,alpha1,alpha2,alpha3,target,in1,in2,cring,ring_count,ring_member;
  int nkill,ncarbon,nbonds,iterations = 0;
  unsigned int l,m;
  std::vector<int> vertices,c4;
  double r2,radius,delta;

  assert(nc4 >= 0);
  assert(nrings >= 0);
  assert(nc4rings >= 0);
  assert(attempts > 0);

  do {
    alpha1 = -rs1 + RND->irandom(2*rs1);
    alpha2 = -rs2 + RND->irandom(2*rs2);
    alpha3 = -rs3 + RND->irandom(2*rs3);
    target = index1(alpha1,alpha2,alpha3);
    radius = 2.0 + 2.5*RND->drandom();
    r2 = radius*radius;
    nkill = 0;
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].locale > 2 || nodes[in1].atomic_number == Atom_Type::empty) continue;
          delta = distance(in1,target);
          if (delta <= r2) {
            nodes[in1].atomic_number = Atom_Type::empty;
            nkill++;
          }
        }
      }
    }
    if (nkill == 0) continue;
    // Let's make sure it worked smoothly...
    vertices.clear();
    for(i=0; i<total; ++i) {
      if (nodes[i].atomic_number != Atom_Type::empty || nodes[i].locale == 5) vertices.push_back(i);
    }
    for(m=0; m<vertices.size(); ++m) {
      nbonds = 0;
      for(l=0; l<nodes[vertices[m]].neighbours.size(); ++l) {
        in2 = nodes[vertices[m]].neighbours[l];
        if (nodes[in2].atomic_number != Atom_Type::empty || nodes[in2].locale == 5) nbonds++;
      }
      if (nbonds == 0) nodes[vertices[m]].atomic_number = Atom_Type::empty;
    }
    assert(connect_pharmacophores());
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].atomic_number == Atom_Type::silver) nodes[in1].atomic_number = Atom_Type::empty;
        }
      }
    }
    // First, see how many carbons with four carbon neighbours we have
    c4.clear();
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].atomic_number == Atom_Type::empty) continue;
          // See if more than neighbour is carbon
          ncarbon = 0;
          for(l=0; l<nodes[in1].neighbours.size(); ++l) {
            in2 = nodes[in1].neighbours[l];
            if (nodes[in2].atomic_number == Atom_Type::carbon) ncarbon++;
          }
          if (ncarbon == 4) c4.push_back(in1);
        }
      }
    }
    if ((signed) c4.size() > nc4) continue;
    // Next, see how many of these four carbon neighbours are in four separate
    // rings
    ring_count = ring_analysis();
    if (ring_count < 0) return false;
    if (ring_count > nrings) continue;
    cring = 0;
    for(l=0; l<c4.size(); ++l) {
      ring_member = 0;
      for(i=0; i<ring_count; ++i) {
        for(q=0; q<6; ++q) {
          if (ring_info[6*i+q] == c4[l]) {
            ring_member++;
            break;
          }
        }
      }
      if (ring_member == 4) cring++;
    }
    if (cring <= nc4rings) return true;
    iterations++;
  } while(iterations < attempts);

  return false;
}

bool Grid::path_selection(bool random)
{
  int i,j,k,inode,temp,pcount,in1,in2,current_pt,next_node,its = 0;
  unsigned int l,m,nbonds;
  bool done;
  std::vector<int> candidate,pathback,winner,vertices,bonds,path_hop;
  std::set<int> convert;
  std::set<int>::const_iterator it;

  for(i=0; i<total; ++i) {
    path_hop.push_back(-1);
  }

  if (random) {
    inode = pharma_nodes[RND->irandom(pharma_nodes.size())];
  }
  else {
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (nodes[in1].atomic_number != Atom_Type::empty) candidate.push_back(in1);
        }
      }
    }
    inode = candidate[RND->irandom(candidate.size())];
  }
  path_hop[inode] = 0;

  do {
    its++;
    convert.clear();
    for(i=-rs1; i<=rs1; ++i) {
      for(j=-rs2; j<=rs2; ++j) {
        for(k=-rs3; k<=rs3; ++k) {
          in1 = index1(i,j,k);
          if (path_hop[in1] < 0) {
            if (nodes[in1].atomic_number != Atom_Type::empty || nodes[in1].locale == 5) {
              for(l=0; l<nodes[in1].neighbours.size(); ++l) {
                in2 = nodes[in1].neighbours[l];
                if (path_hop[in2] >= 0) convert.insert(in1);
              }
            }
          }
        }
      }
    }
    if (convert.empty()) break;
    for(it=convert.begin(); it!=convert.end(); ++it) {
      path_hop[*it] = its;
    }
  } while(true);

  winner.reserve(4);
  for(m=0; m<pharma_nodes.size(); ++m) {
    if (pharma_nodes[m] == inode) continue;
    current_pt = pharma_nodes[m];
    done = false;
    do {
      pcount = path_hop[current_pt];
      winner.clear();
      for(l=0; l<nodes[current_pt].neighbours.size(); ++l) {
        in1 = nodes[current_pt].neighbours[l];
        temp = path_hop[in1];
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
      if (winner.empty()) return false;
      next_node = winner[RND->irandom(winner.size())];
      if (nodes[next_node].locale == 6 || next_node == inode) {
        for(l=0; l<pathback.size(); ++l) {
          if (nodes[pathback[l]].locale == 2) nodes[pathback[l]].locale = 6;
        }
        done = true;
      }
      else {
        if (std::count(pathback.begin(),pathback.end(),next_node) == 0) pathback.push_back(next_node);
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
    if (nbonds == 0) return false;
  }
  // A sanity check...
  if (!connected(bonds)) return false;

  return true;
}

void Grid::write_scaffold(Molecule* output) const
{
  int i,cc,in1,na = 0;
  unsigned int l;
  double xc[3];
  std::vector<int> carbon,atom_index;

  for(i=0; i<total; ++i) {
    atom_index.push_back(-1);
  }

  // First add the atoms...
  for(i=0; i<total; ++i) {
    if (nodes[i].atomic_number == Atom_Type::empty) continue;
    xc[0] = nodes[i].x;
    xc[1] = nodes[i].y;
    xc[2] = nodes[i].z;
    output->add_atom(nodes[i].atomic_number,xc,nodes[i].locale);
    atom_index[i] = na; na++;
  }
  // Now do the bonds...
  for(i=0; i<total; ++i) {
    if (nodes[i].atomic_number == Atom_Type::empty) continue;
    if (nodes[i].atomic_number == Atom_Type::hydrogen) {
      // Hydrogen - among your four neighbours, which contains
      // a heavy atom?
      cc = 0;
      for(l=0; l<nodes[i].neighbours.size(); ++l) {
        in1 = nodes[i].neighbours[l];
        if (nodes[in1].atomic_number == Atom_Type::carbon) {
          cc = in1;
          break;
        }
      }
      output->add_bond(atom_index[i],atom_index[cc],1);
    }
    else {
      // Heavy atom, probably carbon, should be bonded to all
      // of its non-empty neighbours
      carbon.clear();
      for(l=0; l<nodes[i].neighbours.size(); ++l) {
        in1 = nodes[i].neighbours[l];
        if (nodes[in1].atomic_number != Atom_Type::empty) carbon.push_back(in1);
      }
      for(l=0; l<carbon.size(); ++l) {
        output->add_bond(atom_index[i],atom_index[carbon[l]],1);
      }
    }
  }
}






