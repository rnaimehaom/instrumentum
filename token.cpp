#include "token.h"

int Token::rcount = 0;
int Token::nedges = 0;
int** Token::edges;
int** Token::rinfo;

Token::Token()
{
  visited.reserve(200);
}

Token::~Token()
{

}

Token::Token(int sv)
{
  visited.reserve(200);
  visited.push_back(sv);
  hop_count = 0;
}

void Token::initialize(int sv)
{
  visited.push_back(sv);
  hop_count = 0;
}

Token::Token(const std::vector<int>& sv)
{
  visited.reserve(200);
  visited = sv;

  hop_count = visited.size() - 1;
}

void Token::one_hop()
{
  int i,j,starting_pt;
  std::vector<int> temp;
  bool found,fired = false;

  temp.reserve(200);
  starting_pt = visited[hop_count];
  for(i=0; i<nedges; ++i) {
    if (starting_pt == edges[i][0] || starting_pt == edges[i][1]) {
      found = false;
      for(j=0; j<(signed) visited.size(); ++j) {
        if (edges[i][0] == visited[j]) found = true;
      }
      if (!found) {
        if (!fired) {
          visited.push_back(edges[i][0]);
          hop_count++;
          fired = true;
        }
        else {
          temp = visited;
          temp[hop_count] = edges[i][0];
          Token t = Token(temp);
          t.circuit();
        }
        continue;
      }
      found = false;
      for(j=0; j<(signed) visited.size(); ++j) {
        if (edges[i][1] == visited[j]) found = true;
      }
      if (!found) {
        if (!fired) {
          visited.push_back(edges[i][1]);
          hop_count++;
          fired = true;
        }
        else {
          temp = visited;
          temp[hop_count] = edges[i][1];
          Token t = Token(temp);
          t.circuit();
        }
        continue;
      }
      if ((visited[0] == edges[i][0] || visited[0] == edges[i][1]) && hop_count > 2) {
        visited.push_back(visited[0]);
        hop_count++;
        return;
      }
    }
  }
}

bool Token::exists() const
{
  int i,j;
  bool good;

  for(i=0; i<Token::rcount; ++i) {
    good = false;
    for(j=0; j<6; ++j) {
      if (rinfo[i][j] != visited[j]) {
        good = true;
        break;
      }
    }
    if (!good) return true;
  }
  return false;
}

void Token::clear()
{
  visited.clear();
  hop_count = -1;
}

void Token::circuit()
{
  int i,kount = 0;

  do {
    one_hop();
    kount++;
    if (visited[hop_count] == visited[0]) {
      visited.erase(visited.end()-1,visited.end());
      std::sort(visited.begin(),visited.end());
      if (visited.size() == 6) {
        if (!exists()) {
          for(i=0; i<6; ++i) {
            rinfo[Token::rcount][i] = visited[i];
          }
          Token::rcount++;
        }
      }
      break;
    }
    else if (hop_count >= 7) {
      break;
    }
    if (kount > 20) break;
  } while(true);
}
