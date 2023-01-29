#include <fstream>
#include <sstream>
#include <cassert>

#include "ioutil.hpp"
#include "util.hpp"
#include "rel.hpp"

using namespace std;

void ReadBooleanRelation(string fname, vector<vector<bool> > &br, vector<vector<bool> > *&sim, bool fVerbose) {
  ifstream f(fname);
  string line;
  getline(f, line);
  stringstream ss(line);
  int nInputs, nDivisors, nOutputs, nPats;
  ss >> nInputs >> nDivisors >> nOutputs >> nPats;
  assert(nInputs <= 30);
  assert(nOutputs <= 30);
  int nInCombs = 1 << nInputs;
  int nOutCombs = 1 << nOutputs;
  // inputs
  vector<int> pats(nPats);
  getline(f, line);
  for(int i = 0; i < nInputs; i++) {
    getline(f, line);
    for(int j = 0; j < nPats; j++) {
      if(line[j] == '1') {
        pats[j] |= 1 << i;
      }
    }
  }
  if(fVerbose) {
    cout << "Input patterns : " << pats << endl;
  }
  // divisors
  if(nDivisors) {
    sim = new vector<vector<bool> >(nInCombs, vector<bool>(nDivisors));
    getline(f, line);
    for(int i = 0; i < nDivisors; i++) {
      getline(f, line);
      for(int j = 0; j < nPats; j++) {
        if(line[j] == '1') {
          (*sim)[pats[j]][i] = true;
        }
      }
    }
    if(fVerbose) {
      cout << "Divisors :" << endl;
      PrintVecWithIndex(*sim);
    }
  }
  // outputs
  br.resize(nInCombs, vector<bool>(nOutCombs, true));
  getline(f, line);
  for(int i = 0; i < nOutCombs; i++) {
    getline(f, line);
    for(int j = 0; j < nPats; j++) {
      if(line[j] == '0') {
        br[pats[j]][i] = false;
      }
    }
  }
  if(fVerbose) {
    cout << "Boolean relation :" << endl;
    PrintVecWithIndex(br);
  }
}

void WriteBooleanRelation(string fname, vector<vector<bool> > const &br, vector<vector<bool> > const *sim) {
  ofstream f(fname);
  int nInPats = br.size();
  int nOutPats = br[0].size();
  int nInputs = clog2(nInPats);
  int nOutputs = clog2(nOutPats);
  int nDivisors = 0;
  if(sim) {
    nDivisors = (*sim)[0].size();
  }
  f << nInputs << " " << nDivisors << " " << nOutputs << " " << nInPats << endl << endl;
  // inputs
  for(int i = 0; i < nInputs; i++) {
    int dx = 1 << i;
    int nx = nInPats / dx;
    for(int j = 0; j < nx; j++) {
      for(int k = 0; k < dx; k++) {
        f << (j & 1);
      }
    }
    f << endl;
  }
  f << endl;
  // divisors
  if(nDivisors) {
    for(int i = 0; i < nDivisors; i++) {
      for(int j = 0; j < nInPats; j++) {
        f << (*sim)[j][i];
      }
      f << endl;
    }
    f << endl;
  }
  // outputs
  for(int j = 0; j < nOutPats; j++) {
    for(int i = 0; i < nInPats; i++) {
      f << br[i][j];
    }
    f << endl;
  }
}
