#include <iostream>

#include <aig.hpp>

#include "cut.hpp"

using namespace std;

int main(int argc, char **argv) {
  if(argc < 2) {
    return 1;
  }
  string aigname = argv[1];
  aigman aig;
  aig.read(aigname);
  vector<vector<Cut> > cuts;
  CutEnumeration(aig, cuts);
  for(int i = 0; i < aig.nObjs; i++) {
    for(auto &cut: cuts[i]) {
      cout << i << " : ";
      string delim;
      for(int j: cut.leaves) {
        cout << delim << j;
        delim = ", ";
      }
      cout << endl;
    }
  }

  return 0;
}
