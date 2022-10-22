#include <iostream>

#include <map>

#include <fstream>
#include <bitset>
#include <algorithm>

#include <cassert>

#include <aig.hpp>

#include "cut.hpp"
#include "sim.hpp"
#include "synth.hpp"

using namespace std;

void rem(aigman &aig, vector<int> &fos) {
  for(auto it = fos.begin(); it != fos.end();) {
    bool fPo = false;
    for(int i: aig.vPos) {
      if(*it == (i >> 1)) {
        fPo = true;
        break;
      }
    }
    if(fPo) {
      it++;
      continue;
    }
    aig.vValues.clear();
    aig.vValues.resize(aig.nObjs);
    for(int i: aig.vPos) {
      aig.vValues[i >> 1] = 2;
    }
    for(int i: fos) {
      aig.vValues[i] = 1;
    }
    bool fReach = false;
    for(int i: aig.vvFanouts[*it]) {
      if(i > 0 && aig.reach_rec(i)) {
        fReach = true;
        break;
      }
    }
    if(fReach) {
      it++;
      continue;
    }
    it = fos.erase(it);
  }
}

int CountGates(aigman &aig, vector<int> const &inputs, vector<int> const &outputs) {
  aig.vValues.clear();
  aig.vValues.resize(aig.nObjs);
  for(int i: inputs) {
    aig.vValues[i] = 1;
  }
  vector<int> gates;
  for(int i: outputs) {
    aig.getgates_rec(gates, i);
  }
  return gates.size();
}

int main(int argc, char **argv) {
  if(argc < 2) {
    return 1;
  }
  string aigname = argv[1];
  aigman aig;
  aig.read(aigname);
  aig.supportfanouts();
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

  map<vector<int>, vector<int> > m;
  for(int i = 0; i < aig.nObjs; i++) {
    for(auto &cut: cuts[i]) {
      if(cut.leaves.size() != 1) {
        m[cut.leaves].push_back(i);
      }
    }
  }

  for(auto p: m) {
    string delim;
    for(int j: p.first) {
      cout << delim << j;
      delim = ", ";
    }
    cout << " : ";
    delim = "";
    for(int j: p.second) {
      cout << delim << j;
      delim = ", ";
    }
    cout << endl;

    rem(aig, p.second);

    delim = "";
    for(int j: p.first) {
      cout << delim << j;
      delim = ", ";
    }
    cout << " : ";
    delim = "";
    for(int j: p.second) {
      cout << delim << j;
      delim = ", ";
    }
    cout << endl;

    int nGates = CountGates(aig, p.first, p.second);
    cout << "Gates : " << nGates << endl;


    if(nGates > 10) {
      continue;
    }

    vector<vector<bool> > br;
    GetBooleanRelation(aig, p.first, p.second, br);
    // for(int i = 0; i < (int)br.size(); i++) {
    //   cout << i << " : ";
    //   delim = "";
    //   for(int j = 0; j < (int)br[i].size(); j++) {
    //     cout << delim << br[i][j];
    //     delim = ", ";
    //   }
    //   cout << endl;
    // }

    ExMan exman(br);
    aigman *aig2;
    assert(aig2 = exman.Solve(nGates));
    delete aig2;
    if((aig2 = exman.Solve(nGates - 1))) {
      std::vector<int> outputs_shift;
      for(int i: p.second) {
        outputs_shift.push_back(i << 1);
      }
      aig.import(aig2, p.first, outputs_shift);
      delete aig2;
      break;
    }
    continue;

    // ofstream f("tmp.pla");
    // f << ".i " << p.first.size() << endl;
    // f << ".o " << p.second.size() << endl;
    // for(int i = 0; i < (int)br.size(); i++) {
    //   string fival = bitset<32>(i).to_string();
    //   reverse(fival.begin(), fival.end());
    //   fival.resize(p.first.size());
    //   f << fival << " ";
    //   for(int j = 0; j < (int)br[i].size(); j++) {
    //     if(br[i][j]) {
    //       string foval = bitset<32>(j).to_string();
    //       reverse(foval.begin(), foval.end());
    //       foval.resize(p.second.size());
    //       f << foval << endl;
    //       break;
    //     }
    //   }
    // }
    // f << ".e" << endl;
    // string cmd = "abc -q \"read tmp.pla; fx; strash; write_aiger tmp.aig\"";
    // system(cmd.c_str());
    // aigman aig2;
    // aig2.read("tmp.aig");
    // std::vector<int> outputs_shift;
    // for(int i: p.second) {
    //   outputs_shift.push_back(i << 1);
    // }
    // aig.import(&aig2, p.first, outputs_shift);
    // break;
  }

  aig.write("z.aig");

  return 0;
}
