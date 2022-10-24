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

void rem(aigman &aig, vector<int> &gates) {
  aig.vValues.clear();
  aig.vValues.resize(aig.nObjs);
  for(int i: gates) {
    for(int ii = i + i; ii <= i + i + 1; ii++) {
      aig.vValues[aig.vObjs[ii] >> 1]++;
    }
  }
  for(auto it = gates.begin(); it != gates.end();) {
    if((int)aig.vvFanouts[*it].size() == aig.vValues[*it]) {
      it = gates.erase(it);
      continue;
    }
    it++;
  }
}

int main(int argc, char **argv) {
  if(argc < 2) {
    return 1;
  }
  string aigname = argv[1];
  aigman aig(aigname);
  aig.supportfanouts();
  string delim;
  bool fSynthesized = true;
  while(fSynthesized) {
    fSynthesized = false;
    // cut enumeration
    vector<vector<Cut> > cuts;
    CutEnumeration(aig, cuts);
    /*
    for(int i = 0; i < aig.nObjs; i++) {
      for(auto &cut: cuts[i]) {
        cout << i << " : ";
        delim = "";
        for(int j: cut.leaves) {
          cout << delim << j;
          delim = ", ";
        }
        cout << endl;
      }
    }
    */
    // cut to gates
    map<vector<int>, vector<int> > m;
    for(int i = 0; i < aig.nObjs; i++) {
      for(auto &cut: cuts[i]) {
        if(cut.leaves.size() != 1) {
          m[cut.leaves].push_back(i);
        }
      }
    }
    for(auto it = m.begin(); it != m.end(); it++) {
      for(auto it2 = m.begin(); it2 != m.end(); it2++) {
        if(it == it2) {
          continue;
        }
        if(includes(it->first.begin(), it->first.end(), it2->first.begin(), it2->first.end())) {
          vector<int> dest(it->second.size() + it2->second.size());
          auto it3 = set_union(it->second.begin(), it->second.end(), it2->second.begin(), it2->second.end(), dest.begin());
          dest.resize(it3 - dest.begin());
          it->second = dest;
        }
      }
    }
    // for each cut
    for(auto p: m) {
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
      // remove internal gates
      auto gates = p.second;
      int nGates = gates.size();
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
      // switch to sub-cut mode
      if(nGates > 7 || aig.reach(p.second, p.first)) {
        for(int i: p.second) {
          cout << "\tFanin cone of " << i << endl;
          if(aig.reach(vector<int>{i}, p.first)) {
            cout << "\tSkipped due to potential loops" << endl;
            continue;
          }
          for(auto cut: cuts[i]) {
            if(cut.leaves == p.first || cut.leaves.size() == 1) {
              continue;
            }
            cout << "\t\t";
            delim = "";
            for(int j: cut.leaves) {
              cout << delim << j;
              delim = ", ";
            }
            cout << " : ";
            delim = "";
            for(int j: m[cut.leaves]) {
              cout << delim << j;
              delim = ", ";
            }
            cout << endl;
            nGates = m[cut.leaves].size();
            if(nGates > 7) {
              cout << "\t\tSkipped because the cone has more than 7 gates" << endl;
              continue;
            }
            if(!includes(gates.begin(), gates.end(), m[cut.leaves].begin(), m[cut.leaves].end())) {
              cout << "\t\tSkipped because the cone is not included" << endl;
              continue;
            }
            auto outputs = m[cut.leaves];
            rem(aig, outputs);
            cout << "\t\t";
            delim = "";
            for(int j: cut.leaves) {
              cout << delim << j;
              delim = ", ";
            }
            cout << " : ";
            delim = "";
            for(int j: outputs) {
              cout << delim << j;
              delim = ", ";
            }
            cout << endl;
            if(aig.reach(outputs, p.first)) {
              cout << "\t\tSkipped due to potential loops" << endl;
              continue;
            }
            vector<int> gates_(gates.size());
            gates_.resize(set_difference(gates.begin(), gates.end(), m[cut.leaves].begin(), m[cut.leaves].end(), gates_.begin()) - gates_.begin());
            /*
            cout << "\t\tGates outside of the cone : ";
            delim = "";
            for(int j: gates_) {
              cout << delim << j;
              delim = ", ";
            }
            cout << endl;
            */
            for(auto it = gates_.begin(); it != gates_.end();) {
              if(aig.reach(outputs, vector<int>{*it})) {
                it = gates_.erase(it);
                continue;
              }
              it++;
            }
            cout << "\t\tExtra inputs : ";
            delim = "";
            for(int j: gates_) {
              cout << delim << j;
              delim = ", ";
            }
            cout << endl;
            vector<vector<bool> > br;
            GetBooleanRelation(aig, p.first, outputs, br);
            /*
            for(int k = 0; k < (int)br.size(); k++) {
              cout << "\t\t" << k << " : ";
              delim = "";
              for(int j = 0; j < (int)br[k].size(); j++) {
                cout << delim << br[k][j];
                delim = ", ";
              }
              cout << endl;
            }
            */
            vector<vector<bool> > sim;
            GetSim(aig, p.first, gates_, sim);
            /*
            for(int k = 0; k < (int)sim.size(); k++) {
              cout << "\t\t" << k << " : ";
              for(int j = 0; j < (int)sim[k].size(); j++) {
                cout << sim[k][j];
              }
              cout << endl;
            }
            */
            cout << "\t\tOptimize a cut of " << nGates << " gates" << endl;
            ExMan<KissatSolver> exman(br, &sim);
            aigman *aig2;
            if((aig2 = exman.ExSynth(nGates))) {
              std::cout << "\t\tSynthesized a cut of " << aig2->nGates << " gates" << std::endl;
              fSynthesized = true;
              std::vector<int> outputs_shift;
              for(int j: outputs) {
                outputs_shift.push_back(j << 1);
              }
              int a = aig.nGates;
              gates_.insert(gates_.begin(), p.first.begin(), p.first.end());
              aig.import(aig2, gates_, outputs_shift);
              /*
              aig.write("y.aig");
              string cmd = "abc -q \"read y.aig; print_stats; cec " + aigname + "\"";
              int r = system(cmd.c_str());
              */
              assert(a - aig.nGates >= nGates - aig2->nGates);
              delete aig2;
              aig.renumber();
              break;
            }
          }
          if(fSynthesized) {
            break;
          }
        }
        if(fSynthesized) {
          break;
        }
        continue;
      }
      // get relation
      vector<vector<bool> > br;
      GetBooleanRelation(aig, p.first, p.second, br);
      /*
      for(int i = 0; i < (int)br.size(); i++) {
        cout << i << " : ";
        delim = "";
        for(int j = 0; j < (int)br[i].size(); j++) {
          cout << delim << br[i][j];
          delim = ", ";
        }
        cout << endl;
      }
      */
      // synthesis
      cout << "Optimize a cut of " << nGates << " gates" << endl;
      ExMan<KissatSolver> exman(br);
      aigman *aig2;
      /*
      aig2 = exman.Synth(nGates);
      assert(aig2);
      delete aig2;
      */
      if((aig2 = exman.ExSynth(nGates))) {
        std::cout << "Synthesized a cut of " << aig2->nGates << " gates" << std::endl;
        fSynthesized = true;
        std::vector<int> outputs_shift;
        for(int i: p.second) {
          outputs_shift.push_back(i << 1);
        }
        int a = aig.nGates;
        aig.import(aig2, p.first, outputs_shift);
        /*
        for(int i = 0; i < aig.nObjs; i++) {
          if(aig.vDeads[i]) {
            cout << i << " is dead " << endl;
          }
        }
        */
        /*
        aig.write("y.aig");
        string cmd = "abc -q \"read y.aig; print_stats; cec " + aigname + "\"";
        int r = system(cmd.c_str());
        */
        assert(a - aig.nGates >= nGates - aig2->nGates);
        delete aig2;
        aig.renumber();
        break;
      }
    }
  }
  aig.write("z.aig");
  return 0;
}
