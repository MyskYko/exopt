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
    map<vector<int>, vector<int> > mGates;
    for(int i = 0; i < aig.nObjs; i++) {
      for(auto const &cut: cuts[i]) {
        if(cut.leaves.size() != 1) {
          mGates[cut.leaves].push_back(i);
        }
      }
    }
    for(auto it = mGates.begin(); it != mGates.end(); it++) {
      for(auto it2 = mGates.begin(); it2 != mGates.end(); it2++) {
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
    // windows
    map<vector<int>, tuple<vector<int>, vector<int>, bool> > mWindows;
    for(auto const &p: mGates) {
      aig.vValues.clear();
      aig.vValues.resize(aig.nObjs);
      for(int i: p.second) {
        for(int ii = i + i; ii <= i + i + 1; ii++) {
          aig.vValues[aig.vObjs[ii] >> 1]++;
        }
      }
      auto outputs = p.second;
      for(auto it = outputs.begin(); it != outputs.end();) {
        if((int)aig.vvFanouts[*it].size() == aig.vValues[*it]) {
          it = outputs.erase(it);
          continue;
        }
        it++;
      }
      mWindows[p.first] = make_tuple(p.second, outputs, aig.reach(outputs, p.first));
    }
    // for each window
    for(auto const &p: mWindows) {
      auto const &inputs = p.first;
      auto const &gates = get<0>(p.second);
      auto const &outputs = get<1>(p.second);
      const bool fReach = get<2>(p.second);
      int nGates = gates.size();
      delim = "";
      for(int j: inputs) {
        cout << delim << j;
        delim = ", ";
      }
      cout << " : ";
      delim = "";
      for(int j: gates) {
        cout << delim << j;
        delim = ", ";
      }
      cout << endl;
      delim = "";
      for(int j: inputs) {
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
      // single window mode
      if(!fReach && nGates <= 7) {
        // get relation
        vector<vector<bool> > br;
        GetBooleanRelation(aig, inputs, outputs, br);
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
          for(int i: outputs) {
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
        continue;
      }
      // sub-window mode
      for(int i: outputs) {
        cout << "\tFanin cone of " << i << endl;
        if(aig.reach(vector<int>{i}, inputs)) {
          cout << "\tSkipped due to potential loops" << endl;
          continue;
        }
        for(auto const &cut: cuts[i]) {
          if(cut.leaves == inputs || cut.leaves.size() == 1) {
            continue;
          }
          auto const &p2 = mWindows[cut.leaves];
          auto const &gates2 = get<0>(p2);
          auto const &outputs2 = get<1>(p2);
          const bool fReach2 = get<2>(p2);
          int nGates2 = gates2.size();
          cout << "\t\t";
          delim = "";
          for(int j: cut.leaves) {
            cout << delim << j;
            delim = ", ";
          }
          cout << " : ";
          delim = "";
          for(int j: gates2) {
            cout << delim << j;
            delim = ", ";
          }
          cout << endl;
          if(nGates2 > 7) {
            cout << "\t\tSkipped because the cone has more than 7 gates" << endl;
            continue;
          }
          if(!includes(gates.begin(), gates.end(), gates2.begin(), gates2.end())) {
            cout << "\t\tSkipped because the cone is not included" << endl;
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
          for(int j: outputs2) {
            cout << delim << j;
            delim = ", ";
          }
          cout << endl;
          // should be switched
          if(fReach2 || aig.reach(outputs2, inputs)) {
            cout << "\t\tSkipped due to potential loops" << endl;
            continue;
          }
          vector<int> gates_(nGates);
          gates_.resize(set_difference(gates.begin(), gates.end(), gates2.begin(), gates2.end(), gates_.begin()) - gates_.begin());
          cout << "\t\tGates outside of the cone : ";
          delim = "";
          for(int j: gates_) {
            cout << delim << j;
            delim = ", ";
          }
          cout << endl;
          for(auto it = gates_.begin(); it != gates_.end();) {
            if(aig.reach(outputs2, vector<int>{*it})) {
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
          GetBooleanRelation(aig, inputs, outputs2, br);
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
          GetSim(aig, inputs, gates_, sim);
          /*
          for(int k = 0; k < (int)sim.size(); k++) {
            cout << "\t\t" << k << " : ";
            for(int j = 0; j < (int)sim[k].size(); j++) {
              cout << sim[k][j];
            }
            cout << endl;
          }
          */
          cout << "\t\tOptimize a cut of " << nGates2 << " gates" << endl;
          ExMan<KissatSolver> exman(br, &sim);
          aigman *aig2;
          /*
          aig2 = exman.Synth(nGates2);
          assert(aig2);
          delete aig2;
          */
          if((aig2 = exman.ExSynth(nGates2))) {
            std::cout << "\t\tSynthesized a cut of " << aig2->nGates << " gates" << std::endl;
            fSynthesized = true;
            std::vector<int> outputs_shift;
            for(int j: outputs2) {
              outputs_shift.push_back(j << 1);
            }
            int a = aig.nGates;
            gates_.insert(gates_.begin(), inputs.begin(), inputs.end());
            aig.import(aig2, gates_, outputs_shift);
            /*
            aig.write("y.aig");
            string cmd = "abc -q \"read y.aig; print_stats; cec " + aigname + "\"";
            int r = system(cmd.c_str());
            */
            assert(a - aig.nGates >= nGates2 - aig2->nGates);
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
    }
  }
  aig.write("z.aig");
  return 0;
}
