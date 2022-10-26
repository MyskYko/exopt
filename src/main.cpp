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

template <typename T>
ostream &operator<<(ostream &os, vector<T> const &v) {
  string delim;
  for(T const &i: v) {
    os << delim << i;
    delim = ", ";
  }
  return os;
}

ostream &operator<<(ostream &os, Cut const &cut) {
  os << "{" << cut.leaves << "}";
  return os;
}

template <typename T>
void PrintVecWithIndex(vector<T> const &v, string prefix = "") {
  for(int i = 0; i < (int)v.size(); i++) {
    cout << prefix << i << " : " << v[i] << endl;
  }
}

bool Synthesize(aigman &aig, ExMan<KissatSolver> &exman, int nGates, vector<int> const & inputs, vector<int> const & outputs, string prefix = "") {
  cout << prefix << "Synthesizing with less than " << nGates << " gates" << endl;
  aigman *aig2;
  /*
  aig2 = exman.Synth(nGates);
  assert(aig2);
  delete aig2;
  */
  if((aig2 = exman.ExSynth(nGates))) {
    cout << prefix << "Synthesized with " << aig2->nGates << " gates" << endl;
    vector<int> outputs_shift;
    for(int i: outputs) {
      outputs_shift.push_back(i << 1);
    }
    int nGatesAll = aig.nGates;
    aig.import(aig2, inputs, outputs_shift);
    cout << prefix << "Replaced gates : ";
    string delim;
    for(int i = 0; i < aig.nObjs; i++) {
      if(aig.vDeads[i]) {
        cout << delim << i;
        delim = ", ";
      }
    }
    cout << endl;
    /*
    aig.write("y.aig");
    string cmd = "abc -q \"read y.aig; print_stats; cec " + aigname + "\"";
    int r = system(cmd.c_str());
    */
    assert(nGatesAll - aig.nGates >= nGates - aig2->nGates);
    delete aig2;
    aig.renumber();
    return true;
  }
  cout << prefix << "* Synthesis failed" << endl;
  return false;
}

int main(int argc, char **argv) {
  if(argc < 2) {
    return 1;
  }
  int cutsize = 6;
  int windowsize = 6;
  string aigname = argv[1];
  aigman aig(aigname);
  aig.supportfanouts();
  bool fSynthesized = true;
  while(fSynthesized) {
    fSynthesized = false;
    // cut enumeration
    vector<vector<Cut> > cuts;
    CutEnumeration(aig, cuts, cutsize);
    //PrintVecWithIndex(cuts);
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
    vector<tuple<vector<int>, vector<int>, vector<int> > > vWindows;
    for(auto it = mGates.begin(); it != mGates.end();) {
      if((int)it->second.size() > windowsize) {
        it++;
        continue;
      }
      aig.vValues.clear();
      aig.vValues.resize(aig.nObjs);
      for(int i: it->second) {
        for(int ii = i + i; ii <= i + i + 1; ii++) {
          aig.vValues[aig.vObjs[ii] >> 1]++;
        }
      }
      auto outputs = it->second;
      for(auto it2 = outputs.begin(); it2 != outputs.end();) {
        if((int)aig.vvFanouts[*it2].size() == aig.vValues[*it2]) {
          it2 = outputs.erase(it2);
          continue;
        }
        it2++;
      }
      if(aig.reach(outputs, it->first)) {
        it++;
        continue;
      }
      vWindows.push_back(make_tuple(it->first, it->second, outputs));
      it = mGates.erase(it);
    }
    // for each window
    for(auto const &p: vWindows) {
      auto const &inputs = get<0>(p);
      auto const &gates = get<1>(p);
      auto const &outputs = get<2>(p);
      int nGates = gates.size();
      cout << "Inputs : " << inputs << endl;
      cout << "Gates : " << gates << endl;
      cout << "Outputs : " << outputs << endl;
      // get relation
      vector<vector<bool> > br;
      GetBooleanRelation(aig, inputs, outputs, br);
      //PrintVecWithIndex(br);
      // synthesis
      ExMan<KissatSolver> exman(br);
      fSynthesized = Synthesize(aig, exman, nGates, inputs, outputs);
      if(fSynthesized) {
        break;
      }
      continue;
    }
    if(fSynthesized) {
      continue;
    }
    // for each remaining cut
    for(auto const &p: mGates) {
      auto const &inputs = get<0>(p);
      auto const &gates = get<1>(p);
      int nGates = gates.size();
      cout << "Inputs : " << inputs << endl;
      cout << "Gates : " << gates << endl;
      // for each window
      for(auto const&q: vWindows) {
        auto const &inputs2 = get<0>(q);
        auto const &gates2 = get<1>(q);
        auto const &outputs2 = get<2>(q);
        int nGates2 = gates2.size();
        if(!includes(gates.begin(), gates.end(), gates2.begin(), gates2.end())) {
          continue;
        }
        cout << "\t\tInputs : " << inputs2 << endl;
        cout << "\t\tGates : " <<  gates2 << endl;
        cout << "\t\tOutputs : " <<  outputs2 << endl;
        if(aig.reach(outputs2, inputs)) {
          cout << "\t\t* Skipped due to potential loops" << endl;
          continue;
        }
        vector<int> gates_(nGates);
        gates_.resize(set_difference(gates.begin(), gates.end(), gates2.begin(), gates2.end(), gates_.begin()) - gates_.begin());
        cout << "\t\tOutside gates : " << gates_ << endl;
        for(auto it = gates_.begin(); it != gates_.end();) {
          if(aig.reach(outputs2, vector<int>{*it})) {
            it = gates_.erase(it);
            continue;
          }
          it++;
        }
        cout << "\t\tExtra inputs : " << gates_ << endl;
        vector<vector<bool> > br;
        GetBooleanRelation(aig, inputs, outputs2, br);
        //PrintVecWithIndex(br, "\t\t");
        vector<vector<bool> > sim;
        GetSim(aig, inputs, gates_, sim);
        //PrintVecWithIndex(sim, "\t\t");
        ExMan<KissatSolver> exman(br, &sim);
        gates_.insert(gates_.begin(), inputs.begin(), inputs.end());
        fSynthesized = Synthesize(aig, exman, nGates2, gates_, outputs2, "\t\t");
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
