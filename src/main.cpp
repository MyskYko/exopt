#include <iostream>

#include <map>

#include <fstream>
#include <bitset>
#include <algorithm>
#include <random>

#include <cassert>

#include <argparse/argparse.hpp>

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

bool Synthesize(aigman &aig, SynthMan<KissatSolver> &exman, int nGates, vector<int> const & inputs, vector<int> const & outputs, bool fVerbose, string prefix = "") {
  if(fVerbose) {
    cout << prefix << "Synthesizing with less than " << nGates << " gates" << endl;
  }
  aigman *aig2;
  if((aig2 = exman.ExSynth(nGates))) {
    if(fVerbose) {
      cout << prefix << "Synthesized with " << aig2->nGates << " gates" << endl;
    }
    vector<int> outputs_shift;
    for(int i: outputs) {
      outputs_shift.push_back(i << 1);
    }
    int nGatesAll = aig.nGates;
    aig.import(aig2, inputs, outputs_shift);
    if(fVerbose) {
      cout << prefix << "Replaced gates : ";
      string delim;
      for(int i = 0; i < aig.nObjs; i++) {
        if(aig.vDeads[i]) {
          cout << delim << i;
          delim = ", ";
        }
      }
      cout << endl;
    }
    assert(nGatesAll - aig.nGates >= nGates - aig2->nGates);
    delete aig2;
    aig.renumber();
    return true;
  }
  if(fVerbose) {
    cout << prefix << "* Synthesis failed" << endl;
  }
  return false;
}

template <typename T>
void RemoveIncluded(T &s) {
  for(auto it = s.begin(); it != s.end();) {
    bool fIncluded = false;
    for(auto it2 = s.begin(); it2 != s.end(); it2++){
      if(it == it2) {
        continue;
      }
      if(includes(get<1>(*it2).begin(), get<1>(*it2).end(), get<1>(*it).begin(), get<1>(*it).end())) {
        fIncluded = true;
        break;
      }
    }
    if(fIncluded) {
      it = s.erase(it);
    } else {
      it++;
    }
  }
}

template <typename T>
void RemoveIncluded(T &s, aigman &aig) {
  for(auto it = s.begin(); it != s.end();) {
    bool fIncluded = false;
    for(auto it2 = s.begin(); it2 != s.end(); it2++){
      if(it == it2) {
        continue;
      }
      if(includes(get<1>(*it2).begin(), get<1>(*it2).end(), get<1>(*it).begin(), get<1>(*it).end())) {
        fIncluded = true;
        for(int i: get<2>(*it2)) {
          vector<int> v{i};
          if(!aig.reach(get<2>(*it), v)) {
            fIncluded = false;
            break;
          }
        }
        break;
      }
    }
    if(fIncluded) {
      it = s.erase(it);
    } else {
      it++;
    }
  }
}

bool Run(aigman &aig, vector<tuple<vector<int>, vector<int>, vector<int> > > vWindows, bool fVerbose) {
  RemoveIncluded(vWindows);
  for(auto const &p: vWindows) {
    auto const &inputs = get<0>(p);
    auto const &gates = get<1>(p);
    auto const &outputs = get<2>(p);
    int nGates = gates.size();
    if(fVerbose) {
      cout << "Inputs : " << inputs << endl;
      cout << "Gates : " << gates << endl;
      cout << "Outputs : " << outputs << endl;
    }
    // get relation
    vector<vector<bool> > br;
    GetBooleanRelation(aig, inputs, outputs, br);
    //PrintVecWithIndex(br);
    // synthesis
    SynthMan<KissatSolver> exman(br);
    if(Synthesize(aig, exman, nGates, inputs, outputs, fVerbose)) {
      return true;
    }
  }
  return false;
}

int main(int argc, char **argv) {
  argparse::ArgumentParser ap("exopt");
  ap.add_argument("input");
  ap.add_argument("output");
  ap.add_argument("-k", "--cutsize").default_value(8).scan<'i', int>();
  ap.add_argument("-n", "--windowsize").default_value(6).scan<'i', int>();
  ap.add_argument("-a", "--alldivisors").default_value(false).implicit_value(true);
  ap.add_argument("-r", "--numrounds").default_value(10).scan<'i', int>();
  ap.add_argument("-v", "--verbose").default_value(false).implicit_value(true);
  try {
    ap.parse_args(argc, argv);
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    cerr << ap;
    return 1;
  }
  string aigname = ap.get<string>("input");
  string outname = ap.get<string>("output");
  int cutsize = ap.get<int>("--cutsize");
  int windowsize = ap.get<int>("--windowsize");
  bool fAllDivisors = ap.get<bool>("--alldivisors");
  int numrounds = ap.get<int>("--numrounds");
  bool fVerbose = ap.get<bool>("--verbose");
  mt19937 rg;
  aigman aig_orig(aigname);
  aig_orig.supportfanouts();
  aigman aigout = aig_orig;
  for(int round = 0; round < numrounds; round++) {
  aigman aig = aig_orig;
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
    if(round) {
      shuffle(vWindows.begin(), vWindows.end(), rg);
    }
    // for each window
    fSynthesized = Run(aig, vWindows, fVerbose);
    if(fSynthesized) {
      continue;
    }
    // large cuts
    if(fAllDivisors) {
      mGates.clear();
      vector<int> inputs;
      for(int i = 0; i < aig.nPis; i++) {
        inputs.push_back(i + 1);
      }
      vector<int> gates;
      for(int i = aig.nPis + 1; i < aig.nObjs; i++) {
        gates.push_back(i);
      }
      mGates[inputs] = gates;
    } else {
      RemoveIncluded(mGates);
    }
    vector<pair<vector<int>, vector<int> > > vGates(mGates.begin(), mGates.end());
    if(round) {
      shuffle(vGates.begin(), vGates.end(), rg);
    }
    // for each large cut
    for(auto const &p: vGates) {
      auto const &inputs = get<0>(p);
      auto const &gates = get<1>(p);
      int nGates = gates.size();
      if(fVerbose) {
        cout << "Inputs : " << inputs << endl;
        cout << "Gates : " << gates << endl;
      }
      // for each window
      vector<tuple<vector<int>, vector<int>, vector<int> > > vWindows2;
      for(auto const&q: vWindows) {
        auto const &gates2 = get<1>(q);
        auto const &outputs2 = get<2>(q);
        if(!includes(gates.begin(), gates.end(), gates2.begin(), gates2.end())) {
          continue;
        }
        if(aig.reach(outputs2, inputs)) {
          continue;
        }
        vWindows2.push_back(q);
      }
      RemoveIncluded(vWindows2, aig);
      for(auto const&q: vWindows2) {
        auto const &inputs2 = get<0>(q);
        auto const &gates2 = get<1>(q);
        auto const &outputs2 = get<2>(q);
        int nGates2 = gates2.size();
        if(fVerbose) {
          cout << "\t\tInputs : " << inputs2 << endl;
          cout << "\t\tGates : " <<  gates2 << endl;
          cout << "\t\tOutputs : " <<  outputs2 << endl;
        }
        vector<int> extra(nGates);
        extra.resize(set_difference(gates.begin(), gates.end(), gates2.begin(), gates2.end(), extra.begin()) - extra.begin());
        if(fVerbose) {
          cout << "\t\tOutside gates : " << extra << endl;
        }
        for(auto it = extra.begin(); it != extra.end();) {
          if(aig.reach(outputs2, vector<int>{*it})) {
            it = extra.erase(it);
            continue;
          }
          it++;
        }
        if(fVerbose) {
          cout << "\t\tExtra inputs : " << extra << endl;
        }
        vector<vector<bool> > br;
        GetBooleanRelation(aig, inputs, outputs2, br);
        //PrintVecWithIndex(br, "\t\t");
        vector<vector<bool> > sim;
        GetSim(aig, inputs, extra, sim);
        //PrintVecWithIndex(sim, "\t\t");
        SynthMan<KissatSolver> exman(br, &sim);
        extra.insert(extra.begin(), inputs.begin(), inputs.end());
        fSynthesized = Synthesize(aig, exman, nGates2, extra, outputs2, fVerbose, "\t\t");
        if(fSynthesized) {
          break;
        }
      }
      if(fSynthesized) {
        break;
      }
    }
  }
  cout << aig.nGates << endl;
  if(aig.nGates < aigout.nGates) {
    aigout = aig;
  }
  }
  aigout.write(outname);
  return 0;
}
