#include <algorithm>
#include <bitset>

#include "cut.hpp"

using namespace std;

bool Dominate(Cut const &a, Cut const &b) {
  if(a.leaves.size() > b.leaves.size()) {
    return false;
  }
  if((a.signature & b.signature) != a.signature) {
    return false;
  }
  if (a.leaves.size() == b.leaves.size()) {
    return equal(a.leaves.begin(), a.leaves.end(), b.leaves.begin());
  }
  return includes(b.leaves.begin(), b.leaves.end(), a.leaves.begin(), a.leaves.end());
}

void CutEnumeration(aigman const &aig, vector<vector<Cut> > &cuts, unsigned cutsize) {
  cuts.resize(aig.nObjs);
  for(int i = 0; i < aig.nPis; i++) {
    cuts[i + 1].emplace_back(i + 1);
  }
  for(int i = aig.nPis + 1; i < aig.nObjs; i++) {
    int i0 = aig.vObjs[i + i] >> 1;
    int i1 = aig.vObjs[i + i + 1] >> 1;
    for(auto &cut0: cuts[i0]) {
      for(auto &cut1: cuts[i1]) {
        Cut new_cut;
        auto &leaves = new_cut.leaves;
        auto &signature = new_cut.signature;
        signature = cut0.signature | cut1.signature;
        // merge
        if(cut0.leaves.size() + cut1.leaves.size() > cutsize) {
          bitset<64> bs(signature);
          if(bs.count() > cutsize) {
            continue;
          }
        }
        leaves.resize(cut0.leaves.size() + cut1.leaves.size());
        leaves.resize(std::set_union(cut0.leaves.begin(), cut0.leaves.end(), cut1.leaves.begin(), cut1.leaves.end(), leaves.begin()) - leaves.begin());
        if(leaves.size() > cutsize) {
          continue;
        }
        // skip if dominated
        bool dominated = false;
        for(auto &cut: cuts[i]) {
          if(Dominate(cut, new_cut)) {
            dominated = true;
            break;
          }
        }
        if(dominated) {
          continue;
        }
        // insert
        cuts[i].resize(std::stable_partition(cuts[i].begin(), cuts[i].end(), [&](Cut const &cut) { return !Dominate(new_cut, cut); }) - cuts[i].begin());
        cuts[i].push_back(new_cut);
      }
    }
    // unit cut
    cuts[i].emplace_back(i);
  }
}
