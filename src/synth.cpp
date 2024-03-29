#include <cassert>
#include <algorithm>

#include "util.hpp"
#include "synth.hpp"

using namespace std;

template <class T>
SynthMan<T>::SynthMan(vector<vector<bool> > const &br, vector<vector<bool> > const *sim): br(br), sim(sim) {
  nInputs = clog2(br.size());
  nOutputs = clog2(br[0].size());
  if(sim) {
    assert(nInputs == clog2(sim->size()));
    nExtraInputs = (*sim)[0].size();
  } else {
    nExtraInputs = 0;
  }
}

template <class T>
void SynthMan<T>::GenSels() {
  negs.clear();
  negs.resize(nGates * 3);
  for(int i = 0; i < nGates; i++) {
    negs[i*3] = S->NewVar();
    negs[i*3+1] = S->NewVar();
    negs[i*3+2] = S->NewVar();
    S->AddClause(negs[i*3], negs[i*3+1], negs[i*3+2]);
    S->AddClause(-negs[i*3], negs[i*3+1], -negs[i*3+2]);
    S->AddClause(negs[i*3], -negs[i*3+1], -negs[i*3+2]);
    S->AddClause(-negs[i*3], -negs[i*3+1], negs[i*3+2]);
  }
  sels.clear();
  sels.resize(nGates * 2);
  for(int i = 0; i < nGates * 2; i++) {
    sels[i].resize(nInputs + nExtraInputs + i/2 - 1);
    for(int j = 0; j < nInputs + nExtraInputs + i/2 - 1; j++) {
      sels[i][j] = S->NewVar();
    }
    S->Onehot(sels[i]);
  }
  ponegs.clear();
  ponegs.resize(nOutputs);
  posels.clear();
  posels.resize(nOutputs);
  for(int i = 0; i < nOutputs; i++) {
    ponegs[i] = S->NewVar();
    posels[i].resize(nInputs + nExtraInputs + nGates);
    for(int j = 0; j < nInputs + nExtraInputs + nGates; j++) {
      posels[i][j] = S->NewVar();
    }
    S->AMO(posels[i]);
  }
}

template <class T>
void SynthMan<T>::SortSels() {
  // fanin constraints for each gate
  for(int i = 0; i < nGates; i++) {
    for(int j = 0; j < nInputs + nExtraInputs + i - 2; j++) {
      vector<int> vLits(2);
      vLits[0] = -sels[i + i][j];
      for(int k = j + 1; k < nInputs + nExtraInputs + i - 1; k++) {
        vLits[1] = -sels[i + i + 1][k];
        S->AddClause(vLits);
      }
    }
  }
  // fanin constraints for consecutive gates
  for(int i = 0; i < nGates - 1; i++) {
    for(int j = 0; j < nInputs + nExtraInputs + i - 2; j++) {
      vector<int> vLits(2);
      vLits[0] = -sels[i + i + 2][j];
      for(int k = j + 1; k < nInputs + nExtraInputs + i - 1; k++) {
        vLits[1] = -sels[i + i][k];
        S->AddClause(vLits);
      }
    }
  }
  for(int i = 0; i < nGates - 1; i++) {
    for(int j = 1; j < nInputs + nExtraInputs + i - 1; j++) {
      for(int k = 0; k < j; k++) {
        vector<int> vLits(4);
        vLits[0] = -sels[i + i + 3][k];
        vLits[1] = -sels[i + i][j];
        vLits[2] = -sels[i + i + 2][j];
        for(int l = k + 1; l <= j; l++) {
          vLits[3] = -sels[i + i + 1][l];
          S->AddClause(vLits);
        }
      }
    }
  }
  // gate irredundance constraints
  for(int i = 1; i < nGates; i++) {
    for(int j = i - 1; j >= 0; j--) {
      for(int k = 0; k < nInputs + nExtraInputs + j - 1; k++) {
        vector<int> vLits(3);
        vLits[0] = -sels[i + i][nInputs + nExtraInputs + j - 1];
        vLits[1] = -sels[j + j][k];
        vLits[2] = -sels[i + i + 1][k + 1];
        S->AddClause(vLits);
      }
      for(int k = 0; k < nInputs + nExtraInputs + j - 1; k++) {
        vector<int> vLits(3);
        vLits[0] = -sels[i + i][nInputs + nExtraInputs + j - 1];
        vLits[1] = -sels[j + j + 1][k];
        vLits[2] = -sels[i + i + 1][k];
        S->AddClause(vLits);
      }
    }
  }
}

template <class T>
void SynthMan<T>::GenOne(vector<int> cands, vector<int> const &pos) {
  cands.resize(nInputs + nExtraInputs + nGates);
  for(int i = 0; i < nGates; i++) {
    vector<int> fis(2);
    for(int k = 0; k <= 1; k++) {
      fis[k] = S->NewVar();
      for(int j = 0; j < nInputs + nExtraInputs + i - 1; j++) {
        vector<int> vLits(3);
        vLits[0] = -sels[i + i + k][j];
        for(int cand = 0; cand < 2; cand++) {
          vLits[1] = cand? -cands[j + 1 - k]: cands[j + 1 - k];
          vLits[2] = cand? fis[k]: -fis[k];
          S->AddClause(vLits);
        }
      }
    }
    cands[nInputs + nExtraInputs + i] = S->NewVar();
    for(int k = 0; k < 4; k++) {
      vector<int> vLits(4);
      vLits[0] = (k & 1)? -fis[0]: fis[0];
      vLits[1] = (k >> 1)? -fis[1]: fis[1];
      vLits[2] = k? negs[i*3 + k-1]: S->zero;
      vLits[3] = -cands[nInputs + nExtraInputs + i];
      S->AddClause(vLits);
      if(k) {
        vLits[2] = -vLits[2];
        vLits[3] = -vLits[3];
        S->AddClause(vLits);
      }
    }
  }
  for(int i = 0; i < nOutputs; i++) {
    for(int j = 0; j < nInputs + nExtraInputs + nGates; j++) {
      vector<int> vLits(4);
      vLits[0] = -posels[i][j];
      for(int neg = 0; neg < 2; neg++) {
        for(int cand = 0; cand < 2; cand++) {
          vLits[1] = neg? -ponegs[i]: ponegs[i];
          vLits[2] = cand? -cands[j]: cands[j];
          vLits[3] = (neg ^ cand)? pos[i]: -pos[i];
          S->AddClause(vLits);
        }
      }
    }
    vector<int> vLits(nInputs + nExtraInputs + nGates + 2);
    for(int j = 0; j < nInputs + nExtraInputs + nGates; j++) {
      vLits[2 + j] = posels[i][j];
    }
    for(int neg = 0; neg < 2; neg++) {
      vLits[0] = neg? -ponegs[i]: ponegs[i];
      vLits[1] = neg? pos[i]: -pos[i];
      S->AddClause(vLits);
    }
  }
}

template <class T>
aigman *SynthMan<T>::GetAig() {
  aigman *aig = new aigman(nInputs + nExtraInputs, 0);
  vector<int> cands(nInputs + nExtraInputs + nGates);
  for(int i = 0; i < nInputs + nExtraInputs; i++) {
    cands[i] = (i + 1) << 1;
  }
  for(int i = 0; i < nGates; i++) {
    vector<int> fis(2);
    for(int k = 0; k <= 1; k++) {
      int j = 0;
      for(; j < nInputs + nExtraInputs + i - 1; j++) {
        if(S->Value(sels[i + i + k][j])) {
          break;
        }
      }
      assert(j < nInputs + nExtraInputs + i - 1);
      fis[k] = cands[j + 1 - k];
    }
    if(S->Value(negs[i*3]) && S->Value(negs[i*3 + 1]) && S->Value(negs[i*3 + 2])) {
      cands[nInputs + nExtraInputs + i] = (aig->newgate(fis[0] ^ 1, fis[1] ^ 1) << 1) ^ 1;
    } else {
      fis[0] ^= S->Value(negs[i*3 + 1]);
      fis[1] ^= S->Value(negs[i*3]);
      cands[nInputs + nExtraInputs + i] = aig->newgate(fis[0], fis[1]) << 1;
    }
  }
  aig->nPos = nOutputs;
  aig->vPos.resize(nOutputs);
  for(int i = 0; i < nOutputs; i++) {
    int j = 0;
    for(; j < nInputs + nExtraInputs + nGates; j++) {
      if(S->Value(posels[i][j])) {
        break;
      }
    }
    assert(j <= nInputs + nExtraInputs + nGates);
    int val = (j == nInputs + nExtraInputs + nGates)? 0: cands[j];
    aig->vPos[i] = val ^ (int)S->Value(ponegs[i]);
  }
  return aig;
}

template <class T>
aigman *SynthMan<T>::Synth(int nGates_) {
  nGates = nGates_;
  S = new T;
  GenSels();
  SortSels();
  for(int i = 0; i < (int)br.size(); i++) {
    int nOnes = count(br[i].begin(), br[i].end(), true);
    if(nOnes == (int)br[i].size()) {
      continue;
    }
    vector<int> pis(nInputs);
    for(int k = 0; k < nInputs; k++) {
      pis[k] = (i >> k) & 1? S->one: S->zero;
    }
    vector<int> exins(nExtraInputs);
    for(int k = 0; k < nExtraInputs; k++) {
      exins[k] = (*sim)[i][k]? S->one: S->zero;
    }
    vector<int> pos(nOutputs);
    if(nOnes == 1) {
      for(int j = 0; j < (int)br[i].size(); j++) {
        if(br[i][j]) {
          for(int k = 0; k < nOutputs; k++) {
            pos[k] = (j >> k) & 1? S->one: S->zero;
          }
          break;
        }
      }
    } else {
      for(int k = 0; k < nOutputs; k++) {
        pos[k] = S->NewVar();
      }
      vector<int> tmps;
      for(int j = 0; j < (int)br[i].size(); j++) {
        if(br[i][j]) {
          vector<int> vLits(nOutputs);
          for(int k = 0; k < nOutputs; k++) {
            vLits[k] = (j >> k) & 1? pos[k]: -pos[k];
          }
          tmps.push_back(S->AndN(vLits));
        }
      }
      S->AddClause(tmps);
    }
    pis.insert(pis.end(), exins.begin(), exins.end());
    GenOne(pis, pos);
  }
  aigman *aig = NULL;
  if(S->Solve() == 1) {
    aig = GetAig();
  }
  delete S;
  return aig;
}

template <class T>
aigman *SynthMan<T>::ExSynth(int nGates_) {
  assert(nGates_>= 0);
  aigman *aig = NULL;
  while(--nGates_ >= 0) {
    aigman *aig2 = Synth(nGates_);
    if(!aig2) {
      break;
    }
    if(aig) {
      delete aig;
    }
    aig = aig2;
  }
  return aig;
}

void Enumerate(vector<int> &v, int i, vector<vector<int> > &all) {
  if(i+i >= (int)v.size()) {
    all.push_back(vector<int>(v.begin() + 2, v.end()));
    return;
  }
  v[i+i] = v[i+i-2];
  if(v[i+i] == 0) {
    v[i+i+1] = 0;
    Enumerate(v, i+1, all);
  }
  for(int k = v[i+i-1]; k < v[i+i]; k++) {
    if(!k || (v[v[i+i]+v[i+i]] != k && v[v[i+i]+v[i+i]+1] != k)) {
      v[i+i+1] = k;
      Enumerate(v, i+1, all);
    }
  }
  for(int j = v[i+i-2] + 1; j < i; j++) {
    v[i+i] = j;
    for(int k = 0; k < j; k++) {
      if(!k || (v[v[i+i]+v[i+i]] != k && v[v[i+i]+v[i+i]+1] != k)) {
        v[i+i+1] = k;
        Enumerate(v, i+1, all);
      }
    }
  }
}

template <class T>
void SynthMan<T>::GenSels(vector<int> const &assignment) {
  negs.clear();
  negs.resize(nGates * 2);
  sels.clear();
  sels.resize(nGates * 2);
  for(int i = 0; i < nGates * 2; i++) {
    negs[i] = S->NewVar();
    if(!assignment[i]) {
      sels[i].resize(nInputs + nExtraInputs + (i%2) - 1);
      for(int j = 0; j < nInputs + nExtraInputs + (i%2) - 1; j++) {
        sels[i][j] = S->NewVar();
      }
      S->Onehot(sels[i]);
    }
  }
  ponegs.clear();
  ponegs.resize(nOutputs);
  posels.clear();
  posels.resize(nOutputs);
  for(int i = 0; i < nOutputs; i++) {
    ponegs[i] = S->NewVar();
    posels[i].resize(nInputs + nExtraInputs + nGates);
    for(int j = 0; j < nInputs + nExtraInputs + nGates; j++) {
      posels[i][j] = S->NewVar();
    }
    S->AMO(posels[i]);
  }
}

template <class T>
void SynthMan<T>::SortSels(vector<int> const &assignment) {
  for(int i = 0; i < nGates; i++) {
    if(!assignment[i+i] && !assignment[i+i+1]) {
      S->AddClause(-sels[i+i+1][nInputs + nExtraInputs - 1]);
      for(int j = 0; j < nInputs + nExtraInputs - 2; j++) {
        vector<int> vLits(j + 2);
        for(int k = 0; k <= j; k++) {
          vLits[k] = sels[i + i + 1][k];
        }
        vLits[j + 1] = -sels[i + i][j];
        S->AddClause(vLits);
      }
    }
  }
  for(int i = 0; i < nGates - 1; i++) {
    if(!assignment[i+i+2]) {
      for(int j = 0; j < nInputs + nExtraInputs - 2; j++) {
        vector<int> vLits(j + 2);
        for(int k = 0; k <= j; k++) {
          vLits[k] = sels[i + i][k];
        }
        vLits[j + 1] = -sels[i + i + 2][j];
        S->AddClause(vLits);
      }
    }
  }
  for(int i = 0; i < nGates - 1; i++) {
    if(!assignment[i+i+3] && assignment[i+i] == assignment[i+i+2]) {
      assert(!assignment[i+i+1]);
      if(!assignment[i+i]) {
        for(int j = 1; j < nInputs + nExtraInputs - 1; j++) {
          for(int k = 0; k < j; k++) {
            vector<int> vLits(k + 4);
            for(int l = 0; l <= k; l++) {
              vLits[l] = sels[i + i + 1][l];
            }
            vLits[k + 1] = -sels[i + i + 3][k];
            vLits[k + 2] = -sels[i + i][j];
            vLits[k + 3] = -sels[i + i + 2][j];
            S->AddClause(vLits);
          }
        }
      } else {
        for(int j = 0; j < nInputs + nExtraInputs - 1; j++) {
          vector<int> vLits(j + 2);
          for(int k = 0; k <= j; k++) {
            vLits[k] = sels[i + i + 1][k];
          }
          vLits[j + 1] = -sels[i + i + 3][j];
          S->AddClause(vLits);
        }
      }
    }
  }
  for(int i = 1; i < nGates; i++) {
    if(assignment[i+i] && !assignment[i+i+1]) {
      int j = assignment[i+i]-1;
      if(!assignment[j+j]) {
        for(int k = 0; k < nInputs + nExtraInputs - 1; k++) {
          vector<int> vLits(2);
          vLits[0] = -sels[j + j][k];
          vLits[1] = -sels[i + i + 1][k + 1];
          S->AddClause(vLits);
        }
      }
      if(!assignment[j+j+1]) {
        for(int k = 0; k < nInputs + nExtraInputs; k++) {
          vector<int> vLits(2);
          vLits[0] = -sels[j + j + 1][k];
          vLits[1] = -sels[i + i + 1][k];
          S->AddClause(vLits);
        }
      }
    }
  }
}

template <class T>
void SynthMan<T>::GenOne(vector<int> cands, vector<int> const &pos, vector<int> const&assignment) {
  cands.resize(nInputs + nExtraInputs + nGates);
  for(int i = 0; i < nGates; i++) {
    vector<int> fis(2);
    for(int k = 0; k <= 1; k++) {
      int r;
      if(!assignment[i+i+k]) {
        vector<int> tmps(nInputs + nExtraInputs + k - 1);
        for(int j = 0; j < nInputs + nExtraInputs + k - 1; j++) {
          tmps[j] = S->And2(cands[j + 1 - k], sels[i + i + k][j]);
        }
        r = S->OrN(tmps);
      } else {
        r = cands[nInputs + nExtraInputs + assignment[i+i+k] - 1];
      }
      fis[k] = S->Xor2(r, negs[i + i + k]);
    }
    cands[nInputs + nExtraInputs + i] = S->And2(fis[0], fis[1]);
  }
  for(int i = 0; i < nOutputs; i++) {
    vector<int> tmps(nInputs + nExtraInputs + nGates);
    for(int j = 0; j < nInputs + nExtraInputs + nGates; j++) {
      tmps[j] = S->And2(cands[j], posels[i][j]);
    }
    int r = S->OrN(tmps);
    S->Xor2(r, ponegs[i], pos[i]);
  }
}

template <class T>
aigman *SynthMan<T>::GetAig(vector<int> const &assignment) {
  aigman *aig = new aigman(nInputs + nExtraInputs, 0);
  vector<int> cands(nInputs + nExtraInputs + nGates);
  for(int i = 0; i < nInputs + nExtraInputs; i++) {
    cands[i] = i + 1;
  }
  for(int i = 0; i < nGates; i++) {
    vector<int> fis(2);
    for(int k = 0; k <= 1; k++) {
      int r;
      if(!assignment[i+i+k]) {
        int j = 0;
        for(; j < nInputs + nExtraInputs + k - 1; j++) {
          if(S->Value(sels[i + i + k][j])) {
            break;
          }
        }
        assert(j < nInputs + nExtraInputs + k - 1);
        r = cands[j + 1 - k];
      } else {
        r = cands[nInputs + nExtraInputs + assignment[i+i+k] - 1];
      }
      fis[k] = S->Value(negs[i + i + k])? (r << 1) ^ 1: r << 1;
    }
    cands[nInputs + nExtraInputs + i] = aig->newgate(fis[0], fis[1]);
  }
  aig->nPos = nOutputs;
  aig->vPos.resize(nOutputs);
  for(int i = 0; i < nOutputs; i++) {
    int j = 0;
    for(; j < nInputs + nExtraInputs + nGates; j++) {
      if(S->Value(posels[i][j])) {
        break;
      }
    }
    assert(j <= nInputs + nExtraInputs + nGates);
    int val = (j == nInputs + nExtraInputs + nGates)? 0: cands[j] << 1;
    aig->vPos[i] = val ^ (int)S->Value(ponegs[i]);
  }
  return aig;
}

template <class T>
aigman *SynthMan<T>::EnumSynth(int nGates_) {
  nGates = nGates_;
  vector<vector<int> > all;
  vector<int> tmp(nGates * 2 + 2);
  Enumerate(tmp, 1, all);
  for(auto const &v: all) {
    S = new T;
    GenSels(v);
    SortSels(v);
    for(int i = 0; i < (int)br.size(); i++) {
      bool fDc = true;
      for(int j = 0; j < (int)br[i].size(); j++) {
        if(!br[i][j]) {
          fDc = false;
          break;
        }
      }
      if(fDc) {
        continue;
      }
      vector<int> pis(nInputs);
      for(int k = 0; k < nInputs; k++) {
        pis[k] = (i >> k) & 1? S->one: S->zero;
      }
      vector<int> exins(nExtraInputs);
      for(int k = 0; k < nExtraInputs; k++) {
        exins[k] = (*sim)[i][k]? S->one: S->zero;
      }
      vector<int> pos(nOutputs);
      for(int k = 0; k < nOutputs; k++) {
        pos[k] = S->NewVar();
      }
      vector<int> tmps;
      for(int j = 0; j < (int)br[i].size(); j++) {
        if(br[i][j]) {
          vector<int> vLits(nOutputs);
          for(int k = 0; k < nOutputs; k++) {
            vLits[k] = (j >> k) & 1? pos[k]: -pos[k];
          }
          tmps.push_back(S->AndN(vLits));
        }
      }
      S->AddClause(tmps);
      pis.insert(pis.end(), exins.begin(), exins.end());
      GenOne(pis, pos, v);
    }
    if(S->Solve() == 1) {
      aigman *aig = GetAig(v);
      delete S;
      return aig;
    }
    delete S;
  }
  return NULL;
}

template <class T>
aigman *SynthMan<T>::ExEnumSynth(int nGates_) {
  assert(nGates_>= 0);
  aigman *aig = NULL;
  while(--nGates_ >= 0) {
    aigman *aig2 = EnumSynth(nGates_);
    if(!aig2) {
      break;
    }
    if(aig) {
      delete aig;
    }
    aig = aig2;
  }
  return aig;
}

template <class T>
aigman *SynthMan<T>::EnumSynth2(int nGates_) {
  nGates = nGates_;
  S = new T;
  GenSels();
  SortSels();
  for(int i = 0; i < (int)br.size(); i++) {
    bool fDc = true;
    for(int j = 0; j < (int)br[i].size(); j++) {
      if(!br[i][j]) {
        fDc = false;
        break;
      }
    }
    if(fDc) {
      continue;
    }
    vector<int> pis(nInputs);
    for(int k = 0; k < nInputs; k++) {
      pis[k] = (i >> k) & 1? S->one: S->zero;
    }
    vector<int> exins(nExtraInputs);
    for(int k = 0; k < nExtraInputs; k++) {
      exins[k] = (*sim)[i][k]? S->one: S->zero;
    }
    vector<int> pos(nOutputs);
    for(int k = 0; k < nOutputs; k++) {
      pos[k] = S->NewVar();
    }
    vector<int> tmps;
    for(int j = 0; j < (int)br[i].size(); j++) {
      if(br[i][j]) {
        vector<int> vLits(nOutputs);
        for(int k = 0; k < nOutputs; k++) {
          vLits[k] = (j >> k) & 1? pos[k]: -pos[k];
        }
        tmps.push_back(S->AndN(vLits));
      }
    }
    S->AddClause(tmps);
    pis.insert(pis.end(), exins.begin(), exins.end());
    GenOne(pis, pos);
  }
  vector<vector<int> > all;
  vector<int> tmp(nGates * 2 + 2);
  Enumerate(tmp, 1, all);
  vector<set<int> > cores;
  for(auto const &v: all) {
    vector<int> assumption;
    for(int i = 0; i < nGates; i++) {
      for(int k = 0; k < 2; k++) {
        if(v[i+i+k]) {
          int j2 = nInputs + nExtraInputs + v[i+i+k] + k - 2;
          assumption.push_back(sels[i + i + k][j2]);
          for(int j = 0; j < nInputs + nExtraInputs + i - 1; j++) {
            if(j != j2) {
              assumption.push_back(-sels[i + i + k][j]);
            }
          }
        } else {
          for(int j = nInputs + nExtraInputs + k - 1; j < nInputs + nExtraInputs + i - 1; j++) {
            assumption.push_back(-sels[i + i + k][j]);
          }
        }
      }
    }
    set<int> assumption_set(assumption.begin(), assumption.end());
    bool fIncluded = false;
    for(auto const &core: cores) {
      if(includes(assumption_set.begin(), assumption_set.end(), core.begin(), core.end())) {
        fIncluded = true;
      }
    }
    if(fIncluded) {
      continue;
    }
    set<int> core;
    int res = S->Solve(assumption, core);
    if(res == 1) {
      aigman *aig = GetAig();
      delete S;
      return aig;
    }
    if(res == -1) {
      cores.push_back(core);
    }
  }
  delete S;
  return NULL;
}

template <class T>
aigman *SynthMan<T>::ExEnumSynth2(int nGates_) {
  assert(nGates_>= 0);
  aigman *aig = NULL;
  while(--nGates_ >= 0) {
    aigman *aig2 = EnumSynth2(nGates_);
    if(!aig2) {
      break;
    }
    if(aig) {
      delete aig;
    }
    aig = aig2;
  }
  return aig;
}
