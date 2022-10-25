#include <cassert>

#include "util.hpp"
#include "synth.hpp"

using namespace std;

template <class T>
ExMan<T>::ExMan(vector<vector<bool> > const &br, vector<vector<bool> > const *sim): br(br), sim(sim) {
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
void ExMan<T>::GenSels() {
  negs.clear();
  negs.resize(nGates * 2);
  sels.clear();
  sels.resize(nGates * 2);
  for(int i = 0; i < nGates * 2; i++) {
    negs[i] = S->NewVar();
    for(int j = 0; j < nInputs + nExtraInputs + i/2 - 1; j++) {
      sels[i].push_back(S->NewVar());
    }
    S->Onehot(sels[i]);
  }
  ponegs.clear();
  ponegs.resize(nOutputs);
  posels.clear();
  posels.resize(nOutputs);
  for(int i = 0; i < nOutputs; i++) {
    ponegs[i] = S->NewVar();
    for(int j = 0; j < nInputs + nExtraInputs + nGates; j++) {
      posels[i].push_back(S->NewVar());
    }
    S->Onehot(posels[i]);
  }
}

template <class T>
void ExMan<T>::SortSels() {
  for(int i = 0; i < nGates; i++) {
    for(int j = 0; j < nInputs + nExtraInputs + i - 2; j++) {
      vector<int> vLits;
      for(int k = j; k >= 0; k--) {
        vLits.push_back(sels[i + i + 1][k]);
      }
      vLits.push_back(-sels[i + i][j]);
      S->AddClause(vLits);
    }
  }
  for(int i = 0; i < nGates - 1; i++) {
    for(int j = 0; j < nInputs + nExtraInputs + i - 2; j++) {
    vector<int> vLits;
      for(int k = j; k >= 0; k--) {
        vLits.push_back(sels[i + i][k]);
      }
      vLits.push_back(-sels[i + i + 2][j]);
      S->AddClause(vLits);
    }
  }
  for(int i = 0; i < nGates - 1; i++) {
    for(int j = 1; j < nInputs + nExtraInputs + i - 1; j++) {
      for(int k = j - 1; k >= 0; k--) {
        vector<int> vLits;
        vLits.push_back(-sels[i + i][j]);
        vLits.push_back(-sels[i + i + 2][j]);
        for(int l = k; l >= 0; l--) {
          vLits.push_back(sels[i + i + 1][l]);
        }
        vLits.push_back(-sels[i + i + 3][k]);
        S->AddClause(vLits);
      }
    }
  }
}

template <class T>
void ExMan<T>::GenOne(vector<int> cands, vector<int> const &pos) {
  for(int i = 0; i < nGates; i++) {
    vector<int> fis;
    for(int k = 0; k <= 1; k++) {
      vector<int> tmps;
      for(int j = 0; j < (int)cands.size() - 1; j++) {
        tmps.push_back(S->NewVar());
        S->And2(cands[j + 1 - k], sels[i + i + k][j], tmps.back());
      }
      int r = S->NewVar();
      S->OrN(tmps, r);
      fis.push_back(S->NewVar());
      S->Xor2(r, negs[i + i + k], fis.back());
    }
    cands.push_back(S->NewVar());
    S->And2(fis[0], fis[1], cands.back());
  }
  for(int i = 0; i < nOutputs; i++) {
    vector<int> tmps;
    for(int j = 0; j < (int)cands.size(); j++) {
      tmps.push_back(S->NewVar());
      S->And2(cands[j], posels[i][j], tmps.back());
    }
    int r = S->NewVar();
    S->OrN(tmps, r);
    S->Xor2(r, ponegs[i], pos[i]);
  }
}

template <class T>
aigman *ExMan<T>::GetAig() {
  aigman *aig = new aigman(nInputs + nExtraInputs, 0);
  vector<int> cands(nInputs + nExtraInputs);
  for(int i = 0; i < nInputs + nExtraInputs; i++) {
    cands[i] = i + 1;
  }
  for(int i = 0; i < nGates; i++) {
    vector<int> fis;
    for(int k = 0; k <= 1; k++) {
      int j = 0;
      for(; j < (int)cands.size() - 1; j++) {
        if(S->Value(sels[i + i + k][j])) {
          break;
        }
      }
      assert(j < (int)cands.size() - 1);
      fis.push_back(S->Value(negs[i + i + k])? (cands[j + 1 - k] << 1) ^ 1: cands[j + 1 - k] << 1);
    }
    cands.push_back(aig->newgate(fis[0], fis[1]));
  }
  for(int i = 0; i < nOutputs; i++) {
    int j = 0;
    for(; j < (int)cands.size(); j++) {
      if(S->Value(posels[i][j])) {
        break;
      }
    }
    assert(j < (int)cands.size());
    aig->vPos.push_back(S->Value(ponegs[i])? (cands[j] << 1) ^ 1: cands[j] << 1);
    aig->nPos++;
  }
  return aig;
}

template <class T>
aigman *ExMan<T>::Synth(int nGates_) {
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
        tmps.push_back(S->NewVar());
        S->AndN(vLits, tmps.back());
      }
    }
    S->AddClause(tmps);
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
aigman *ExMan<T>::ExSynth(int nGates_) {
  assert(nGates_);
  aigman *aig = NULL;
  while(--nGates_) {
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
