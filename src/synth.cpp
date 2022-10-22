#include <iostream>
#include <cmath>
#include <cassert>

#include "util.hpp"
#include "synth.hpp"

using namespace std;

template <class T>
ExMan<T>::ExMan(vector<vector<bool> > const &br): br(br) {
  nInputs = clog2(br.size());
  nOutputs = clog2(br[0].size());
}

template <class T>
void ExMan<T>::GenSels() {
  negs.clear();
  negs.resize(nGates * 2);
  sels.clear();
  sels.resize(nGates * 2);
  for(int i = 0; i < nGates * 2; i++) {
    negs[i] = S->NewVar();
    for(int j = 0; j < nInputs + i/2; j++) {
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
    for(int j = 0; j < nInputs + nGates; j++) {
      posels[i].push_back(S->NewVar());
    }
    S->Onehot(posels[i]);
  }
}

template <class T>
void ExMan<T>::GenOne(vector<int> cands, vector<int> const &pos) {
  for(int i = 0; i < nGates; i++) {
    vector<int> fis;
    for(int k = i + i; k <= i + i + 1; k++) {
      vector<int> tmps;
      for(int j = 0; j < (int)cands.size(); j++) {
        tmps.push_back(S->NewVar());
        S->And2(cands[j], sels[k][j], tmps.back());
      }
      int r = S->NewVar();
      S->OrN(tmps, r);
      fis.push_back(S->NewVar());
      S->Xor2(r, negs[k], fis.back());
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
  aigman *aig = new aigman(nInputs, 0);
  vector<int> cands(nInputs);
  for(int i = 0; i < nInputs; i++) {
    cands[i] = i + 1;
  }
  for(int i = 0; i < nGates; i++) {
    vector<int> fis;
    for(int k = i + i; k <= i + i + 1; k++) {
      int j = 0;
      for(; j < (int)cands.size(); j++) {
        if(S->Value(sels[k][j])) {
          break;
        }
      }
      assert(j < (int)cands.size());
      if(S->Value(negs[k])) {
        aig->vObjs.push_back((cands[j] << 1) ^ 1);
      } else {
        aig->vObjs.push_back(cands[j] << 1);
      }
    }
    cands.push_back(aig->nObjs++);
    aig->nGates++;
  }
  for(int i = 0; i < nOutputs; i++) {
    int j = 0;
    for(; j < (int)cands.size(); j++) {
      if(S->Value(posels[i][j])) {
        break;
      }
    }
    assert(j < (int)cands.size());
    if(S->Value(ponegs[i])) {
      aig->vPos.push_back((cands[j] << 1) ^ 1);
    } else {
      aig->vPos.push_back(cands[j] << 1);
    }
    aig->nPos++;
  }
  
  return aig;
}

template <class T>
aigman *ExMan<T>::Solve(int nGates_) {
  nGates = nGates_;
  S = new T;
  GenSels();
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
      pis[k] = S->NewVar();
    }
    vector<int> pos(nOutputs);
    for(int k = 0; k < nOutputs; k++) {
      pos[k] = S->NewVar();
    }
    for(int k = 0; k < nInputs; k++) {
      if((i >> k) & 1) {
        S->AddClause(pis[k]);
      } else {
        S->AddClause(-pis[k]);
      }
    }
    vector<int> tmps;
    for(int j = 0; j < (int)br[i].size(); j++) {
      if(br[i][j]) {
        vector<int> vLits(nOutputs);
        for(int k = 0; k < nOutputs; k++) {
          if((j >> k) & 1) {
            vLits[k] = pos[k];
          } else {
            vLits[k] = -pos[k];
          }
        }
        tmps.push_back(S->NewVar());
        S->AndN(vLits, tmps.back());
      }
    }
    if(tmps.size() > 1) {
      int r = S->NewVar();
      S->OrN(tmps, r);
      S->AddClause(r);
    } else {
      S->AddClause(tmps[0]);
    }
    GenOne(pis, pos);
  }
  aigman *aig = NULL;
  if(S->Solve() == 1) {
    aig = GetAig();
  }
  delete S;
  return aig;
}
