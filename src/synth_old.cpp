#include <cassert>
#include <algorithm>

#include "util.hpp"
#include "synth.hpp"

using namespace std;

template <class T>
void SynthMan<T>::GenSelsOld() {
  negs.clear();
  negs.resize(nGates * 2);
  sels.clear();
  sels.resize(nGates * 2);
  for(int i = 0; i < nGates * 2; i++) {
    negs[i] = S->NewVar();
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
void SynthMan<T>::SortSelsOld() {
  for(int i = 0; i < nGates; i++) {
    for(int j = 0; j < nInputs + nExtraInputs + i - 2; j++) {
      vector<int> vLits(j + 2);
      for(int k = 0; k <= j; k++) {
        vLits[k] = sels[i + i + 1][k];
      }
      vLits[j + 1] = -sels[i + i][j];
      S->AddClause(vLits);
    }
  }
  for(int i = 0; i < nGates - 1; i++) {
    for(int j = 0; j < nInputs + nExtraInputs + i - 2; j++) {
      vector<int> vLits(j + 2);
      for(int k = 0; k <= j; k++) {
        vLits[k] = sels[i + i][k];
      }
      vLits[j + 1] = -sels[i + i + 2][j];
      S->AddClause(vLits);
    }
  }
  for(int i = 0; i < nGates - 1; i++) {
    for(int j = 1; j < nInputs + nExtraInputs + i - 1; j++) {
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
  }
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
void SynthMan<T>::GenOneOld(vector<int> cands, vector<int> const &pos) {
  cands.resize(nInputs + nExtraInputs + nGates);
  for(int i = 0; i < nGates; i++) {
    vector<int> fis(2);
    for(int k = 0; k <= 1; k++) {
      vector<int> tmps(nInputs + nExtraInputs + i - 1);
      for(int j = 0; j < nInputs + nExtraInputs + i - 1; j++) {
        tmps[j] = S->And2(cands[j + 1 - k], sels[i + i + k][j]);
      }
      int r = S->OrN(tmps);
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
aigman *SynthMan<T>::GetAigOld() {
  aigman *aig = new aigman(nInputs + nExtraInputs, 0);
  vector<int> cands(nInputs + nExtraInputs + nGates);
  for(int i = 0; i < nInputs + nExtraInputs; i++) {
    cands[i] = i + 1;
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
      fis[k] = S->Value(negs[i + i + k])? (cands[j + 1 - k] << 1) ^ 1: cands[j + 1 - k] << 1;
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
