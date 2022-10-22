#include <iostream>
#include <cmath>
#include <cassert>

#include "synth.hpp"

using namespace std;

static inline int clog2(int n) {
  int t = 1;
  int count = 0;
  while(n > t) {
    t = t << 1;
    count++;
  }
  return count;
}

ExMan::ExMan(vector<vector<bool> > const &br): br(br) {
  nVars = 0;
  nInputs = clog2(br.size());
  nOutputs = clog2(br[0].size());
}
ExMan::~ExMan() {
}

void ExMan::AddClause(vector<int> const &vLits) {
  // string delim;
  // for(int i = 0; i < (int)vLits.size(); i++) {
  //   cout << delim << vLits[i];
  //   delim = " ";
  // }
  // cout << " 0" << endl;
  for(int i = 0; i < (int)vLits.size(); i++) {
    kissat_add(S, vLits[i]);
  }
  kissat_add(S, 0);
}
void ExMan::AddClause(int a) {
  // cout << a << " 0" << endl;
  kissat_add(S, a), kissat_add(S, 0);
}
void ExMan::AddClause(int a, int b) {
  // cout << a << " " << b << " 0" << endl;
  kissat_add(S, a), kissat_add(S, b), kissat_add(S, 0);
}
void ExMan::AddClause(int a, int b, int c) {
  // cout << a << " " << b << " " << c << " 0" << endl;
  kissat_add(S, a), kissat_add(S, b), kissat_add(S, c), kissat_add(S, 0);
}

void ExMan::Bimander(std::vector<int> const &vLits, int nbim) {
  std::vector<int> vLits2;
  int m = vLits.size() / nbim + vLits.size() % nbim;
  int nb = clog2(m);
  int c = nVars + 1;
  nVars += nb;
  for(int i = 0; i < m; i++) {
    vLits2.clear();
    for(int j = 0; j < nbim && i*nbim + j < (int)vLits.size(); j++) {
      vLits2.push_back(vLits[i*nbim + j]);
    }
    if(vLits2.size() > 1) {
      for(int p = 0; p < (int)vLits2.size(); p++) {
        for(int q = p+1; q < (int)vLits2.size(); q++) {
          AddClause(-vLits2[p], -vLits2[q]);
        }
      }
    }
    for(int k = 0; k < nb; k++) {
      int b = 1 << k;
      if(i & b) {
        for(int j = 0; j < (int)vLits2.size(); j++) {
          AddClause(-vLits2[j], c + k);
        }
      } else {
        for(int j = 0; j < (int)vLits2.size(); j++) {
          AddClause(-vLits2[j], -(c + k));
        }
      }
    }
  }
}
void ExMan::Onehot(vector<int> const &vLits) {
  Bimander(vLits);
  AddClause(vLits);
}

void ExMan::And2(int a, int b, int c) {
  AddClause(a, -c), AddClause(b, -c), AddClause(-a, -b, c);
}
void ExMan::Xor2(int a, int b, int c) {
  AddClause(a, b, -c), AddClause(-a, -b, -c), AddClause(-a, b, c), AddClause(a, -b, c);
}
void ExMan::AndN(std::vector<int> vLits, int r) {
  for(int i = 0; i < (int)vLits.size(); i++) {
    AddClause(vLits[i], -r);
  }
  for(int i = 0; i < (int)vLits.size(); i++) {
    vLits[i] = -vLits[i];
  }
  vLits.push_back(r);
  AddClause(vLits);
}
void ExMan::OrN(std::vector<int> vLits, int r) {
  for(int i = 0; i < (int)vLits.size(); i++) {
    AddClause(-vLits[i], r);
  }
  vLits.push_back(-r);
  AddClause(vLits);
}

void ExMan::GenSels() {
  negs.clear();
  negs.resize(nGates * 2);
  sels.clear();
  sels.resize(nGates * 2);
  for(int i = 0; i < nGates * 2; i++) {
    negs[i] = ++nVars;
    for(int j = 0; j < nInputs + i/2; j++) {
      sels[i].push_back(++nVars);
    }
    Onehot(sels[i]);
  }
  ponegs.clear();
  ponegs.resize(nOutputs);
  posels.clear();
  posels.resize(nOutputs);
  for(int i = 0; i < nOutputs; i++) {
    ponegs[i] = ++nVars;
    for(int j = 0; j < nInputs + nGates; j++) {
      posels[i].push_back(++nVars);
    }
    Onehot(posels[i]);
  }
}

void ExMan::GenOne(vector<int> cands, vector<int> const &pos) {
  for(int i = 0; i < nGates; i++) {
    vector<int> fis;
    for(int k = i + i; k <= i + i + 1; k++) {
      vector<int> tmps;
      for(int j = 0; j < (int)cands.size(); j++) {
        tmps.push_back(++nVars);
        And2(cands[j], sels[k][j], tmps.back());
      }
      int r = ++nVars;
      OrN(tmps, r);
      fis.push_back(++nVars);
      Xor2(r, negs[k], fis.back());
    }
    cands.push_back(++nVars);
    And2(fis[0], fis[1], cands.back());
  }
  for(int i = 0; i < nOutputs; i++) {
    vector<int> tmps;
    for(int j = 0; j < (int)cands.size(); j++) {
      tmps.push_back(++nVars);
      And2(cands[j], posels[i][j], tmps.back());
    }
    int r = ++nVars;
    OrN(tmps, r);
    Xor2(r, ponegs[i], pos[i]);
  }
}

aigman *ExMan::GetAig() {
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
        if(kissat_value(S, sels[k][j]) > 0) {
          break;
        }
      }
      assert(j < (int)cands.size());
      if(kissat_value(S, negs[k]) > 0) {
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
      if(kissat_value(S, posels[i][j]) > 0) {
        break;
      }
    }
    assert(j < (int)cands.size());
    if(kissat_value(S, ponegs[i]) > 0) {
      aig->vPos.push_back((cands[j] << 1) ^ 1);
    } else {
      aig->vPos.push_back(cands[j] << 1);
    }
    aig->nPos++;
  }
  
  return aig;
}

aigman *ExMan::Solve(int nGates_) {
  nGates = nGates_;
  S = kissat_init();

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
      pis[k] = ++nVars;
    }
    vector<int> pos(nOutputs);
    for(int k = 0; k < nOutputs; k++) {
      pos[k] = ++nVars;
    }
    for(int k = 0; k < nInputs; k++) {
      if((i >> k) & 1) {
        AddClause(pis[k]);
      } else {
        AddClause(-pis[k]);
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
        tmps.push_back(++nVars);
        AndN(vLits, tmps.back());
      }
    }
    if(tmps.size() > 1) {
      int r = ++nVars;
      OrN(tmps, r);
      AddClause(r);
    } else {
      AddClause(tmps[0]);
    }

    GenOne(pis, pos);
  }

  int res = kissat_solve(S);
  aigman *aig = NULL;
  if(res == 10) {
    aig = GetAig();
  }

  kissat_release(S);

  return aig;
}
