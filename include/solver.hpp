#pragma once

#include <vector>

class Solver {
private:
  int nVars;

  void Comparator(int a, int b, int c, int d);
  void PwSplit(std::vector<int> const &a, std::vector<int> &b, std::vector<int> &c);
  void PwMerge(std::vector<int> const &a, std::vector<int> const &b, std::vector<int> &c);
  void PwSort(std::vector<int> const &a, std::vector<int> &d);

  static bool PreferDirectMerge(unsigned n, unsigned k);
  void DirectMerge(std::vector<int> const &in1, std::vector<int> const &in2, std::vector<int> &outvars, int k);
  void DirectCardClauses(std::vector<int> const &invars, int start, int pos, int j, std::vector<int> &args);
  void DirectNetwork(std::vector<int> const &invars, std::vector<int> &outvars, int k);
  void DirectCombine4(std::vector<int> const &x, std::vector<int> const &y, std::vector<int>& outvars, int k);
  void Comparator2(int x1, int x2, int y1, int y2);
  void OddEvenCombine(std::vector<int> const &in1, std::vector<int> const &in2, std::vector<int> &outvars, int k);
  void OddEvenMerge4(std::vector<int> const in[], std::vector<int> &outvars, int k);

  void PairwiseMerge(std::vector<int> const &x, std::vector<int> const &y, std::vector<int> &outvars, int k);
  void DirectPairwiseMerge(std::vector<int> const &in1, std::vector<int> const &in2, std::vector<int> &outvars, int k);

protected:
  bool fDirect;

  Solver(): nVars(0), fDirect(true), zero(0x7fffffff), one(-0x7fffffff) {}

  void Bimander(std::vector<int> const &vLits, int nbim);

  void PwNet(std::vector<int> vLits, std::vector<int> & res);
  void OddEvenSel4(std::vector<int> const &invars, std::vector<int> &outvars, int k);
  void PairwiseSel(std::vector<int> const &invars, std::vector<int> &outvars, int k);

  virtual void AddClause_(std::vector<int> const &vLits) = 0;
  virtual bool Value_(int i) = 0;
  virtual void AMO_(std::vector<int> const &vLits) = 0;
  virtual void AMK_(std::vector<int> const &vLits, int k) = 0;

public:
  const int zero;
  const int one;

  virtual ~Solver() {}

  virtual int Solve() = 0;

  inline int NewVar();

  inline void AddClause(std::vector<int> const &vLits);
  inline void AddClause(int a);
  inline void AddClause(int a, int b);
  inline void AddClause(int a, int b, int c);

  inline bool Value(int i);

  inline void AMO(std::vector<int> const &vLits);
  inline void Onehot(std::vector<int> const &vLits);
  inline void AMK(std::vector<int> const &vLits, int k);
  
  inline void And2(int a, int b, int c);
  inline void Xor2(int a, int b, int c);
  inline void AndN(std::vector<int> vLits, int r);
  inline void OrN(std::vector<int> vLits, int r);
};

int Solver::NewVar() {
  return ++nVars;
}

void Solver::AddClause(std::vector<int> const &vLits) {
  std::vector<int> vLits2(vLits.size());
  int j = 0;
  for(int i = 0; i < (int)vLits.size(); i++) {
    if(vLits[i] == one) {
      return;
    }
    if(vLits[i] == zero) {
      continue;
    }
    vLits2[j++] = vLits[i];
  }
  vLits2.resize(j);
  AddClause_(vLits2);
}
void Solver::AddClause(int a) {
  AddClause(std::vector<int>{a});
}
void Solver::AddClause(int a, int b) {
  AddClause(std::vector<int>{a, b});
}
void Solver::AddClause(int a, int b, int c) {
  AddClause(std::vector<int>{a, b, c});
}

bool Solver::Value(int i) {
  if(i == zero) {
    return false;
  }
  if(i == one) {
    return true;
  }
  return Value_(i);
}

void Solver::AMO(std::vector<int> const &vLits) {
  std::vector<int> vLits2(vLits.size());
  int j = 0;
  bool fOne = false;
  for(int i = 0; i < (int)vLits.size(); i++) {
    if(vLits[i] == one) {
      if(fOne) {
        AddClause_(std::vector<int>());
        return;
      }
      fOne = true;
      continue;
    }
    if(vLits[i] == zero) {
      continue;
    }
    vLits2[j++] = vLits[i];
  }
  if(fOne) {
    for(int i = 0; i < j; i++) {
      AddClause_(std::vector<int>{-vLits2[i]});
    }
    return;
  }
  vLits2.resize(j);
  AMO_(vLits2);
}
void Solver::Onehot(std::vector<int> const &vLits) {
  AMO(vLits);
  AddClause(vLits);
}
void Solver::AMK(std::vector<int> const &vLits, int k) {
  if(k < 0) {
    AddClause_(std::vector<int>());
    return;
  }
  std::vector<int> vLits2(vLits.size());
  int j = 0;
  for(int i = 0; i < (int)vLits.size(); i++) {
    if(vLits[i] == one) {
      if(!k) {
        AddClause_(std::vector<int>());
        return;
      }
      k--;
      continue;
    }
    if(vLits[i] == zero) {
      continue;
    }
    vLits2[j++] = vLits[i];
  }
  if(j <= k) {
    return;
  }
  if(!k) {
    for(int i = 0; i < j; i++) {
      AddClause_(std::vector<int>{-vLits2[i]});
    }
    return;
  }
  vLits2.resize(j);
  if(k == 1) {
    AMO_(vLits2);
  } else {
    AMK_(vLits2, k);
  }
}

void Solver::And2(int a, int b, int c) {
  AddClause(a, -c), AddClause(b, -c), AddClause(-a, -b, c);
}
void Solver::Xor2(int a, int b, int c) {
  AddClause(a, b, -c), AddClause(-a, -b, -c), AddClause(-a, b, c), AddClause(a, -b, c);
}
void Solver::AndN(std::vector<int> vLits, int r) {
  for(int i = 0; i < (int)vLits.size(); i++) {
    AddClause(vLits[i], -r);
  }
  for(int i = 0; i < (int)vLits.size(); i++) {
    vLits[i] = -vLits[i];
  }
  vLits.push_back(r);
  AddClause(vLits);
}
void Solver::OrN(std::vector<int> vLits, int r) {
  for(int i = 0; i < (int)vLits.size(); i++) {
    AddClause(-vLits[i], r);
  }
  vLits.push_back(-r);
  AddClause(vLits);
}
