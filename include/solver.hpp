#pragma once

#include <vector>

class Solver {
private:
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

  Solver(): fDirect(true) {}

  void Bimander(std::vector<int> const &vLits, int nbim);

  void PwNet(std::vector<int> vLits, std::vector<int> & res);
  void OddEvenSel4(std::vector<int> const &invars, std::vector<int> &outvars, int k);
  void PairwiseSel(std::vector<int> const &invars, std::vector<int> &outvars, int k);

public:
  virtual ~Solver() {}

  virtual int NewVar() = 0;
  virtual void AddClause(std::vector<int> const &vLits) = 0;
  virtual int Solve() = 0;
  virtual bool Value(int i) = 0;
  
  virtual void AMO(std::vector<int> const &vLits) = 0;
  virtual void Onehot(std::vector<int> const &vLits) = 0;
  virtual void AMK(std::vector<int> const &vLits, int k) = 0;

  void AddClause(int a);
  void AddClause(int a, int b);
  void AddClause(int a, int b, int c);

  void And2(int a, int b, int c);
  void Xor2(int a, int b, int c);
  void AndN(std::vector<int> vLits, int r);
  void OrN(std::vector<int> vLits, int r);
};
