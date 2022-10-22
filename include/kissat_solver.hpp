#pragma once

extern "C" {
  #include <kissat.h>
}

#include "solver.hpp"

class KissatSolver: public Solver {
private:
  int nVars;
  kissat *S;

public:
  KissatSolver(): nVars(0) {
    S = kissat_init();
  }
  ~KissatSolver() {
    kissat_release(S);
  }

  int NewVar() {
    return ++nVars;
  }
  using Solver::AddClause;
  void AddClause(std::vector<int> const &vLits) {
    for(int i = 0; i < (int)vLits.size(); i++) {
      kissat_add(S, vLits[i]);
    }
    kissat_add(S, 0);
  }
  int Solve() {
    int res = kissat_solve(S);
    return res == 10? 1: res == 20? -1: 0;
  }
  bool Value(int i) {
    return kissat_value(S, i) > 0;
  }

  void AMO(std::vector<int> const &vLits) {
    Bimander(vLits, 2);
  }
  void Onehot(std::vector<int> const &vLits) {
    AMO(vLits);
    AddClause(vLits);
  }
  void AMK(std::vector<int> const &vLits, int k) {
    if((int)vLits.size() <= k) {
      return;
    }
    if(k == 1) {
      AMO(vLits);
      return;
    }
    std::vector<int> res;
    OddEvenSel4(vLits, res, k + 1);
    AddClause(-res[k]);
  }
};
