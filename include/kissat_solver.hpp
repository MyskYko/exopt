#pragma once

#include <iostream>
#include <cstdlib>

extern "C" {
  #include <kissat.h>
}

#include "solver.hpp"

class KissatSolver: public Solver {
private:
  kissat *S;
  int nClauses;

  void AddClause_(std::vector<int> const &vLits) {
    for(int i = 0; i < (int)vLits.size(); i++) {
      kissat_add(S, vLits[i]);
    }
    kissat_add(S, 0);
    nClauses++;
  }

  bool Value_(int i) {
    return kissat_value(S, i) > 0;
  }

  void AMO_(std::vector<int> const &vLits) {
    Bimander(vLits, 2);
  }

  void AMK_(std::vector<int> const &vLits, int k) {
    std::vector<int> res;
    OddEvenSel4(vLits, res, k + 1);
    AddClause(-res[k]);
  }

public:
  KissatSolver(): S(kissat_init()), nClauses(0) {}
  ~KissatSolver() {
    kissat_release(S);
  }

  int Solve() {
    int res = kissat_solve(S);
    return res == 10? 1: res == 20? -1: 0;
  }

  int Solve(std::vector<int> const &assumption, std::set<int> &core) {
    (void)assumption;
    (void)core;
    std::abort();
  }

  void PrintStat() {
    std::cout << "nVars: " << nVars << " nClauses: " << nClauses << std::endl;
  }
};
