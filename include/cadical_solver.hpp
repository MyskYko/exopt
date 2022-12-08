#pragma once

#include <cadical.hpp>

#include "solver.hpp"

class CadicalSolver: public Solver {
private:
  CaDiCaL::Solver *S;

  void AddClause_(std::vector<int> const &vLits) {
    for(int i = 0; i < (int)vLits.size(); i++) {
      S->add(vLits[i]);
    }
    S->add(0);
  }

  bool Value_(int i) {
    return S->val(i) > 0;
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
  CadicalSolver(): S(new CaDiCaL::Solver) {}
  ~CadicalSolver() {
    delete S;
  }

  int Solve() {
    int res = S->solve();
    return res == 10? 1: res == 20? -1: 0;
  }

  int Solve(std::vector<int> const &assumption, std::set<int> &core) {
    for(int i: assumption) {
      S->assume(i);
    }
    int res = S->solve();
    for(int i: assumption) {
      if(S->failed(i)) {
        core.insert(i);
      }
    }
    return res == 10? 1: res == 20? -1: 0;
  }
};
