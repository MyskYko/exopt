#pragma once

#include <aig.hpp>

#include "kissat_solver.hpp"

template <class T>
class ExMan {
private:
  T *S;
  std::vector<std::vector<bool> > const &br;
  int nInputs;
  int nOutputs;
  int nGates;
  std::vector<int> negs;
  std::vector<std::vector<int> > sels;
  std::vector<int> ponegs;
  std::vector<std::vector<int> > posels;

  void GenSels();
  void GenOne(std::vector<int> cands, std::vector<int> const &pos);
  aigman *GetAig();

public:
  ExMan(std::vector<std::vector<bool> > const &br);

  aigman *Solve(int nGates_);
};

template class ExMan<KissatSolver>;
