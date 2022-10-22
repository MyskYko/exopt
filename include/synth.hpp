#pragma once

extern "C" {
  #include <kissat.h>
}

#include <aig.hpp>

class ExMan {
private:
  kissat *S;
  int nVars;

  void AddClause(std::vector<int> const &vLits);
  void AddClause(int a);
  void AddClause(int a, int b);
  void AddClause(int a, int b, int c);

  void Bimander(std::vector<int> const &vLits, int nbim = 2);
  void Onehot(std::vector<int> const &vLits);

  void And2(int a, int b, int c);
  void Xor2(int a, int b, int c);
  void AndN(std::vector<int> vLits, int r);
  void OrN(std::vector<int> vLits, int r);

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
  ~ExMan();
  
  aigman *Solve(int nGates_);
};
