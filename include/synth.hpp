#pragma once

#include <aig.hpp>

#include "kissat_solver.hpp"
#include "cadical_solver.hpp"

template <class T>
class SynthMan {
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

  std::vector<std::vector<bool> > const *sim;
  int nExtraInputs;

  void GenSels();
  void SortSels();
  void GenOne(std::vector<int> cands, std::vector<int> const &pos);
  aigman *GetAig();

  void GenSelsOld();
  void SortSelsOld();
  void GenOneOld(std::vector<int> cands, std::vector<int> const &pos);
  aigman *GetAigOld();

  void GenSels(std::vector<int> const &assignment);
  void SortSels(std::vector<int> const &assignment);
  void GenOne(std::vector<int> cands, std::vector<int> const &pos, std::vector<int> const &assignment);
  aigman *GetAig(std::vector<int> const &assignment);

public:
  SynthMan(std::vector<std::vector<bool> > const &br, std::vector<std::vector<bool> > const *sim = NULL);

  aigman *Synth(int nGates_);
  aigman *ExSynth(int nGates_);

  aigman *EnumSynth(int nGates_);
  aigman *ExEnumSynth(int nGates_);

  aigman *EnumSynth2(int nGates_);
  aigman *ExEnumSynth2(int nGates_);
};

template class SynthMan<KissatSolver>;
template class SynthMan<CadicalSolver>;
