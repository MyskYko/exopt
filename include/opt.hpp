#pragma once

#include <random>

#include "synth.hpp"

class OptMan {
private:
  aigman &aig;
  int cutsize;
  int windowsize;
  bool fAllDiv;
  bool fVerbose;
  std::mt19937 rg;

  int *nProblems;

  std::vector<std::pair<std::vector<int>, std::vector<int> > > vLarge;
  std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int> > > vWindows;

  template <typename T> void RemoveIncluded(T &s);
  template <typename T> void RemoveIncluded(T &s, bool fNoNewFo);
  bool Synthesize(SynthMan<KissatSolver> &synthman, int nGates, std::vector<int> const &inputs, std::vector<int> const &outputs, std::string prefix = "");

public:
  OptMan(aigman &aig, int cutsize, int windowsize, bool fAllDiv, int seed, bool fVerbose, int *nProblems = NULL);

  void Randomize();
  bool OptWindows();
  bool OptLarge();
};


