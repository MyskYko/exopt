#pragma once

#include <aig.hpp>

struct Cut {
  std::vector<int> leaves;
  unsigned long long signature;

  Cut() {}
  Cut(int i): leaves(std::vector<int>{i}), signature(1ull << (i % 64)) {}
};

void CutEnumeration(aigman const &aig, std::vector<std::vector<Cut> > &cuts, unsigned cutsize = 6);
