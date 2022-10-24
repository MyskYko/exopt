#pragma once

#include <aig.hpp>

void GetBooleanRelation(aigman &aig, std::vector<int> const &inputs, std::vector<int> const &outputs, std::vector<std::vector<bool> > &br);

void GetSim(aigman &aig, std::vector<int> const &inputs, std::vector<int> const &outputs, std::vector<std::vector<bool> > &sim);
