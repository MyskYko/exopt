#pragma once

#include <string>
#include <vector>

void ReadBooleanRelation(std::string fname, std::vector<std::vector<bool> > &br, std::vector<std::vector<bool> > *&sim, bool fVerbose);

void WriteBooleanRelation(std::string fname, std::vector<std::vector<bool> > const &br, std::vector<std::vector<bool> > const *sim);
