#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "cut.hpp"

template <typename T>
std::ostream &operator<<(std::ostream &os, std::vector<T> const &v) {
  std::string delim;
  for(T const &i: v) {
    os << delim << i;
    delim = ", ";
  }
  return os;
}

std::ostream &operator<<(std::ostream &os, Cut const &cut);

template <typename T>
void PrintVecWithIndex(std::vector<T> const &v, std::string prefix = "") {
  for(int i = 0; i < (int)v.size(); i++) {
    std::cout << prefix << i << " : " << v[i] << std::endl;
  }
}
