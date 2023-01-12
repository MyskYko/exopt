#include "ioutil.hpp"

std::ostream &operator<<(std::ostream &os, Cut const &cut) {
  os << "{" << cut.leaves << "}";
  return os;
}
