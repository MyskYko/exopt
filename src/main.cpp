#include <argparse/argparse.hpp>

#include <aig.hpp>

#include "opt.hpp"

using namespace std;

int main(int argc, char **argv) {
  argparse::ArgumentParser ap("exopt");
  ap.add_argument("input");
  ap.add_argument("output");
  ap.add_argument("-k", "--cutsize").default_value(8).scan<'i', int>();
  ap.add_argument("-n", "--windowsize").default_value(6).scan<'i', int>();
  ap.add_argument("-a", "--alldivisors").default_value(false).implicit_value(true);
  ap.add_argument("-r", "--numrounds").default_value(10).scan<'i', int>();
  ap.add_argument("-v", "--verbose").default_value(false).implicit_value(true);
  try {
    ap.parse_args(argc, argv);
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    cerr << ap;
    return 1;
  }
  string aigname = ap.get<string>("input");
  string outname = ap.get<string>("output");
  int cutsize = ap.get<int>("--cutsize");
  int windowsize = ap.get<int>("--windowsize");
  bool fAllDivisors = ap.get<bool>("--alldivisors");
  int numrounds = ap.get<int>("--numrounds");
  bool fVerbose = ap.get<bool>("--verbose");
  mt19937 rg;
  aigman aig_orig(aigname);
  aig_orig.supportfanouts();
  aigman aigout = aig_orig;
  for(int round = 0; round < numrounds; round++) {
    aigman aig = aig_orig;
    while(true) {
      bool fFirst;
      if(round == 0) {
        fFirst = true;
      } else if(round == 1) {
        fFirst = false;
      } else {
        fFirst = rg() % 2;
      }
      OptMan opt(aig, cutsize, windowsize, fAllDivisors, round, fVerbose);
      if(round > 1) {
        opt.Randomize();
      }
      if(fFirst && opt.OptWindows()) {
        continue;
      }
      if(opt.OptLarge()) {
        continue;
      }
      if(!fFirst && opt.OptWindows()) {
        continue;
      }
      break;
    }
    cout << aig.nGates << endl;
    if(aig.nGates < aigout.nGates) {
      aigout = aig;
    }
  }
  aigout.write(outname);
  return 0;
}
