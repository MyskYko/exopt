#include <argparse/argparse.hpp>

#include "opt.hpp"
#include "rel.hpp"

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
  ap.add_argument("-g", "--numgates").scan<'i', int>();
  ap.add_argument("-d", "--dump").default_value(false).implicit_value(true);
  try {
    ap.parse_args(argc, argv);
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    cerr << ap;
    return 1;
  }
  string inname = ap.get<string>("input");
  string outname = ap.get<string>("output");
  int cutsize = ap.get<int>("--cutsize");
  int windowsize = ap.get<int>("--windowsize");
  bool fAllDivisors = ap.get<bool>("--alldivisors");
  int numrounds = ap.get<int>("--numrounds");
  bool fDump = ap.get<bool>("--dump");
  bool fVerbose = ap.get<bool>("--verbose");
  mt19937 rg;
  if(inname.substr(inname.find_last_of(".") + 1) == "rel") {
    int nGates = ap.get<int>("--numgates");
    vector<vector<bool> > br;
    vector<vector<bool> > *sim = NULL;
    ReadBooleanRelation(inname, br, sim, fVerbose);
    SynthMan<KissatSolver> synthman(br, sim);
    cout << "Synthesizing with at most " << nGates << " gates" << endl;
    aigman *aig = synthman.ExSynth(nGates + 1);
    if(aig) {
      aig->write(outname);
      cout << "Synthesized with " << aig->nGates << " gates" << endl;
      delete aig;
    } else {
      cout << "* Synthesis failed" << endl;
    }
    return 0;
  }
  aigman aig_orig(inname);
  aig_orig.supportfanouts();
  aigman aigout = aig_orig;
  int *nProblems = NULL;
  if(fDump) {
    nProblems = new int;
    *nProblems = 0;
  }
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
      OptMan opt(aig, cutsize, windowsize, fAllDivisors, round, fVerbose, nProblems);
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
    if(aig.nGates == aig_orig.nGates) {
      break;
    }
    if(aig.nGates < aigout.nGates) {
      aigout = aig;
    }
  }
  if(fDump) {
    delete nProblems;
  }
  cout << aigout.nGates << endl;
  aigout.write(outname);
  return 0;
}
