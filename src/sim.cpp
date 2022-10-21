#include <cassert>

#include "sim.hpp"

using namespace std;

void GetBooleanRelation(aigman &aig, vector<int> const &inputs, vector<int> const &outputs, vector<vector<bool> > &br) {
  assert(aig.nPis <= 16);
  assert(inputs.size() <= 30);
  assert(outputs.size() <= 16);
  // set all relation acceptable
  int nfipats = 1 << inputs.size();
  int nfopats = 1 << outputs.size();
  br.resize(nfipats);
  for(int i = 0; i < nfipats; i++) {
    br[i].resize(nfopats, true);
  }
  // get focone
  vector<int> focone;
  aig.getfocone(outputs, focone);
  // generate PI patterns
  vector<unsigned long long> inpats;
  inpats.push_back(0xaaaaaaaaaaaaaaaaull);
  inpats.push_back(0xccccccccccccccccull);
  inpats.push_back(0xf0f0f0f0f0f0f0f0ull);
  inpats.push_back(0xff00ff00ff00ff00ull);
  inpats.push_back(0xffff0000ffff0000ull);
  inpats.push_back(0xffffffff00000000ull);
  inpats.resize(aig.nPis);
  int ninpats = aig.nPis <= 6? 1: (1 << (aig.nPis - 6));
  for(int i = 0; i < ninpats; i++) {
    for(int j = 0; j < aig.nPis - 6; j++) {
      if((i >> j) & 1) {
        inpats[j + 6] = 0xffffffffffffffffull;
      } else {
        inpats[j + 6] = 0ull;
      }
    }
    // run simulation
    aig.simulate(inpats);
    // get PO patterns
    vector<unsigned long long> outpats(aig.nPos);
    for(int j = 0; j < aig.nPos; j++) {
      outpats[j] = aig.getsim(aig.vPos[j]);
    }
    // get FI values
    vector<int> fivals(64);
    for(int k = 0; k < (int)inputs.size(); k++) {
      for(int j = 0; j < 64; j++) {
        fivals[j] |= ((aig.vSims[inputs[k]] >> j) & 1) << k;
      }
    }
    // generate FO patterns
    for(int k = 0; k < nfopats; k++) {
      for(int j = 0; j < (int)outputs.size(); j++) {
        if((k >> j) & 1) {
          aig.vSims[outputs[j]] = 0xffffffffffffffffull;
        } else {
          aig.vSims[outputs[j]] = 0ull;
        }
      }
      // run simulation from FOs
      aig.resimulate(focone);
      // exclude unnacceptable relation
      unsigned long long diff = 0ull;
      for(int j = 0; j < aig.nPos; j++) {
        diff |= outpats[j] ^ aig.getsim(aig.vPos[j]);
      }
      for(int j = 0; j < 64; j++) {
        if((diff >> j) & 1) {
          br[fivals[j]][k] = false;
        }
      }
    }
  }
}
