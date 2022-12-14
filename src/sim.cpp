#include <cassert>

#include "sim.hpp"

using namespace std;

static const unsigned long long basepats[] = {0xaaaaaaaaaaaaaaaaull,
                                              0xccccccccccccccccull,
                                              0xf0f0f0f0f0f0f0f0ull,
                                              0xff00ff00ff00ff00ull,
                                              0xffff0000ffff0000ull,
                                              0xffffffff00000000ull};

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
  vector<unsigned long long> inpats(basepats, basepats + 6);
  inpats.resize(aig.nPis);
  int ninpats = 1 << aig.nPis;
  int nsims = ninpats >> 6;
  if(!nsims) {
    nsims = 1;
  }
  for(int i = 0; i < nsims; i++) {
    for(int j = 0; j < aig.nPis - 6; j++) {
      inpats[j + 6] = (i >> j) & 1? 0xffffffffffffffffull: 0ull;
    }
    // run simulation
    aig.simulate(inpats);
    // get PO patterns
    vector<unsigned long long> outpats(aig.nPos);
    for(int j = 0; j < aig.nPos; j++) {
      outpats[j] = aig.getsim(aig.vPos[j]);
    }
    // get FI values
    vector<int> fivals(min(64, ninpats));
    for(int k = 0; k < (int)inputs.size(); k++) {
      for(int j = 0; j < 64 && j < ninpats; j++) {
        fivals[j] |= ((aig.vSims[inputs[k]] >> j) & 1) << k;
      }
    }
    // generate FO patterns
    for(int k = 0; k < nfopats; k++) {
      for(int j = 0; j < (int)outputs.size(); j++) {
        aig.vSims[outputs[j]] = (k >> j) & 1? 0xffffffffffffffffull: 0ull;
      }
      // run simulation from FOs
      aig.resimulate(focone);
      // exclude unnacceptable relation
      unsigned long long diff = 0ull;
      for(int j = 0; j < aig.nPos; j++) {
        diff |= outpats[j] ^ aig.getsim(aig.vPos[j]);
      }
      for(int j = 0; j < 64 && j < ninpats; j++) {
        br[fivals[j]][k] = (diff >> j) & 1? false: br[fivals[j]][k];
      }
    }
  }
}

void GetSim(aigman &aig, vector<int> const &inputs, vector<int> const &outputs, vector<vector<bool> > &sim) {
  assert(inputs.size() <= 16);
  aig.vSims.resize(aig.nObjs);
  // allocate
  int nfipats = 1 << inputs.size();
  sim.resize(nfipats);
  for(int i = 0; i < nfipats; i++) {
    sim[i].resize(outputs.size());
  }
  // generate input patterns
  for(int i = 0; i < (int)inputs.size() && i < 6; i++) {
    aig.vSims[inputs[i]] = basepats[i];
  }
  int nsims = nfipats >> 6;
  if(!nsims) {
    nsims = 1;
  }
  for(int i = 0; i < nsims; i++) {
    for(int j = 0; j < (int)inputs.size() - 6; j++) {
      aig.vSims[inputs[j + 6]] = (i >> j) & 1? 0xffffffffffffffffull: 0ull;
    }
    // run simulation
    aig.resimulate(outputs);
    // get output values
    for(int k = 0; k < (int)outputs.size(); k++) {
      for(int j = 0; j < 64 && j < nfipats; j++) {
        sim[j + i * 64][k] = (aig.vSims[outputs[k]] >> j) & 1;
      }
    }
  }
}
