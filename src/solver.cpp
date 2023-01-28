#include <cassert>

#include "util.hpp"
#include "solver.hpp"

using namespace std;

void Solver::Pairwise(vector<int> const &vLits) {
  for(int i = 1; i < (int)vLits.size(); i++) {
    for(int j = 0; j < i; j++) {
      AddClause_(vector<int>{-vLits[i], -vLits[j]});
    }
  }
}

void Solver::Bimander(vector<int> const &vLits, int nbim) {
  vector<int> vLits2;
  int m = vLits.size() / nbim + vLits.size() % nbim;
  int nb = clog2(m);
  vector<int> cv(nb);
  for(int i = 0; i < nb; i++) {
    cv[i] = NewVar();
  }
  for(int i = 0; i < m; i++) {
    vLits2.clear();
    for(int j = 0; j < nbim && i*nbim + j < (int)vLits.size(); j++) {
      vLits2.push_back(vLits[i*nbim + j]);
    }
    if(vLits2.size() > 1) {
      for(int p = 0; p < (int)vLits2.size(); p++) {
        for(int q = p+1; q < (int)vLits2.size(); q++) {
          AddClause_(vector<int>{-vLits2[p], -vLits2[q]});
        }
      }
    }
    for(int k = 0; k < nb; k++) {
      int b = 1 << k;
      if(i & b) {
        for(int j = 0; j < (int)vLits2.size(); j++) {
          AddClause_(vector<int>{-vLits2[j], cv[k]});
        }
      } else {
        for(int j = 0; j < (int)vLits2.size(); j++) {
          AddClause_(vector<int>{-vLits2[j], -cv[k]});
        }
      }
    }
  }
}

void Solver::Comparator(int a, int b, int c, int d) {
  AddClause(-c, a, b), AddClause(c, -a), AddClause(c, -b);
  AddClause(d, -a, -b), AddClause(-d, a), AddClause(-d, b);
}
void Solver::PwSplit(vector<int> const &a, vector<int> &b, vector<int> &c) {
  int n = a.size() / 2;
  b.resize(n), c.resize(n);
  for(int i = 0; i < n; i++) {
    b[i] = NewVar(), c[i] = NewVar();
    Comparator(a[i + i], a[i + i + 1], b[i], c[i]);
  }
}
void Solver::PwMerge(std::vector<int> const &a, std::vector<int> const &b, std::vector<int> &c) {
  std::vector<int> a_next, b_next, d, e;
  int n = a.size();
  if(n == 1) {
    c.push_back(a[0]), c.push_back(b[0]);
    return;
  }
  a_next.resize(n / 2), b_next.resize(n / 2);
  for(int i = 0; i < n / 2; i++) {
    a_next[i] = a[i + i], b_next[i] = b[i + i];
  }
  PwMerge(a_next, b_next, d);
  for(int i = 0; i < n / 2; i++) {
    a_next[i] = a[i + i + 1], b_next[i] = b[i + i + 1];
  }
  PwMerge(a_next, b_next, e);
  c.resize(n + n);
  c[0] = d[0];
  for(int i = 0; i < n - 1; i++) {
    c[i + i + 1] = NewVar(), c[i + i + 2] = NewVar();
    Comparator(e[i], d[i + 1], c[i + i + 1], c[i + i + 2]);
  }
  c[n + n - 1] = e[n - 1];
}
void Solver::PwSort(std::vector<int> const &a, std::vector<int> &d) {
  if(a.size() == 1) {
    d.push_back(a[0]);
    return;
  }
  std::vector<int> b, c, b_next, c_next;
  PwSplit(a, b, c);
  PwSort(b, b_next);
  b.clear();
  PwSort(c, c_next);
  c.clear();
  PwMerge(b_next, c_next, d);
}
void Solver::PwNet(std::vector<int> vLits, std::vector<int> &res) {
  int n = pow2roundup(vLits.size());
  if((int)vLits.size() != n) {
    vLits.resize(n, zero);
  }
  PwSort(vLits, res);
}

bool Solver::PreferDirectMerge(unsigned n, unsigned k) {
  static const unsigned minTest = 94, maxTest = 183;
  static const unsigned short nBound[] = {94+171, 95+150, 96+177, 97+156, 98+135, 99+126,
                                          100+141,101+128,102+119,103+110,104+121,105+112,106+103,107+98, 108+109,109+100,
                                          110+95, 111+90, 112+97, 113+92, 114+84, 115+82, 116+89, 117+84, 118+76, 119+74,
                                          120+81, 121+76, 122+71, 123+70, 124+73, 125+69, 126+67, 127+63, 128+69, 129+64,
                                          130+63, 131+58, 132+62, 133+60, 134+56, 135+54, 136+57, 137+56, 138+52, 139+50,
                                          140+53, 141+52, 142+48, 143+46, 144+49, 145+48, 146+44, 147+43, 148+45, 149+44,
                                          150+41, 151+39, 152+42, 153+40, 154+39, 155+38, 156+38, 157+37, 158+35, 159+34,
                                          160+37, 161+33, 162+32, 163+31, 164+33, 165+32, 166+28, 167+27, 168+29, 169+28,
                                          170+27, 171+26, 172+26, 173+24, 174+23, 175+22, 176+22, 177+21, 178+20, 179+19,
                                          180+21, 181+20, 182+16, 183+15, 184+17, 185+16, 186+15, 187+14, 188+13, 189+12,
                                          190+11, 191+10, 192+10, 193+9,  194+8,  195+7,  196+6,  197+5,  198+4,  199+3,
                                          200+2,  201+1};
  if(k > n) {
    k = n;
  }
  return k == 1 || (k >= 4 && k < minTest && n >= 10) || (k >= minTest && k <= maxTest && n < nBound[k-minTest]);
}
void Solver::DirectMerge(vector<int> const &in1, vector<int> const &in2, std::vector<int> &outvars, int k) {
  int a = min(k, (int)in1.size());
  int b = min(k, (int)in2.size());
  int c = min(k, a + b);
  if(b == 0) {
    for(int i = 0; i < c; i++) {
      outvars.push_back(in1[i]);
    }
  } else if(a == 0) {
    for(int i = 0; i < c; i++) {
      outvars.push_back(in2[i]);
    }
  } else {
    for(int i = 0; i < c; i++) {
      outvars.push_back(NewVar());
    }
    for(int i = 0; i < a; i++) {
      AddClause(-in1[i], outvars[i]);
    }
    for(int i = 0; i < b; i++) {
      AddClause(-in2[i], outvars[i]);
    }
    for(int j = 0; j < b; j++) {
      for(int i = 0; i < min(a, c-j-1); i++) {
        AddClause(-in1[i], -in2[j], outvars[i+j+1]);
      }
    }
  }
}
void Solver::DirectCardClauses(vector<int> const &invars, int start, int pos, int j, vector<int> &args) {
  int n = invars.size();
  if(pos == j) {
    AddClause(args);
  } else {
    for(int i = start; i <= n - (j-pos); i++) {
      args[pos] = -invars[i];
      DirectCardClauses(invars, i+1, pos+1, j, args);
    }
  }
}
void Solver::DirectNetwork(vector<int> const &invars, vector<int> &outvars, int k) {
  assert(outvars.empty());
  int n = invars.size();
  if(k == 0 || k > n) {
    k = n;
  }
  for(int i = 0; i < k; i++) {
    outvars.push_back(NewVar());
  }
  for(int j = 1; j <= k; j++) {
    vector<int> args(j);
    args.push_back(outvars[j-1]);
    DirectCardClauses(invars, 0, 0, j, args);
  }
}
void Solver::DirectCombine4(std::vector<int> const &x, std::vector<int> const &y, std::vector<int>& outvars, int k) {
  int a = x.size(), b = y.size();
  assert(a >= b), assert(a <= b+4), assert(a >= 2), assert(b >= 1);
  if(k > a + b) {
    k = a + b;
  }
  outvars.push_back(x[0]);
  int last = (k < a + b || k % 2 == 1 || a == b + 2)? k: k-1;
  for(int i = 0, j = 1; j < last; j++, i = j/2) {
    int ret = NewVar();
    outvars.push_back(ret);
    if(j % 2 == 0) {
      if(i + 1 < a && i < b + 2) {
        if(i >= 2) {
          AddClause(-x[i+1], -y[i-2], ret);
        } else {
          AddClause(-x[i+1], ret);
        }
      }
      if(i < a && i < b + 1) {
        AddClause(-x[i], -y[i-1], ret);
      }
    } else {
      if(i > 0 && i + 2 < a) {
        AddClause(-x[i+2], ret);
      }
      if(i < b) {
        AddClause(-y[i], ret);
      }
      if(i + 1 < a && i < b + 1) {
        if(i > 0) {
          AddClause(-x[i+1], -y[i-1], ret);
        } else {
          AddClause(-x[i+1], ret);
        }
      }
    }
  }
  if(k == a + b && k % 2 == 0 && a != b + 2) {
    outvars.push_back(a == b ? y[b-1] : x[a-1]);
  }
  if(k < a + b) {
    AddClause(-x[a-1], -y[b-1]);
    if(k + 1 < a + b) {
      AddClause(-x[a-2], -y[b-1]);
      if(b >= 2) {
        AddClause(-x[a-1], -y[b-2]);
      }
    }
  }
}
void Solver::Comparator2(int x1, int x2, int y1, int y2) {
  AddClause(-x1, y1), AddClause(-x2, y1), AddClause(-x1, -x2, y2);
}
void Solver::OddEvenCombine(vector<int> const &in1, vector<int> const &in2, std::vector<int> &outvars, int k) {
  int a = in1.size(), b = in2.size();
  if(k > a + b) {
    k = a + b;
  }
  outvars.push_back(in1[0]);
  for(int i = 0; i < (k-1)/2; i++) {
    outvars.push_back(NewVar());
    outvars.push_back(NewVar());
    Comparator2(in2[i], in1[i+1], outvars[i*2+1], outvars[i*2+2]);
  }
  if(k % 2 == 0) {
    if(k < a + b) {
      int ret = NewVar();
      outvars.push_back(ret);
      AddClause(-in2[k/2-1], ret);
      AddClause(-in1[k/2], ret);
    } else if(a == b) {
      outvars.push_back(in2[k/2-1]);
    } else {
      outvars.push_back(in1[k/2]);
    }
  }
  if(k < a + b) {
    AddClause(-in1[a-1], -in2[b-1]);
  }
}
void Solver::OddEvenMerge4(vector<int> const in[], vector<int> &outvars, int k) {
  int nn[4];
  nn[0] = in[0].size(), nn[1] = in[1].size(), nn[2] = in[2].size(), nn[3] = in[3].size();
  assert(nn[0] > 0); assert(nn[0] >= nn[1]); assert(nn[1] >= nn[2]); assert(nn[2] >= nn[3]);
  k = min(k, nn[0] + nn[1] + nn[2] + nn[3]);
  for(int j = 0; j < 4; j++) {
    if(nn[j] > k){
      nn[j] = k;
    }
  }
  if(nn[1] == 0) {
    for(int i = 0; i < nn[0]; i++) {
      outvars.push_back(in[0][i]);
    }
  } else if(nn[0] == 1) {
    vector<int> invars;
    for(int j = 0; j < 4; j++) {
      if(nn[j] > 0) {
        invars.push_back(in[j][0]);
      }
    }
    DirectNetwork(invars, outvars, k);
  } else {
    vector<int> even_odd[2][4], x, y;
    for(int j = 0; j < 4; j++) {
      for(int i = 0; i < nn[j]; i++) {
        even_odd[i%2][j].push_back(in[j][i]);
      }
    }
    OddEvenMerge4(even_odd[0], x, k/2+2);
    OddEvenMerge4(even_odd[1], y, k/2);
    if(nn[2] > 0) {
      DirectCombine4(x, y, outvars, k);
    } else {
      OddEvenCombine(x, y, outvars, k);
    }
  }
}
void Solver::OddEvenSel4(std::vector<int> const &invars, std::vector<int> &outvars, int k) {
  int n = invars.size();
  assert(k <= n);
  if(k == 0) {
    for(int i = 0; i < n; i++) {
      AddClause(-invars[i]);
    }
  } else if(n == 1) {
    outvars.push_back(invars[0]);
  } else if(n > 1) {
    if(n <= 4 || (fDirect && (k <= 1 || (k == 2 && n <= 9) || n <= 6))) {
      DirectNetwork(invars, outvars,k);
    } else {
      int nn[4], kk[4];
      int p2 = pow2roundup((k+5)/6);
      if(n >= 8 && 4 * p2 <= n)  {
        nn[1] = nn[2] = nn[3] = p2;
      } else if(n < 8 || k == n) {
        nn[1] = (n+2)/4, nn[2] = (n+1)/4, nn[3] = n/4;
      } else {
        nn[1] = nn[2] = nn[3] = k/4;
      }
      nn[0] = n - nn[1] - nn[2] - nn[3];
      vector<int> in[4], out[4];
      for(int base = 0, j = 0; j < 4; base += nn[j], j++) {
        for(int i = 0; i < nn[j]; i++) {
          in[j].push_back(invars[base+i]);
        }
      }
      for(int j = 0; j < 4; j++) {
        kk[j] = min(k, nn[j]);
        OddEvenSel4(in[j], out[j], kk[j]);
      }
      if(fDirect && PreferDirectMerge(kk[0]+kk[1]+kk[2]+kk[3], k)) {
        vector<int> out1,out2;
        DirectMerge(out[0], out[1], out1, min(kk[0]+kk[1], k));
        DirectMerge(out[2], out[3], out2, min(kk[2]+kk[3], k));
        DirectMerge(out1, out2, outvars, k);
      } else {
        OddEvenMerge4(out, outvars, k);
      }
    }
  }
}

void Solver::DirectPairwiseMerge(vector<int> const &in1, vector<int> const &in2, vector<int> &outvars, int k) {
  int a = min(k, (int)in1.size()), b = min(k, (int)in2.size()), c = min(k, a + b);
  if(b == 0) {
    for(int i = 0; i < c; i++) {
      outvars.push_back(in1[i]);
    }
  } else {
    for(int i = 0; i < c; i++) {
      outvars.push_back(NewVar());
    }
    for(int i = 0; i < a; i++) {
      AddClause(-in1[i], outvars[i]);
    }
    for(int i = 0; i < min(b, c/2); i++) {
      AddClause(-in2[i], outvars[2*i+1]);
    }
    for(int j = 0; j < b; j++) {
      for(int i = j+1; i < min(a, c-j-1); i++) {
        AddClause(-in1[i], -in2[j], outvars[i+j+1]);
      }
    }
  }
}
void Solver::PairwiseMerge(vector<int> const &x, vector<int> const &y, vector<int> &outvars, int k) {
  int n1 = x.size(), n2 = y.size();
  vector<int> xi = x, yi = y;
  int h = pow2roundup(n1);
  for(; n2 < k/2; n2++) {
    yi.push_back(0);
  }
  while(h > 1) {
    h = h/2;
    for(int j = 0; j < n2; j++) {
      if(j + h < n1) {
        int xout, yout;
        if(xi[j+h] == 0) {
          yout = yi[j], xout = xi[j+h];
        } else if(yi[j] == 0) {
          yout = xi[j+h], xout = yi[j];
        } else {
          xout = NewVar(), yout = NewVar();
          Comparator2(xi[j+h], yi[j], yout, xout);
        }
        xi[j+h] = xout, yi[j] = yout;
      }
    }
  }
  for(int j = 0; j < k; j++) {
    outvars.push_back(j % 2 == 0? xi[j/2]: yi[j/2]);
    assert(outvars[j]);
  }
  for(int j = (k+1)/2; j < n1; j++) {
    if(xi[j] != 0) {
      AddClause(-xi[j]);
    }
  }
}
void Solver::PairwiseSel(vector<int> const &invars, vector<int> &outvars, int k) {
  assert(outvars.empty());
  int n = invars.size();
  k = min(k, n);
  if(k == 0) {
    for(int i = 0; i < n; i++) {
      AddClause(-invars[i]);
    }
  } else if(n == 1) {
    outvars.push_back(invars[0]);
  } else if (n <= 2 || (fDirect && (k <= 1 || (k == 2 && n <= 9) || n <= 6))) {
    DirectNetwork(invars, outvars, k);
  } else {
    int n1, n2;
    int p2 = pow2roundup((k+2)/3);
    n2 = (n <= 7 ? n/2 : p2 <= k/2)? p2: k - p2;
    n1 = n - n2;
    vector<int> x, y;
    for(int i = 0; i < n2; i++) {
      x.push_back(NewVar());
      y.push_back(NewVar());
      Comparator2(invars[2*i], invars[2*i+1], x[i], y[i]);
    }
    for(int i = n2; i < n1; i++) {
      x.push_back(invars[n2+i]);
    }
    vector<int> x_, y_;
    PairwiseSel(x, x_, min(k, n1));
    PairwiseSel(y, y_, min(k/2, n2));
    if(fDirect && PreferDirectMerge(x_.size() + y_.size(), k)) {
      DirectPairwiseMerge(x_, y_, outvars, k);
    } else {
      PairwiseMerge(x_, y_, outvars, k);
    }
  }
}
