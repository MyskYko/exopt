#pragma once

static inline int clog2(unsigned x) {
  static const unsigned t[5] = {0xffff0000,
                                0x0000ff00,
                                0x000000f0,
                                0x0000000c,
                                0x00000002};
  int y = (((x & (x - 1)) == 0)? 0: 1);
  int j = 16;
  for(int i = 0; i < 5; i++) {
    int k = (((x & t[i]) == 0)? 0: j);
    y += k;
    x >>= k;
    j >>= 1;
  }
  return y;
}

static inline int pow2roundup (int x) {
  if(x == 0) {
    return 0;
  }
  x--;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return x + 1;
}
