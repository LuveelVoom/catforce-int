// LifeAPI provide comfortable functions (API) to manipulate, iterate, evolve,
// compare and report Life objects. This is mainly done in order to provide fast
// (using C) but still comfortable search utility. Contributor Chris Cain.
// Written by Michael Simkin 2014

// Modified by pzq to support INT rules (albeit slower)
// June, 2023
// Rule is now selected at runtime via the "rule" line of the .in file.

#include <algorithm>
#include <array>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <random>
#include <string>
#include <chrono>
#include <cstdio>
#include <ctype.h>

// Apple silicon / ARM: use sse2neon instead of immintrin.h (x86 keeps immintrin.h)
#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
#include "sse2neon.h"
#else
#include <immintrin.h>
#endif

// Torus size (N x N). The default build is 64; `make CatForce128` builds a
// 128x128 version (keyword "torus 128" in the .in file switches to it).
#ifndef CATFORCE_N
#define CATFORCE_N 64
#endif
#define N CATFORCE_N

#if N == 64
typedef uint64_t row_t;
#elif N == 128
typedef unsigned __int128 row_t;
#else
#error "CATFORCE_N must be 64 or 128"
#endif

#define SUCCESS 1
#define FAIL 0

// GCC
#ifdef __GNUC__
#ifndef __clang__
// #include <x86intrin.h>

uint64_t reverse_uint64_t(uint64_t x) {
  const uint64_t h1 = 0x5555555555555555ULL;
  const uint64_t h2 = 0x3333333333333333ULL;
  const uint64_t h4 = 0x0F0F0F0F0F0F0F0FULL;
  const uint64_t v1 = 0x00FF00FF00FF00FFULL;
  const uint64_t v2 = 0x0000FFFF0000FFFFULL;
  x = ((x >> 1) & h1) | ((x & h1) << 1);
  x = ((x >> 2) & h2) | ((x & h2) << 2);
  x = ((x >> 4) & h4) | ((x & h4) << 4);
  x = ((x >> 8) & v1) | ((x & v1) << 8);
  x = ((x >> 16) & v2) | ((x & v2) << 16);
  x = (x >> 32) | (x << 32);
  return x;
}

#define __builtin_rotateleft64(x, r) (((x) << ((r) & 63)) | ((x) >> ((64 - ((r) & 63)) & 63)))
#define __builtin_rotateright64(x, r) (((x) >> ((r) & 63)) | ((x) << ((64 - ((r) & 63)) & 63)))
#define __builtin_bitreverse64 ::reverse_uint64_t
#endif
#endif

// MSVC
#ifdef __MSC_VER
#include <intrin.h>
#define __builtin_popcount __popcnt64
#define __builtin_rotateleft64 _rotl64
#define __builtin_rotateright64 _rotr64

inline int __builtin_ctzll(uint64_t x) {
  unsigned long log2;
  _BitScanReverse64(&log2, x);
  return log2;
}
#endif

// Row helpers: a row (one column of the torus) is N bits.
#if N == 64
inline row_t RowRotl(row_t x, unsigned k) { return __builtin_rotateleft64(x, k); }
inline row_t RowRotr(row_t x, unsigned k) { return __builtin_rotateright64(x, k); }
inline int RowPop(row_t x) { return __builtin_popcountll(x); }
inline int RowCtz(row_t x) { return __builtin_ctzll(x); }
inline int RowClz(row_t x) { return __builtin_clzll(x); }
inline row_t RowReverse(row_t x) { return __builtin_bitreverse64(x); }
inline uint64_t RowFold(row_t x) { return x; }
#else
inline row_t RowRotl(row_t x, unsigned k) {
  k &= N - 1;
  return k ? (x << k) | (x >> (N - k)) : x;
}
inline row_t RowRotr(row_t x, unsigned k) {
  k &= N - 1;
  return k ? (x >> k) | (x << (N - k)) : x;
}
inline int RowPop(row_t x) { return __builtin_popcountll((uint64_t)x) + __builtin_popcountll((uint64_t)(x >> 64)); }
inline int RowCtz(row_t x) { // x != 0
  uint64_t lo = (uint64_t)x;
  return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((uint64_t)(x >> 64));
}
inline int RowClz(row_t x) { // x != 0
  uint64_t hi = (uint64_t)(x >> 64);
  return hi ? __builtin_clzll(hi) : 64 + __builtin_clzll((uint64_t)x);
}
inline row_t RowReverse(row_t x) {
  return ((row_t)__builtin_bitreverse64((uint64_t)x) << 64) | __builtin_bitreverse64((uint64_t)(x >> 64));
}
inline uint64_t RowFold(row_t x) { return (uint64_t)x ^ ((uint64_t)(x >> 64) * 0x9E3779B97F4A7C15ULL); }
#endif
static const row_t ROW_ALL = ~(row_t)0;
// bit 4k set for every k (the 4-row groups of the table stepper)
static const row_t ROW_EVERY4 = (row_t)0x1111111111111111ULL | ((row_t)0x1111111111111111ULL << (N - 64));

inline row_t convolve_row(row_t x, row_t y) {
  if(y == 0)
    return 0;

  row_t result = 0;
  while (x != 0) {
    int lsb = RowCtz(x);
    result |= RowRotl(y, lsb);
    x &= ~(((row_t)1) << lsb);
  }
  return result;
}

namespace PRNG {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<uint64_t> dist(std::llround(std::pow(2,61)), std::llround(std::pow(2,62)));
// Public domain PRNG by Sebastian Vigna 2014, see http://xorshift.di.unimi.it

uint64_t s[16] = {0x12345678};
int p = 0;

uint64_t rand64() {
  uint64_t s0 = s[p];
  uint64_t s1 = s[p = (p + 1) & 15];
  s1 ^= s1 << 31; // a
  s1 ^= s1 >> 11; // b
  s0 ^= s0 >> 30; // c
  return (s[p] = s0 ^ s1) * 1181783497276652981ULL;
}

} // namespace PRNG

// void fastMemcpy(void *pvDest, void *pvSrc, size_t nBytes) {
//   assert(nBytes % 32 == 0);
//   assert((intptr_t(pvDest) & 31) == 0);
//   assert((intptr_t(pvSrc) & 31) == 0);
//   const __m256i *pSrc = reinterpret_cast<const __m256i*>(pvSrc);
//   __m256i *pDest = reinterpret_cast<__m256i*>(pvDest);
//   int64_t nVects = nBytes / sizeof(*pSrc);
//   for (; nVects > 0; nVects--, pSrc++, pDest++) {
//     const __m256i loaded = _mm256_stream_load_si256(pSrc);
//     _mm256_stream_si256(pDest, loaded);
//   }
//   _mm_sfence();
// }

enum CopyType { COPY, OR, XOR, AND, ANDNOT, ORNOT };

enum SymmetryTransform {
  Identity,
  ReflectAcrossXEven,
  ReflectAcrossX,
  ReflectAcrossYEven,
  ReflectAcrossY,
  Rotate90Even,
  Rotate90,
  Rotate270Even,
  Rotate270,
  Rotate180OddBoth,
  Rotate180EvenHorizontal,
  Rotate180EvenVertical,
  Rotate180EvenBoth,
  ReflectAcrossYeqX,
  ReflectAcrossYeqNegX,
  // reflect across y = -x+3/2, fixing (0,0), instead of y=-x+1/2,
  // sending (0,0) to (-1,-1). Needed for D4x_1 symmetry.
  ReflectAcrossYeqNegXP1
};

enum StaticSymmetry {
  C1,
  D2AcrossX,
  D2AcrossXEven,
  D2AcrossY,
  D2AcrossYEven,
  D2negdiagodd,
  D2diagodd,
  C2,
  C2even,
  C2verticaleven,
  C2horizontaleven,
  C4,
  C4even,
  D4,
  D4even,
  D4verticaleven,
  D4horizontaleven,
  D4diag,
  D4diageven,
  D8,
  D8even,
};

inline row_t RotateLeft(row_t x, unsigned int k) { return RowRotl(x, k); }
inline row_t RotateRight(row_t x, unsigned int k) { return RowRotr(x, k); }

inline row_t RotateLeft(row_t x) { return RotateLeft(x, 1); }
inline row_t RotateRight(row_t x) { return RotateRight(x, 1); }

class LifeTarget;

class LifeState {
public:
  row_t state[N];

  int min;
  int max;
  int gen;

  LifeState() : state{0}, min(0), max(N - 1), gen(0) {}

  void Set(int x, int y) { state[x] |= ((row_t)1 << (y)); }
  void Erase(int x, int y) { state[x] &= ~((row_t)1 << (y)); }
  int Get(int x, int y) const { return (int)((state[x] >> y) & 1); }
  void SetCell(int x, int y, int val) {
    if (val == 1) {
      Set((x + N) % N, (y + N) % N);
    }
    if (val == 0)
      Erase((x + N) % N, (y + N) % N);
  }
  int GetCell(int x, int y) const {
    return Get((x + N) % N, (y + N) % N);
  }
  uint64_t GetHash() const {
    row_t result = 0;

    for (int i = 0; i < N; i++) {
      result += RotateLeft(state[i], (int)(i / 2));
    }

    return RowFold(result);
  }

#if defined(__AVX2__) && N == 64
  void RecalculateMinMax() {
    min = 0;
    max = N - 1;

    const char *p = (const char *)state;
    size_t len = 8*N;
    const char *p_init = p;
    const char *endp = p + len;
    do {
      __m256i v1 = _mm256_loadu_si256((const __m256i*)p);
      __m256i v2 = _mm256_loadu_si256((const __m256i*)(p+32));
      __m256i vor = _mm256_or_si256(v1,v2);
      if (!_mm256_testz_si256(vor, vor)) {
        min = (p-p_init)/8;
        break;
      }
      p += 64;
    } while(p < endp);

    p = endp-64;
    do {
      __m256i v1 = _mm256_loadu_si256((const __m256i*)p);
      __m256i v2 = _mm256_loadu_si256((const __m256i*)(p+32));
      __m256i vor = _mm256_or_si256(v1,v2);
      if (!_mm256_testz_si256(vor, vor)) {
        max = (p-p_init)/8 + 7;
        break;
      }
      p -= 64;
    } while(p >= p_init);
  }
#else
  void RecalculateMinMax() {
    min = 0;
    max = N - 1;

    for (int i = 0; i < N; i++) {
      if (state[i] != 0) {
        min = i;
        break;
      }
    }

    for (int i = N - 1; i >= 0; i--) {
      if (state[i] != 0) {
        max = i;
        break;
      }
    }
  }
#endif

public:
  void Print() const;

  void Copy(const LifeState &delta, CopyType op) {
    if (op == COPY) {
      for (int i = 0; i < N; i++)
        state[i] = delta.state[i];

      min = delta.min;
      max = delta.max;
      gen = delta.gen;
      return;
    }
    if (op == OR) {
      for (int i = 0; i < N; i++)
        state[i] |= delta.state[i];
      min = std::min(min, delta.min);
      max = std::max(max, delta.max);
      return;
    }
    if (op == AND) {
      for (int i = 0; i < N; i++)
        state[i] &= delta.state[i];
    }
    if (op == ANDNOT) {
      for (int i = 0; i < N; i++)
        state[i] &= ~delta.state[i];
    }
    if (op == ORNOT) {
      for (int i = 0; i < N; i++)
        state[i] |= ~delta.state[i];
      RecalculateMinMax();
    }
    if (op == XOR) {
      for (int i = 0; i < N; i++)
        state[i] ^= delta.state[i];
      RecalculateMinMax();
    }
  }

  void Copy(const LifeState &delta) { Copy(delta, COPY); }

  inline void Copy(const LifeState &delta, int x, int y) {
    row_t temp1[N] = {0};

    if (x < 0)
      x += N;
    if (y < 0)
      y += N;

    for (int i = delta.min; i <= delta.max; i++)
      temp1[i] = RotateLeft(delta.state[i], y);

    memmove(state, temp1 + (N - x), x * sizeof(row_t));
    memmove(state + x, temp1, (N - x) * sizeof(row_t));

    min = 0;
    max = N - 1;
  }

  void Join(const LifeState &delta) { Copy(delta, OR); }

  void Join(const LifeState &delta, int x, int y) {
    row_t temp[2*N] = {0};

    if (x < 0)
      x += N;
    if (y < 0)
      y += N;

    for (int i = 0; i < N; i++) {
      temp[i]   = RotateLeft(delta.state[i], y);
      temp[i+N] = RotateLeft(delta.state[i], y);
    }

    const int shift = N - x;
    for (int i = 0; i < N; i++) {
      state[i] |= temp[i+shift];
    }

    min = 0;
    max = N - 1;
  }

  void JoinWSymChain(const LifeState &state, int x, int y,
                     const std::vector<SymmetryTransform> &symChain) {
    // instead of passing in the symmetry group {id, g_1, g_2,...g_n} and
    // applying each to default orientation we pass in a "chain" of symmetries
    // {h_1, ...h_n-1} that give the group when "chained together": g_j =
    // product of h_1 thru h_j that way, we don't need to initialize a new
    // LifeState for each symmetry.

    LifeState transformed = state;
    transformed.Move(x, y);

    for (auto sym : symChain) {
      LifeState soFar = transformed;
      soFar.Transform(sym);
      transformed.Join(soFar);
    }
    Join(transformed);
  }

  void JoinWSymChain(const LifeState &state,
                     const std::vector<SymmetryTransform> &symChain) {
    LifeState transformed = state;

    for (auto sym : symChain) {
      LifeState soFar = transformed;
      soFar.Transform(sym);
      transformed.Join(soFar);
    }
    Join(transformed);
  }

  unsigned GetPop() const {
    unsigned pop = 0;

    for (int i = min; i <= max; i++) {
      pop += RowPop(state[i]);
    }

    return pop;
  }

  // bool IsEmpty() const {
  //   for (int i = 0; i < N; i++) {
  //     if(state[i] != 0)
  //       return false;
  //   }

  //   return true;
  // }

  bool IsEmpty() const {
    row_t all = 0;
    for (int i = 0; i < N; i++) {
      all |= state[i];
    }

    return all == 0;
  }

  void Inverse() {
    for (int i = 0; i < N; i++) {
      state[i] = ~state[i];
    }
  }

  bool operator==(const LifeState &b) const {
    for (int i = 0; i < N; i++)
      if (state[i] != b.state[i])
        return false;

    return true;
  }

  LifeState operator~() const {
    LifeState result;
    for (int i = 0; i < N; i++) {
      result.state[i] = ~state[i];
    }
    return result;
  }

  LifeState operator&(const LifeState &other) const {
    LifeState result;
    for (int i = 0; i < N; i++) {
      result.state[i] = state[i] & other.state[i];
    }
    return result;
  }

  LifeState& operator&=(const LifeState &other) {
    for (int i = 0; i < N; i++) {
      state[i] = state[i] & other.state[i];
    }
    return *this;
  }

  LifeState operator|(const LifeState &other) const {
    LifeState result;
    for (int i = 0; i < N; i++) {
      result.state[i] = state[i] | other.state[i];
    }
    return result;
  }

  LifeState& operator|=(const LifeState &other) {
    for (int i = 0; i < N; i++) {
      state[i] = state[i] | other.state[i];
    }
    return *this;
  }

  LifeState operator^(const LifeState &other) const {
    LifeState result;
    for (int i = 0; i < N; i++) {
      result.state[i] = state[i] ^ other.state[i];
    }
    return result;
  }

  LifeState& operator^=(const LifeState &other) {
    for (int i = 0; i < N; i++) {
      state[i] = state[i] ^ other.state[i];
    }
    return *this;
  }

  inline bool AreDisjoint(const LifeState &pat) const {
    int min = 0;
    int max = N - 1;

    row_t differences = 0;
    #pragma clang loop vectorize(enable)
    for (int i = min; i <= max; i++) {
      row_t difference = (~state[i] & pat.state[i]) ^ (pat.state[i]);
      differences |= difference;
    }

    return differences == 0;
  }

  inline bool Contains(const LifeState &pat) const {
    int min = 0;
    int max = N - 1;

    row_t differences = 0;
    #pragma clang loop vectorize(enable)
    for (int i = min; i <= max; i++) {
      row_t difference = (state[i] & pat.state[i]) ^ (pat.state[i]);
      differences |= difference;
    }

    return differences == 0;
  }

  bool Contains(const LifeState &pat, int targetDx, int targetDy) const {
    int min = pat.min;
    int max = pat.max;

    int dy = (targetDy + N) % N;

    for (int i = min; i <= max; i++) {
      int curX = (N + i + targetDx) % N;

      if ((RotateRight(state[curX], dy) & pat.state[i]) != (pat.state[i]))
        return false;
    }
    return true;
  }

  bool AreDisjoint(const LifeState &pat, int targetDx, int targetDy) const {
    int min = pat.min;
    int max = pat.max;
    int dy = (targetDy + N) % N;

    for (int i = min; i <= max; i++) {
      int curX = (N + i + targetDx) % N;

      if (((~RotateRight(state[curX], dy)) & pat.state[i]) != pat.state[i])
        return false;
    }

    return true;
  }

  inline bool Contains(const LifeTarget &target, int dx, int dy) const;
  inline bool Contains(const LifeTarget &target) const;

  void Reverse(int idxS, int idxE) {
    for (int i = 0; idxS + 2*i < idxE; i++) {
      int l = idxS + i;
      int r = idxE - i;

      row_t temp = state[l];
      state[l] = state[r];
      state[r] = temp;
    }
  }

  void Move(int x, int y) {
    row_t temp[2*N] = {0};

    if (x < 0)
      x += N;
    if (y < 0)
      y += N;

    for (int i = 0; i < N; i++) {
      temp[i]   = RotateLeft(state[i], y);
      temp[i+N] = RotateLeft(state[i], y);
    }

    const int shift = N - x;
    for (int i = 0; i < N; i++) {
      state[i] = temp[i+shift];
    }

    if ((min + x) % N < (max + x) % N) {
      min = (min + x) % N;
      max = (max + x) % N;
    } else {
      min = 0;
      max = N - 1;
    }
  }

  void BitReverse() {
    for (int i = 0; i < N; i++) {
      state[i] = RowReverse(state[i]);
    }
  }

  void FlipY() { // even reflection across y-axis, ie (0,0) maps to (0, -1)
    Reverse(0, N - 1);
  }

  void Transpose(bool whichDiagonal) {
    int j, k;
    row_t m, t;

    for (j = N / 2, m = ROW_ALL >> (N / 2); j; j >>= 1, m ^= m << j) {
      for (k = 0; k < N; k = ((k | j) + 1) & ~j) {
        if (whichDiagonal) {
          t = (state[k] ^ (state[k | j] >> j)) & m;
          state[k] ^= t;
          state[k | j] ^= (t << j);
        } else {
          t = (state[k] >> j ^ (state[k | j])) & m;
          state[k] ^= (t << j);
          state[k | j] ^= t;
        }
      }
    }
    min = 0;
    max = N - 1;
  }

  void Transpose() { Transpose(true); }

  // even reflection across x-axis, ie (0,0) maps to (0, -1)
  void FlipX() { BitReverse(); }

  void Transform(SymmetryTransform transf);

  void Transform(int dx, int dy, SymmetryTransform transf) {
    Move(dx, dy);
    Transform(transf);
  }

  LifeState ZOI() const {
    LifeState temp;
    LifeState boundary;
    for (int i = 0; i < N; i++) {
      row_t col = state[i];
      temp.state[i] = col | RotateLeft(col) | RotateRight(col);
    }

    boundary.state[0] = temp.state[N - 1] | temp.state[0] | temp.state[1];

    for (int i = 1; i < N - 1; i++)
      boundary.state[i] = temp.state[i - 1] | temp.state[i] | temp.state[i + 1];

    boundary.state[N - 1] = temp.state[N - 2] | temp.state[N - 1] | temp.state[0];

    boundary.RecalculateMinMax();
    return boundary;
  }

  LifeState GetBoundary() const {
    LifeState boundary = ZOI();
    boundary.Copy(*this, ANDNOT);
    return boundary;
  }

  LifeState BigZOI() const {
    LifeState b;
    b.state[0] = state[0] | RotateLeft(state[0]) | RotateRight(state[0]) |
                 state[N - 1] | state[0 + 1];
    for (int i = 1; i < N-1; i++) {
      b.state[i] = state[i] | RotateLeft(state[i]) | RotateRight(state[i]) | state[i-1] | state[i+1];
    }
    b.state[N-1] = state[N-1] | RotateLeft(state[N-1]) | RotateRight(state[N-1]) |
                 state[N-1 - 1] | state[0];

    LifeState c;
    c.state[0] = b.state[0] | b.state[N - 1] | b.state[0 + 1];
    for (int i = 1; i < N - 1; i++) {
      c.state[i] = b.state[i] | b.state[i - 1] | b.state[i + 1];
    }
    c.state[N - 1] = b.state[N - 1] | b.state[N - 1 - 1] | b.state[0];

    LifeState zoi;

    zoi.state[0] =
      c.state[0] | RotateLeft(c.state[0]) | RotateRight(c.state[0]);
    for (int i = 1; i < N - 1; i++) {
      zoi.state[i] =
        c.state[i] | RotateLeft(c.state[i]) | RotateRight(c.state[i]);
    }
    zoi.state[N - 1] = c.state[N - 1] | RotateLeft(c.state[N - 1]) |
      RotateRight(c.state[N - 1]);

    zoi.RecalculateMinMax();
    return zoi;
  }

  // OR of y rotated by 0, 1, ..., len-1 (= convolve_row(run of len ones, y))
  static inline row_t Smear(row_t y, unsigned len) {
    row_t acc = y;
    for (unsigned covered = 1; covered < len;) {
      unsigned step = std::min(covered, len - covered);
      acc |= RowRotl(acc, step);
      covered += step;
    }
    return acc;
  }

  // Every cell of `other` shifted by every live cell of this (cell (j, b) of
  // this places other's cell (c, r) at (c + j, r + b)). Only the non-empty
  // columns of `other` are visited, and each run of live cells in a column of
  // this is applied at once.
  LifeState Convolve(const LifeState &other) const {
    LifeState result;
    int cols[N];
    int nc = 0;
    for (int c = 0; c < N; c++)
      if (other.state[c] != 0)
        cols[nc++] = c;
    if (nc == 0)
      return result;

    for (unsigned j = 0; j < N; j++) {
      row_t x = state[j];
      if (x == 0)
        continue;
      if (x == ROW_ALL) {
        for (int t = 0; t < nc; t++)
          result.state[(cols[t] + j) & (N - 1)] = ROW_ALL;
        continue;
      }
      while (x != 0) {
        unsigned postshift;
        row_t shifted;
        if ((x & 1) == 0) { // possibly wrapped
          int lsb = RowCtz(x);
          shifted = RowRotr(x, lsb);
          postshift = lsb;
        } else {
          int lead = RowClz(~x);
          shifted = RowRotl(x, lead);
          postshift = N - lead;
        }
        unsigned runlength = RowCtz(~shifted);
        row_t run = ((row_t)1 << runlength) - 1;
        for (int t = 0; t < nc; t++) {
          const int c = cols[t];
          result.state[(c + j) & (N - 1)] |= RowRotl(Smear(other.state[c], runlength), postshift);
        }
        x &= ~RowRotl(run, postshift);
      }
    }

    result.min = 0;
    result.max = N - 1;
    return result;
  }

  void Clear() {
    for (int i = 0; i < N; i++)
      state[i] = 0;

    min = 0;
    max = 0;
    gen = 0;
  }

  LifeState MatchLive(const LifeState &live) const {
    LifeState invThis = ~*this;
    LifeState flipLive = live;
    flipLive.Transform(Rotate180OddBoth);
    return ~invThis.Convolve(flipLive);
  }

  LifeState MatchLiveAndDead(const LifeState &live, const LifeState &dead) const {
    LifeState invThis = ~*this;
    LifeState flipLive = live;
    flipLive.Transform(Rotate180OddBoth);
    LifeState flipDead = dead;
    flipDead.Transform(Rotate180OddBoth);
    return ~invThis.Convolve(flipLive) & ~Convolve(flipDead);
  }

  LifeState Match(const LifeState &live) const {
    return MatchLiveAndDead(live, live.GetBoundary());
  }

  LifeState Match(const LifeTarget &target) const;

private:
  void inline Add(row_t &b1, row_t &b0, const row_t &val) {
    b1 |= b0 & val;
    b0 ^= val;
  }

  void inline Add(row_t &b2, row_t &b1, row_t &b0,
                  const row_t &val) {
    row_t t_b2 = b0 & val;

    b2 |= t_b2 & b1;
    b1 ^= t_b2;
    b0 ^= val;
  }

  row_t inline Evolve(const row_t &temp, const row_t &bU0,
                         const row_t &bU1, const row_t &bB0,
                         const row_t &bB1) {
    row_t sum0 = RotateLeft(temp);

    row_t sum1 = 0;
    Add(sum1, sum0, RotateRight(temp));
    Add(sum1, sum0, bU0);

    row_t sum2 = 0;
    Add(sum2, sum1, bU1);
    Add(sum2, sum1, sum0, bB0);
    Add(sum2, sum1, bB1);

    return ~sum2 & sum1 & (temp | sum0);
  }

  // From Page 15 of
  // https://www.gathering4gardner.org/g4g13gift/math/RokickiTomas-GiftExchange-LifeAlgorithms-G4G13.pdf
  row_t inline Rokicki(const row_t &a, const row_t &bU0,
                          const row_t &bU1, const row_t &bB0,
                          const row_t &bB1) const {
    row_t aw = RotateLeft(a);
    row_t ae = RotateRight(a);
    row_t s0 = aw ^ ae;
    row_t s1 = aw & ae;
    row_t ts0 = bB0 ^ bU0;
    row_t ts1 = (bB0 & bU0) | (ts0 & s0);
    return (bB1 ^ bU1 ^ ts1 ^ s1) & ((bB1 | bU1) ^ (ts1 | s1)) &
           ((ts0 ^ s0) | a);
  }

public:
  void Step();
  // Next generation of only the columns selected by `cols` (bit i = column i),
  // written to out[i]; other entries of out are not touched. Equal to Step()'s
  // result for those columns, but costs only as much as the columns requested.
  void StepColumnsInto(row_t cols, row_t *out) const;

  void Step(int numIters) {
    for (int i = 0; i < numIters; i++) {
      Step();
      // RemoveGliders();
    }
  }

  static int Parse(LifeState &state, const char *rle, int starti);

  static int Parse(LifeState &state, const char *rle, int dx, int dy) {
    if (Parse(state, rle, 0) == -1) {
      state.Move(dx, dy);
      return SUCCESS;
    } else {
      return FAIL;
    }
  }

  static int Parse(LifeState &state, const char *rle) {
    return Parse(state, rle, 0, 0);
  }

  static int Parse(LifeState &state, const char *rle, int dx, int dy,
                   SymmetryTransform transf) {
    int result = Parse(state, rle);

    if (result == SUCCESS)
      state.Transform(dx, dy, transf);

    return result;
  }

  static LifeState Parse(const char *rle, int dx, int dy,
                         SymmetryTransform trans) {
    LifeState result;
    Parse(result, rle);
    result.Transform(dx, dy, trans);
    result.RecalculateMinMax();

    return result;
  }

  static LifeState Parse(const char *rle, int dx, int dy) {
    LifeState result;
    Parse(result, rle, dx, dy);

    return result;
  }

  static LifeState Parse(const char *rle) { return Parse(rle, 0, 0); }

  std::string RLE() const;

  static LifeState RandomState() {
    LifeState result;
    for (int i = 0; i < N; i++)
      for (int b = 0; b < N; b += 64)
        result.state[i] |= (row_t)PRNG::dist(PRNG::e2) << b;

    result.RecalculateMinMax();

    return result;
  }

#if defined(__AVX2__) && N == 64
  // https://stackoverflow.com/questions/56153183/is-using-avx2-can-implement-a-faster-processing-of-lzcnt-on-a-word-array
  std::pair<int, int> FirstOn() const
  {
    const char *p = (const char *)state;
    size_t len = 8*N;
    //assert(len % 64 == 0);
    //optimal if p is 64-byte aligned, so we're checking single cache-lines
    const char *p_init = p;
    const char *endp = p + len;
    do {
      __m256i v1 = _mm256_loadu_si256((const __m256i*)p);
      __m256i v2 = _mm256_loadu_si256((const __m256i*)(p+32));
      __m256i vor = _mm256_or_si256(v1,v2);
      if (!_mm256_testz_si256(vor, vor)) {        // find the first non-zero cache line
        __m256i v1z = _mm256_cmpeq_epi32(v1, _mm256_setzero_si256());
        __m256i v2z = _mm256_cmpeq_epi32(v2, _mm256_setzero_si256());
        uint32_t zero_map = _mm256_movemask_ps(_mm256_castsi256_ps(v1z));
        zero_map |= _mm256_movemask_ps(_mm256_castsi256_ps(v2z)) << 8;

        unsigned idx = __builtin_ctz(~zero_map);  // Use ctzll for GCC, because GCC is dumb and won't optimize away a movsx
        uint32_t nonzero_chunk;
        memcpy(&nonzero_chunk, p+4*idx, sizeof(nonzero_chunk));  // aliasing / alignment-safe load
        if(idx % 2 == 0) {
          return std::make_pair((p-p_init + 4*idx)/8, __builtin_ctz(nonzero_chunk));
        } else {
          return std::make_pair((p-p_init + 4*(idx-1))/8, __builtin_ctz(nonzero_chunk) + 32);
        }
      }
      p += 64;
    } while(p < endp);
    return std::make_pair(-1, -1);
  }
#else
  std::pair<int, int> FirstOn() const {
    int foundq = N;
    for (int x = 0; x < N; x+=4) {
      if (state[x] != 0 ||
          state[x+1] != 0 ||
          state[x+2] != 0 ||
          state[x+3] != 0) {
        foundq = x;
      }
    }
    if (foundq == N) {
      return std::make_pair(0, 0);
    }

    int foundx = foundq + 3;
    if (state[foundq] != 0) {
      foundx = foundq;
    } else if (state[foundq + 1] != 0) {
      foundx = foundq + 1;
    } else if (state[foundq + 2] != 0) {
      foundx = foundq + 2;
    }
    return std::make_pair(foundx, RowCtz(state[foundx]));
  }
#endif

  static LifeState SolidRect(int x, int y, int w, int h) {
    row_t column;
    if (h < N)
      column = RotateLeft(((row_t)1 << h) - 1, y);
    else
      column = ROW_ALL;

    unsigned start, end;
    if (w < N) {
      start = (x + N) % N;
      end = (x + w + N) % N;
    } else {
      start = 0;
      end = N;
    }

    LifeState result;
    if (end > start) {
      for (unsigned int i = start; i < end; i++)
        result.state[i] = column;
    } else {
      for (unsigned int i = 0; i < end; i++)
        result.state[i] = column;
      for (unsigned int i = start; i < N; i++)
        result.state[i] = column;
    }
    return result;
  }

  std::array<int, 4> XYBounds() const {
    int minCol = -N / 2;
    int maxCol = N / 2 - 1;

    for (int i = -N / 2; i <= N / 2 - 1; i++) {
      if (state[(i + N) % N] != 0) {
        minCol = i;
        break;
      }
    }

    for (int i = N / 2 - 1; i >= -N / 2; i--) {
      if (state[(i + N) % N] != 0) {
        maxCol = i;
        break;
      }
    }

    row_t orOfCols(0);
    for (int i = minCol; i <= maxCol; ++i) {
      orOfCols = orOfCols | state[(i + N) % N];
    }
    if (orOfCols == 0) {
      return std::array<int, 4>({0, 0, 0, 0});
    }
    orOfCols = RowRotr(orOfCols, N / 2);
    int topMargin = RowCtz(orOfCols);
    int bottomMargin = RowClz(orOfCols);
    return std::array<int, 4>(
        {minCol, topMargin - N / 2, maxCol, N / 2 - 1 - bottomMargin});
  }
};

// ===========================================================================
// Rule support
//
// The rule is set at runtime (see SetRule) from the "rule" line of the .in
// file. Accepted formats:
//   * Life-like / INT (Hensel) strings:  B3/S23, B2ei3aeikqr/S2-a3, 23/3, ...
//   * Golly MAP strings:                 MAPARYXfhZofugWaH7oaIDogBZofuhog...
// B3/S23 (Conway) automatically uses the fast bit-parallel stepper; every
// other rule uses a table-driven stepper that only visits cells near live
// cells. Rules with B0 are not supported (infinite background would flip).
// ===========================================================================

// Table index (matches the Step loop below): bit 3*(1-dx) + (dy+1) is the
// state of the cell at (x+dx, y+dy); bit 4 is the cell itself. So bits 6..8 are
// the west column (x-1) and bits 0..2 the east column (x+1). Getting this
// wrong mirrors non-isotropic (MAP) rules left-right.
__attribute__((aligned(64))) static bool table[512];
// 3-column x 6-row window table: maps (6 rows of left, centre, right column)
// to the next state of the 4 middle rows, so one lookup updates 4 cells.
// Rebuilt by SetRule.
static uint8_t winTable[1u << 18];

// Rule-specific compiled kernels (see RuleCompiler.h and --compile-rule). When
// the current rule has one built into this binary and it passes its self-check,
// Step() uses it instead of the table stepper.
typedef void (*CompiledStepFn)(row_t *state);
typedef void (*CompiledColsFn)(const row_t *state, row_t cols, row_t *out);
static CompiledStepFn compiledStep = nullptr;
static CompiledColsFn compiledCols = nullptr;
static std::string compiledInfo;
static void SelectCompiledKernel();
static bool useFastLife = true;     // true iff the current rule is B3/S23
static std::string currentRuleName = "B3/S23";

// Hensel isotropic neighbourhood classes. Each mask is an 8-bit neighbourhood,
// MSB first, row-major:  NW=0x80 N=0x40 NE=0x20 W=0x10 E=0x08 SW=0x04 S=0x02 SE=0x01
// (class definitions taken from the ca-rules crate, which follows Golly/Hensel)
static const struct { char n; char l; uint8_t masks[8]; int cnt; } HENSEL_CLASSES[] = {
  {'0','c',{0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},1},
  {'1','c',{0x01,0x04,0x20,0x80,0x00,0x00,0x00,0x00},4},
  {'1','e',{0x02,0x08,0x10,0x40,0x00,0x00,0x00,0x00},4},
  {'2','c',{0x05,0x21,0x84,0xa0,0x00,0x00,0x00,0x00},4},
  {'2','e',{0x0a,0x12,0x48,0x50,0x00,0x00,0x00,0x00},4},
  {'2','k',{0x0c,0x11,0x22,0x30,0x41,0x44,0x82,0x88},8},
  {'2','a',{0x03,0x06,0x09,0x14,0x28,0x60,0x90,0xc0},8},
  {'2','i',{0x18,0x42,0x00,0x00,0x00,0x00,0x00,0x00},2},
  {'2','n',{0x24,0x81,0x00,0x00,0x00,0x00,0x00,0x00},2},
  {'3','c',{0x25,0x85,0xa1,0xa4,0x00,0x00,0x00,0x00},4},
  {'3','e',{0x1a,0x4a,0x52,0x58,0x00,0x00,0x00,0x00},4},
  {'3','k',{0x32,0x4c,0x51,0x8a,0x00,0x00,0x00,0x00},4},
  {'3','a',{0x0b,0x16,0x68,0xd0,0x00,0x00,0x00,0x00},4},
  {'3','i',{0x07,0x29,0x94,0xe0,0x00,0x00,0x00,0x00},4},
  {'3','n',{0x0d,0x15,0x23,0x61,0x86,0xa8,0xb0,0xc4},8},
  {'3','y',{0x31,0x45,0x8c,0xa2,0x00,0x00,0x00,0x00},4},
  {'3','q',{0x26,0x2c,0x34,0x64,0x83,0x89,0x91,0xc1},8},
  {'3','j',{0x0e,0x13,0x2a,0x49,0x54,0x70,0x92,0xc8},8},
  {'3','r',{0x19,0x1c,0x38,0x43,0x46,0x62,0x98,0xc2},8},
  {'4','c',{0xa5,0x00,0x00,0x00,0x00,0x00,0x00,0x00},1},
  {'4','e',{0x5a,0x00,0x00,0x00,0x00,0x00,0x00,0x00},1},
  {'4','k',{0x33,0x4d,0x55,0x71,0x8e,0xaa,0xb2,0xcc},8},
  {'4','a',{0x0f,0x17,0x2b,0x69,0x96,0xd4,0xe8,0xf0},8},
  {'4','i',{0x1d,0x63,0xb8,0xc6,0x00,0x00,0x00,0x00},4},
  {'4','n',{0x27,0x2d,0x87,0x95,0xa9,0xb4,0xe1,0xe4},8},
  {'4','y',{0x35,0x65,0x8d,0xa3,0xa6,0xac,0xb1,0xc5},8},
  {'4','q',{0x36,0x6c,0x8b,0xd1,0x00,0x00,0x00,0x00},4},
  {'4','j',{0x3a,0x4e,0x53,0x59,0x5c,0x72,0x9a,0xca},8},
  {'4','r',{0x1b,0x1e,0x4b,0x56,0x6a,0x78,0xd2,0xd8},8},
  {'4','t',{0x39,0x47,0x9c,0xe2,0x00,0x00,0x00,0x00},4},
  {'4','w',{0x2e,0x74,0x93,0xc9,0x00,0x00,0x00,0x00},4},
  {'4','z',{0x3c,0x66,0x99,0xc3,0x00,0x00,0x00,0x00},4},
  {'5','c',{0x5b,0x5e,0x7a,0xda,0x00,0x00,0x00,0x00},4},
  {'5','e',{0xa7,0xad,0xb5,0xe5,0x00,0x00,0x00,0x00},4},
  {'5','k',{0x75,0xae,0xb3,0xcd,0x00,0x00,0x00,0x00},4},
  {'5','a',{0x2f,0x97,0xe9,0xf4,0x00,0x00,0x00,0x00},4},
  {'5','i',{0x1f,0x6b,0xd6,0xf8,0x00,0x00,0x00,0x00},4},
  {'5','n',{0x3b,0x4f,0x57,0x79,0x9e,0xdc,0xea,0xf2},8},
  {'5','y',{0x5d,0x73,0xba,0xce,0x00,0x00,0x00,0x00},4},
  {'5','q',{0x3e,0x6e,0x76,0x7c,0x9b,0xcb,0xd3,0xd9},8},
  {'5','j',{0x37,0x6d,0x8f,0xab,0xb6,0xd5,0xec,0xf1},8},
  {'5','r',{0x3d,0x67,0x9d,0xb9,0xbc,0xc7,0xe3,0xe6},8},
  {'6','c',{0x5f,0x7b,0xde,0xfa,0x00,0x00,0x00,0x00},4},
  {'6','e',{0xaf,0xb7,0xed,0xf5,0x00,0x00,0x00,0x00},4},
  {'6','k',{0x77,0x7d,0xbb,0xbe,0xcf,0xdd,0xee,0xf3},8},
  {'6','a',{0x3f,0x6f,0x9f,0xd7,0xeb,0xf6,0xf9,0xfc},8},
  {'6','i',{0xbd,0xe7,0x00,0x00,0x00,0x00,0x00,0x00},2},
  {'6','n',{0x7e,0xdb,0x00,0x00,0x00,0x00,0x00,0x00},2},
  {'7','c',{0x7f,0xdf,0xfb,0xfe,0x00,0x00,0x00,0x00},4},
  {'7','e',{0xbf,0xef,0xf7,0xfd,0x00,0x00,0x00,0x00},4},
  {'8','c',{0xff,0x00,0x00,0x00,0x00,0x00,0x00,0x00},1},
};

static const int NUM_HENSEL_CLASSES = sizeof(HENSEL_CLASSES) / sizeof(HENSEL_CLASSES[0]);

// Convert a 3x3 table index (see above) to the 8-bit Hensel neighbourhood mask.
static inline int TableIndexToMask(int idx) {
  int mask = 0;
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      if (dx == 0 && dy == 0)
        continue;
      if (idx & (1 << (3 * (1 - dx) + (dy + 1)))) {
        int row = dy + 1, col = dx + 1;           // row 0 = top
        // position among the 8 neighbours in row-major order (centre skipped)
        int pos = (row == 0) ? col : (row == 1) ? (col == 0 ? 3 : 4) : 5 + col;
        int bit = 7 - pos;
        mask |= 1 << bit;
      }
    }
  }
  return mask;
}

// Parses a Hensel INT string like "B2ei3aeikqr/S2-a3". Fills birth[256] and
// survive[256] (indexed by 8-bit neighbourhood mask). Returns false + err on
// failure.
static bool ParseHenselRule(const std::string &rule, bool *birth,
                            bool *survive, std::string &err) {
  std::string r;
  for (char ch : rule)
    if (!isspace((unsigned char)ch))
      r += ch;
  if (r.empty()) {
    err = "empty rule";
    return false;
  }

  bool hasBS = false;
  for (char ch : r)
    if (ch == 'B' || ch == 'b' || ch == 'S' || ch == 's')
      hasBS = true;

  // 0 = birth, 1 = survival
  int mode = -1;
  size_t i = 0;
  if (!hasBS) {
    // numeric S/B notation, e.g. 23/3
    size_t slash = r.find('/');
    if (slash == std::string::npos) {
      err = "expected something like B3/S23";
      return false;
    }
    r = "S" + r.substr(0, slash) + "/B" + r.substr(slash + 1);
  }

  bool seenB = false, seenS = false;
  while (i < r.size()) {
    char ch = r[i];
    if (ch == '/') {
      i++;
      continue;
    }
    if (ch == 'B' || ch == 'b') {
      mode = 0;
      seenB = true;
      i++;
      continue;
    }
    if (ch == 'S' || ch == 's') {
      mode = 1;
      seenS = true;
      i++;
      continue;
    }
    if (mode < 0 || ch < '0' || ch > '8') {
      err = std::string("unexpected character '") + ch + "'";
      return false;
    }
    i++;
    bool negate = false;
    if (i < r.size() && r[i] == '-') {
      negate = true;
      i++;
    }
    std::string letters;
    while (i < r.size() && isalpha((unsigned char)r[i]) && r[i] != 'B' &&
           r[i] != 'b' && r[i] != 'S' && r[i] != 's') {
      letters += (char)tolower((unsigned char)r[i]);
      i++;
    }
    bool *dest = (mode == 0) ? birth : survive;

    // validate letters
    for (char l : letters) {
      bool ok = false;
      for (int c = 0; c < NUM_HENSEL_CLASSES; c++)
        if (HENSEL_CLASSES[c].n == ch && HENSEL_CLASSES[c].l == l)
          ok = true;
      if (!ok) {
        err = std::string("letter '") + l + "' is not valid for neighbour count " + ch;
        return false;
      }
    }
    if (negate && letters.empty()) {
      err = std::string("'-' after ") + ch + " must be followed by letters";
      return false;
    }

    for (int c = 0; c < NUM_HENSEL_CLASSES; c++) {
      if (HENSEL_CLASSES[c].n != ch)
        continue;
      bool listed = letters.find(HENSEL_CLASSES[c].l) != std::string::npos;
      bool on;
      if (letters.empty())
        on = true;               // bare digit: every configuration
      else if (negate)
        on = !listed;
      else
        on = listed;
      if (on)
        for (int m = 0; m < HENSEL_CLASSES[c].cnt; m++)
          dest[HENSEL_CLASSES[c].masks[m]] = true;
    }
  }
  if (!seenB && !seenS) {
    err = "rule must contain B and/or S";
    return false;
  }
  return true;
}

static int Base64Val(char c) {
  if (c >= 'A' && c <= 'Z') return c - 'A';
  if (c >= 'a' && c <= 'z') return c - 'a' + 26;
  if (c >= '0' && c <= '9') return c - '0' + 52;
  if (c == '+') return 62;
  if (c == '/') return 63;
  return -1;
}

// Parses a Golly MAP rule (512 bits, base64). Golly's bit index is
// 256*NW + 128*N + 64*NE + 32*W + 16*C + 8*E + 4*SW + 2*S + SE, bits are read
// MSB first from each decoded byte. Fills birth/survive indexed by 8-bit mask.
static bool ParseMapRule(const std::string &rule, bool *birth, bool *survive,
                         std::string &err) {
  std::string r = rule.substr(3);
  while (!r.empty() && (r.back() == '=' || isspace((unsigned char)r.back())))
    r.pop_back();
  std::vector<int> bits;
  for (char ch : r) {
    int v = Base64Val(ch);
    if (v < 0) {
      err = "invalid base64 character in MAP rule";
      return false;
    }
    for (int k = 5; k >= 0; k--)
      bits.push_back((v >> k) & 1);
  }
  if (bits.size() < 512) {
    err = "MAP rule is too short (need 512 bits / 86 base64 characters)";
    return false;
  }
  for (int k = 0; k < 512; k++) {
    if (!bits[k])
      continue;
    int mask = ((k & 0x1e0) >> 1) | (k & 0x0f);
    if (k & 0x10)
      survive[mask] = true;
    else
      birth[mask] = true;
  }
  return true;
}

// Sets the global rule. Returns false and fills err if the rule is invalid.
static bool SetRule(const std::string &ruleIn, std::string &err) {
  std::string rule = ruleIn;
  while (!rule.empty() && isspace((unsigned char)rule.back()))
    rule.pop_back();
  size_t st = 0;
  while (st < rule.size() && isspace((unsigned char)rule[st]))
    st++;
  rule = rule.substr(st);

  bool birth[256] = {false}, survive[256] = {false};
  bool ok;
  if (rule.size() > 3 && (rule.compare(0, 3, "MAP") == 0 ||
                          rule.compare(0, 3, "map") == 0))
    ok = ParseMapRule(rule, birth, survive, err);
  else
    ok = ParseHenselRule(rule, birth, survive, err);
  if (!ok)
    return false;

  if (birth[0]) {
    err = "rules with B0 are not supported";
    return false;
  }

  bool isLife = true;
  for (int idx = 0; idx < 512; idx++) {
    int mask = TableIndexToMask(idx);
    bool alive = (idx >> 4) & 1;
    table[idx] = alive ? survive[mask] : birth[mask];
    int cnt = __builtin_popcount(mask);
    bool lifeVal = alive ? (cnt == 2 || cnt == 3) : (cnt == 3);
    if (table[idx] != lifeVal)
      isLife = false;
  }
  useFastLife = isLife;
  currentRuleName = rule;
  for (unsigned w = 0; w < (1u << 18); w++) {
    unsigned a = w & 0x3F, b = (w >> 6) & 0x3F, c = (w >> 12) & 0x3F;
    unsigned o = 0;
    for (int r = 0; r < 4; r++) {
      unsigned neigh = (((a >> r) & 7) << 6) | (((b >> r) & 7) << 3) | ((c >> r) & 7);
      o |= (unsigned)table[neigh] << r;
    }
    winTable[w] = (uint8_t)o;
  }
  SelectCompiledKernel();
  return true;
}

// Renders the current table back into a canonical Hensel string (for logging).
static std::string CurrentRuleCanonical() {
  bool birth[256], survive[256];
  for (int m = 0; m < 256; m++) {
    birth[m] = survive[m] = false;
  }
  for (int idx = 0; idx < 512; idx++) {
    int mask = TableIndexToMask(idx);
    if ((idx >> 4) & 1)
      survive[mask] = table[idx];
    else
      birth[mask] = table[idx];
  }
  std::string out;
  for (int part = 0; part < 2; part++) {
    out += part == 0 ? "B" : "/S";
    const bool *arr = part == 0 ? birth : survive;
    for (int n = 0; n <= 8; n++) {
      std::string on, off;
      int total = 0, onCount = 0;
      for (int c = 0; c < NUM_HENSEL_CLASSES; c++) {
        if (HENSEL_CLASSES[c].n != '0' + n)
          continue;
        total++;
        if (arr[HENSEL_CLASSES[c].masks[0]]) {
          onCount++;
          on += HENSEL_CLASSES[c].l;
        } else
          off += HENSEL_CLASSES[c].l;
      }
      if (onCount == 0)
        continue;
      out += (char)('0' + n);
      if (onCount == total)
        continue;
      if (on.size() <= off.size())
        out += on;
      else
        out += "-" + off;
    }
  }
  return out;
}

// Rule string for RLE headers: the user's own string, except MAP strings,
// which are written out as a readable Hensel string.
// True if the current rule is invariant under rotations and reflections.
static bool RuleIsIsotropic() {
  if (useFastLife)
    return true;
  for (int idx = 0; idx < 512; idx++) {
    int rot = 0, ref = 0;
    for (int dx = -1; dx <= 1; dx++)
      for (int dy = -1; dy <= 1; dy++)
        if (idx & (1 << (3 * (dx + 1) + (dy + 1)))) {
          rot |= 1 << (3 * (-dy + 1) + (dx + 1));   // rotate 90 degrees
          ref |= 1 << (3 * (-dx + 1) + (dy + 1));   // mirror
        }
    if (table[idx] != table[rot] || table[idx] != table[ref])
      return false;
  }
  return true;
}

// Exact identity of the current rule (all 512 table entries), for cache keys.
static std::string RuleTableKey() {
  std::string key(128, '0');
  for (int i = 0; i < 512; i += 4) {
    int v = table[i] | (table[i + 1] << 1) | (table[i + 2] << 2) | (table[i + 3] << 3);
    key[i / 4] = "0123456789abcdef"[v];
  }
  return key;
}

// Rule string for RLE headers. Hensel/Life-like input is written as given. A MAP
// rule is written as a readable Hensel string when it is isotropic (so it is
// exactly representable), otherwise as the MAP string itself, which Golly reads.
static std::string RuleForHeader() {
  bool isMap = currentRuleName.size() > 3 && (currentRuleName.compare(0, 3, "MAP") == 0 ||
                                              currentRuleName.compare(0, 3, "map") == 0);
  if (isMap && RuleIsIsotropic())
    return CurrentRuleCanonical();
  return currentRuleName;
}

static std::string RuleDescription() {
  if (!RuleIsIsotropic())
    return currentRuleName.substr(0, 40) + (currentRuleName.size() > 40 ? "..." : "") +
           (compiledStep ? " (non-isotropic, compiled kernel" + compiledInfo + ")"
                         : " (non-isotropic, table-driven stepper" + compiledInfo + ")");
  return currentRuleName + " [" + CurrentRuleCanonical() + "]" +
         (useFastLife ? " (Conway: fast stepper)"
                      : compiledStep ? " (compiled kernel" + compiledInfo + ")"
                                     : " (table-driven stepper" + compiledInfo + ")");
}

// Table-driven step (any B1-B8 / S rule). Only columns next to a live cell
// are visited, and within a column only the 4-row groups next to a live cell;
// each group is one lookup in the window table (winTable). (A per-column memo
// cache used to sit in front of this; measured on recorded search states it was
// 1.5-2x slower than looking the windows up directly, so it was removed.)

static inline row_t StepColumnSlow(row_t a, row_t b, row_t c) {
  row_t m = a | b | c;
  row_t d = m | RotateLeft(m) | RotateRight(m);          // rows that can change
  row_t g = d | (d >> 1) | (d >> 2) | (d >> 3);          // bit 4k set if group k (rows 4k..4k+3) may change
  g &= ROW_EVERY4;
  row_t out = 0;
  while (g) {
    int pos = RowCtz(g);                                     // = 4k
    g &= g - 1;
    int sh = pos - 1;                                        // window starts one row above the group
    unsigned idx = (unsigned)(RotateRight(a, sh) & 0x3F) |
                   ((unsigned)(RotateRight(b, sh) & 0x3F) << 6) |
                   ((unsigned)(RotateRight(c, sh) & 0x3F) << 12);
    out |= (row_t)winTable[idx] << pos;
  }
  return out;
}


static void StepTable(row_t *state) {
  // Bit i of nz is set if column i has any live cell; only those columns and
  // their two neighbours can be non-empty in the next generation (every other
  // column is empty and stays empty, since B0 rules are rejected).
  row_t nz = 0;
  for (int i = 0; i < N; i++)
    nz |= (row_t)(state[i] != 0) << i;
  row_t todo = nz | RotateLeft(nz) | RotateRight(nz);

  // Columns are updated in place, in increasing order. Column i needs the old
  // values of columns i-1, i and i+1: i+1 is not updated yet, the old i-1 is
  // remembered, and the two wrap-around neighbours are saved up front.
  const row_t oldFirst = state[0], oldLast = state[N - 1];
  int lastIdx = -2;
  row_t lastOld = 0;
  while (todo) {
    int i = RowCtz(todo);
    todo &= todo - 1;
    row_t left = (i == 0) ? oldLast : (lastIdx == i - 1 ? lastOld : state[i - 1]);
    row_t right = (i == N - 1) ? oldFirst : state[i + 1];
    row_t cur = state[i];
    state[i] = StepColumnSlow(left, cur, right);
    lastIdx = i;
    lastOld = cur;
  }
}

// Bit i of the result is set iff column i has any live cell.
static inline row_t NonEmptyColumns(const row_t *state) {
  row_t nz = 0;
#if defined(__AVX2__) && N == 64
  for (int i = 0; i < N; i += 4) {
    __m256i v = _mm256_loadu_si256((const __m256i *)(state + i));
    unsigned zero = (unsigned)_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(v, _mm256_setzero_si256())));
    nz |= (uint64_t)(~zero & 0xF) << i;
  }
#else
  for (int i = 0; i < N; i++)
    nz |= (row_t)(state[i] != 0) << i;
#endif
  return nz;
}

#if defined(__GNUC__) || defined(__clang__)
// Four columns per SIMD vector (AVX2 on x86; two NEON registers on ARM).
// Kernels shift by RK_TOP (= row width - 1) to rotate a column by one row.
#define RK_TOP (N - 1)
#if N == 64
typedef uint64_t RuleVec __attribute__((vector_size(32)));

// One generation with a compiled column kernel K: the columns next to live
// cells are gathered, stepped four at a time, and written back.
template <RuleVec (*K)(RuleVec, RuleVec, RuleVec)>
static void StepCompiledTmpl(uint64_t *state) {
  uint64_t nz = NonEmptyColumns(state);
  uint64_t todo = nz | RotateLeft(nz) | RotateRight(nz);
  alignas(32) uint64_t L[N + 4], C[N + 4], R[N + 4], O[N + 4];
  int idx[N];
  int n = 0;
  while (todo) {
    int i = __builtin_ctzll(todo);
    todo &= todo - 1;
    idx[n] = i;
    L[n] = state[(i - 1) & (N - 1)];
    C[n] = state[i];
    R[n] = state[(i + 1) & (N - 1)];
    n++;
  }
  for (int j = n; j < ((n + 3) & ~3); j++)
    L[j] = C[j] = R[j] = 0;
  for (int j = 0; j < n; j += 4) {
    RuleVec l, c, r;
    memcpy(&l, L + j, sizeof l);
    memcpy(&c, C + j, sizeof c);
    memcpy(&r, R + j, sizeof r);
    RuleVec o = K(l, c, r);
    memcpy(O + j, &o, sizeof o);
  }
  for (int j = 0; j < n; j++)
    state[idx[j]] = O[j];
}

// Next generation of the selected columns only (see LifeState::StepColumnsInto).
template <RuleVec (*K)(RuleVec, RuleVec, RuleVec)>
static void StepColsCompiledTmpl(const uint64_t *state, uint64_t cols, uint64_t *out) {
  alignas(32) uint64_t L[N + 4], C[N + 4], R[N + 4], O[N + 4];
  int idx[N];
  int n = 0;
  while (cols) {
    int i = __builtin_ctzll(cols);
    cols &= cols - 1;
    idx[n] = i;
    L[n] = state[(i - 1) & (N - 1)];
    C[n] = state[i];
    R[n] = state[(i + 1) & (N - 1)];
    n++;
  }
  for (int j = n; j < ((n + 3) & ~3); j++)
    L[j] = C[j] = R[j] = 0;
  for (int j = 0; j < n; j += 4) {
    RuleVec l, c, r;
    memcpy(&l, L + j, sizeof l);
    memcpy(&c, C + j, sizeof c);
    memcpy(&r, R + j, sizeof r);
    RuleVec o = K(l, c, r);
    memcpy(O + j, &o, sizeof o);
  }
  for (int j = 0; j < n; j++)
    out[idx[j]] = O[j];
}
#else
// Wider torus: a column does not fit a SIMD lane, so the kernel runs on one
// whole column at a time.
typedef row_t RuleVec;

template <RuleVec (*K)(RuleVec, RuleVec, RuleVec)>
static void StepCompiledTmpl(row_t *state) {
  row_t nz = NonEmptyColumns(state);
  row_t todo = nz | RotateLeft(nz) | RotateRight(nz);
  row_t O[N];
  for (row_t bits = todo; bits; bits &= bits - 1) {
    int i = RowCtz(bits);
    O[i] = K(state[(i - 1) & (N - 1)], state[i], state[(i + 1) & (N - 1)]);
  }
  for (row_t bits = todo; bits; bits &= bits - 1) {
    int i = RowCtz(bits);
    state[i] = O[i];
  }
}

template <RuleVec (*K)(RuleVec, RuleVec, RuleVec)>
static void StepColsCompiledTmpl(const row_t *state, row_t cols, row_t *out) {
  for (; cols; cols &= cols - 1) {
    int i = RowCtz(cols);
    out[i] = K(state[(i - 1) & (N - 1)], state[i], state[(i + 1) & (N - 1)]);
  }
}
#endif

struct CompiledRuleEntry {
  const char *key;  // RuleTableKey() of the rule
  const char *rule; // the rule as written when it was compiled
  CompiledStepFn step;
  CompiledColsFn cols;
};

#if defined(__has_include)
#if __has_include("rule_kernels.h")
#include "rule_kernels.h"
#define CATFORCE_HAVE_COMPILED_RULES 1
#endif
#endif
#endif // __GNUC__ || __clang__

// Picks the compiled kernel for the current rule, if this binary has one. It is
// first checked against the table stepper on random patterns (so a stale or
// damaged rule_kernels.h can never change results), then timed against it on
// this machine; the faster one is used.
static void SelectCompiledKernel() {
  compiledStep = nullptr;
  compiledCols = nullptr;
  compiledInfo.clear();
#ifdef CATFORCE_HAVE_COMPILED_RULES
  if (useFastLife)
    return;
  std::string key = RuleTableKey();
  const CompiledRuleEntry *entry = nullptr;
  for (const CompiledRuleEntry &e : kCompiledRules)
    if (key == e.key)
      entry = &e;
  if (!entry) {
    compiledInfo = "; no compiled kernel for this rule, see --compile-rule";
    return;
  }

  // self-check
  uint64_t seed = 0x9E3779B97F4A7C15ULL;
  auto rnd = [&seed]() {
    seed ^= seed << 13;
    seed ^= seed >> 7;
    seed ^= seed << 17;
    return seed;
  };
  for (int t = 0; t < 200; t++) {
    row_t a[N], b[N];
    for (int i = 0; i < N; i++) {
      row_t w = 0;
      for (int k = 0; k < N; k += 64) {
        uint64_t part = rnd() & rnd();
        if (t % 2) part &= rnd();
        w |= (row_t)part << k;
      }
      if (t % 3 == 0 && i > 8 && i < N - 8) w = 0;
      a[i] = b[i] = w;
    }
    for (int s = 0; s < 4; s++) {
      StepTable(a);
      entry->step(b);
      if (memcmp(a, b, sizeof a) != 0) {
        std::cout << "Warning: the compiled kernel for this rule gives different results from the "
                     "table stepper; it is NOT used. Rebuild with --compile-rule." << std::endl;
        compiledInfo = "; compiled kernel failed its self-check";
        return;
      }
    }
  }

  // speed: same sparse soups through both
  auto timeOne = [&](CompiledStepFn f) {
    uint64_t s2 = 12345;
    auto start = std::chrono::steady_clock::now();
    uint64_t sink = 0;
    for (int t = 0; t < 400; t++) {
      row_t a[N] = {0};
      for (int i = 24; i < 40; i++) {
        s2 = s2 * 6364136223846793005ULL + 1442695040888963407ULL;
        uint64_t w = s2 ^ (s2 >> 29);
        a[i] = (row_t)(w & (w >> 7) & 0x00FFFF0000ULL);
      }
      for (int s = 0; s < 8; s++)
        f(a);
      sink += (uint64_t)a[30];
    }
    double ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count();
    return ns + (sink == 42 ? 1 : 0);
  };
  double tTable = timeOne(StepTable), tCompiled = timeOne(entry->step);
  tTable = std::min(tTable, timeOne(StepTable));
  tCompiled = std::min(tCompiled, timeOne(entry->step));
  if (tCompiled < tTable) {
    compiledStep = entry->step;
    compiledCols = entry->cols;
    char buf[64];
    snprintf(buf, sizeof buf, ", %.1fx the table stepper here", tTable / tCompiled);
    compiledInfo = buf;
  } else {
    compiledInfo = "; the compiled kernel is slower than it on this machine";
  }
#else
  if (!useFastLife)
    compiledInfo = "; --compile-rule builds a faster kernel for this rule";
#endif
}

void LifeState::Step() {
  if (!useFastLife) {
    gen++;
    if (compiledStep)
      compiledStep(state);
    else
      StepTable(state);
    return;
  }

  // Only columns that contain live cells, and their two neighbours, can be
  // non-empty next generation; every other column is empty and stays empty.
  row_t nz = NonEmptyColumns(state);
  row_t todo = nz | RotateLeft(nz) | RotateRight(nz);

  row_t newst[N]; // only the entries for columns in `todo` are written/read
  for (row_t bits = todo; bits; bits &= bits - 1) {
    int i = RowCtz(bits);
    row_t up = state[(i - 1) & (N - 1)], dn = state[(i + 1) & (N - 1)];
    // vertical 3-cell sums of the neighbouring columns (zero if they are empty)
    row_t ul = RotateLeft(up), ur = RotateRight(up);
    row_t dl = RotateLeft(dn), dr = RotateRight(dn);
    row_t uxor = ul ^ ur ^ up, uand = ((ul ^ ur) & up) | (ul & ur);
    row_t dxor = dl ^ dr ^ dn, dand = ((dl ^ dr) & dn) | (dl & dr);
    newst[i] = Rokicki(state[i], uxor, uand, dxor, dand);
  }
  for (row_t bits = todo; bits; bits &= bits - 1) {
    int i = RowCtz(bits);
    state[i] = newst[i];
  }

  min = 0;
  max = N - 1;
  gen++;
}

void LifeState::StepColumnsInto(row_t cols, row_t *out) const {
  if (!useFastLife && compiledCols) {
    compiledCols(state, cols, out);
    return;
  }
  if (!useFastLife) {
    for (row_t bits = cols; bits; bits &= bits - 1) {
      int i = RowCtz(bits);
      out[i] = StepColumnSlow(state[(i - 1) & (N - 1)], state[i], state[(i + 1) & (N - 1)]);
    }
    return;
  }
  for (row_t bits = cols; bits; bits &= bits - 1) {
    int i = RowCtz(bits);
    row_t up = state[(i - 1) & (N - 1)], dn = state[(i + 1) & (N - 1)];
    row_t ul = RotateLeft(up), ur = RotateRight(up);
    row_t dl = RotateLeft(dn), dr = RotateRight(dn);
    row_t uxor = ul ^ ur ^ up, uand = ((ul ^ ur) & up) | (ul & ur);
    row_t dxor = dl ^ dr ^ dn, dand = ((dl ^ dr) & dn) | (dl & dr);
    out[i] = Rokicki(state[i], uxor, uand, dxor, dand);
  }
}

void LifeState::Transform(SymmetryTransform transf) {
  switch (transf) {
  case Identity:
    break;
  case ReflectAcrossXEven:
    FlipX();
    break;
  case ReflectAcrossX:
    FlipX();
    Move(0, 1);
    break;
  case ReflectAcrossYEven:
    FlipY();
    break;
  case ReflectAcrossY:
    FlipY();
    Move(1, 0);
    break;
  case Rotate180EvenBoth:
    FlipX();
    FlipY();
    break;
  case Rotate180EvenVertical:
    FlipX();
    FlipY();
    Move(1, 0);
    break;
  case Rotate180EvenHorizontal:
    FlipX();
    FlipY();
    Move(0, 1);
    break;
  case Rotate180OddBoth:
    FlipX();
    FlipY();
    Move(1, 1);
    break;
  case ReflectAcrossYeqX:
    Transpose(false);
    break;
  case ReflectAcrossYeqNegX:
    Transpose(true);
    break;
  case ReflectAcrossYeqNegXP1:
    Transpose(true);
    Move(1, 1);
    break;
  case Rotate90Even:
    FlipX();
    Transpose(false);
    break;
  case Rotate90:
    FlipX();
    Transpose(false);
    Move(1, 0);
    break;
  case Rotate270Even:
    FlipY();
    Transpose(false);
    break;
  case Rotate270:
    FlipY();
    Transpose(false);
    Move(0, 1);
    break;
  }
  min = 0;
  max = N - 1;
}

void LifeState::Print() const {
  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      if (GetCell(i - N / 2, j - N / 2) == 0) {
        int hor = 0;
        int ver = 0;

        if ((j - N / 2) % 10 == 0)
          hor = 1;

        if ((i - N / 2) % 10 == 0)
          ver = 1;

        if (hor == 1 && ver == 1)
          printf("+");
        else if (hor == 1)
          printf("-");
        else if (ver == 1)
          printf("|");
        else
          printf(".");
      } else
        printf("O");
    }
    printf("\n");
  }

  printf("\n\n\n\n\n\n");
}

int LifeState::Parse(LifeState &state, const char *rle, int starti) {
  char ch;
  int cnt, i, j;
  int x, y;
  x = 0;
  y = 0;
  cnt = 0;

  i = starti;

  while ((ch = rle[i]) != '\0') {

    if (ch >= '0' && ch <= '9') {
      cnt *= 10;
      cnt += (ch - '0');
    } else if (ch == 'o') {

      if (cnt == 0)
        cnt = 1;

      for (j = 0; j < cnt; j++) {
        state.SetCell(x, y, 1);
        x++;
      }

      cnt = 0;
    } else if (ch == 'b') {
      if (cnt == 0)
        cnt = 1;

      x += cnt;
      cnt = 0;

    } else if (ch == '$') {
      if (cnt == 0)
        cnt = 1;

      if (cnt == 129)
        return i + 1;

      y += cnt;
      x = 0;
      cnt = 0;
    } else if (ch == '!') {
      break;
    } else {
      return -2;
    }

    i++;
  }

  state.RecalculateMinMax();

  return -1;
}

std::string LifeState::RLE() const {
  std::stringstream result;

  unsigned eol_count = 0;

  for (unsigned j = 0; j < N; j++) {
    bool last_val = GetCell(0 - N / 2, j - N / 2) == 1;
    unsigned run_count = 0;

    for (unsigned i = 0; i < N; i++) {
      bool val = GetCell(i - N / 2, j - N / 2) == 1;

      // Flush linefeeds if we find a live cell
      if (val && eol_count > 0) {
        if (eol_count > 1)
          result << eol_count;

        result << "$";

        eol_count = 0;
      }

      // Flush current run if val changes
      if (val == !last_val) {
        if (run_count > 1)
          result << run_count;

        if (last_val == 1)
          result << "o";
        else
          result << "b";

        run_count = 0;
      }

      run_count++;
      last_val = val;
    }

    // Flush run of live cells at end of line
    if (last_val) {
      if (run_count > 1)
        result << run_count;

      result << "o";

      run_count = 0;
    }

    eol_count++;
  }

  // Flush trailing linefeeds
  if (eol_count > 0) {
    if (eol_count > 1)
      result << eol_count;

    result << "$";

    eol_count = 0;
  }

  return result.str();
}

class LifeTarget {
public:
  LifeState wanted;
  LifeState unwanted;

  LifeTarget() {}
  LifeTarget(LifeState &state) {
    wanted = state;
    unwanted = state.GetBoundary();
  }

  void Transform(SymmetryTransform transf) {
    wanted.Transform(transf);
    unwanted.Transform(transf);
  }

  static int Parse(LifeTarget &target, const char *rle, int x, int y,
                   SymmetryTransform transf) {
    LifeState Temp;
    int result = LifeState::Parse(Temp, rle, x, y, transf);

    if (result == SUCCESS) {
      target.wanted = Temp;
      target.unwanted = Temp.GetBoundary();
      return SUCCESS;
    }

    return FAIL;
  }

  static LifeTarget Parse(const char *rle, int x, int y,
                          SymmetryTransform transf) {
    LifeTarget target;
    Parse(target, rle, x, y, transf);
    return target;
  }

  static LifeTarget Parse(const char *rle, int x, int y) {
    return Parse(rle, x, y, Identity);
  }

  static LifeTarget Parse(const char *rle) { return Parse(rle, 0, 0); }
};

inline bool LifeState::Contains(const LifeTarget &target, int dx,
                                int dy) const {
  return Contains(target.wanted, dx, dy) &&
         AreDisjoint(target.unwanted, dx, dy);
}

inline bool LifeState::Contains(const LifeTarget &target) const {
  return Contains(target.wanted) && AreDisjoint(target.unwanted);
}

inline LifeState LifeState::Match(const LifeTarget &target) const {
  return MatchLiveAndDead(target.wanted, target.unwanted);
}

// typedef struct {
//   int *xList;
//   int *yList;
//   int len;
//   int allocated;
// } Locator;

// Locator *NewLocator() {
//   Locator *result = (Locator *)(malloc(sizeof(Locator)));

//   result->xList = (int *)(malloc(sizeof(int)));
//   result->yList = (int *)(malloc(sizeof(int)));
//   result->len = 0;
//   result->allocated = 1;

//   return result;
// }

// Locator *Realloc(Locator *locator) {
//   if (locator->allocated <= locator->len) {
//     locator->allocated *= 2;
//     locator->xList =
//         (int *)(realloc(locator->xList, locator->allocated * sizeof(int)));
//     locator->yList =
//         (int *)(realloc(locator->yList, locator->allocated * sizeof(int)));
//   }
//   return locator;
// }

// void Add(Locator *locator, int x, int y) {
//   Realloc(locator);

//   locator->xList[locator->len] = x;
//   locator->yList[locator->len] = y;
//   locator->len++;
// }

// Locator *State2Locator(LifeState *state) {
//   Locator *result = NewLocator();

//   for (int j = 0; j < N; j++) {
//     for (int i = 0; i < N; i++) {
//       int val = state->Get(i, j);

//       if (val == 1)
//         Add(result, i, j);
//     }
//   }

//   return result;
// }

// void ClearAtX(LifeState *state, Locator *locator, int x, uint64_t val) {
//   if (val == 0ULL)
//     return;

//   int len = locator->len;
//   int *xList = locator->xList;
//   int *yList = locator->yList;

//   for (int i = 0; i < len; i++) {
//     int idx = (x + xList[i] + N) % N;
//     int circulate = (yList[i] + 64) % 64;

//     state->state[idx] &= ~RotateLeft(val, circulate);
//   }
// }

// uint64_t LocateAtX(LifeState *state, Locator *locator, int x, int negate) {
//   uint64_t result = ~0ULL;
//   int len = locator->len;
//   int *xList = locator->xList;
//   int *yList = locator->yList;

//   for (int i = 0; i < len; i++) {
//     int idx = (x + xList[i] + N) % N;
//     int circulate = (yList[i] + 64) % 64;

//     if (negate == false)
//       result &= RotateRight(state->state[idx], circulate);
//     else
//       result &= ~RotateRight(state->state[idx], circulate);

//     if (result == 0ULL)
//       break;
//   }

//   return result;
// }

// uint64_t LocateAtX(LifeState *state, Locator *onLocator, Locator *offLocator,
//                    int x) {
//   uint64_t onLocate = LocateAtX(state, onLocator, x, false);

//   if (onLocate == 0)
//     return 0;

//   return onLocate & LocateAtX(state, offLocator, x, true);
// }

// void LocateInRange(LifeState *state, Locator *locator, LifeState *result,
//                    int minx, int maxx, int negate) {
//   for (int i = minx; i <= maxx; i++) {
//     result->state[i] = LocateAtX(state, locator, i, negate);
//   }
// }

// void LocateInRange(LifeState *state, Locator *onLocator, Locator *offLocator,
//                    LifeState *result, int minx, int maxx) {
//   for (int i = minx; i <= maxx; i++) {
//     result->state[i] = LocateAtX(state, onLocator, offLocator, i);
//   }
// }

// void Locate(LifeState *state, Locator *locator, LifeState *result) {
//   LocateInRange(state, locator, result, state->min, state->max, false);
// }

// typedef struct {
//   Locator *onLocator;
//   Locator *offLocator;
// } TargetLocator;

// TargetLocator *NewTargetLocator() {
//   TargetLocator *result = (TargetLocator *)(malloc(sizeof(TargetLocator)));

//   result->onLocator = NewLocator();
//   result->offLocator = NewLocator();

//   return result;
// }

// TargetLocator *Target2Locator(LifeTarget *target) {
//   TargetLocator *result = (TargetLocator *)(malloc(sizeof(TargetLocator)));
//   result->onLocator = State2Locator(target->wanted);
//   result->offLocator = State2Locator(target->unwanted);

//   return result;
// }

// uint64_t LocateAtX(LifeState *state, TargetLocator *targetLocator, int x) {
//   return LocateAtX(state, targetLocator->onLocator,
//   targetLocator->offLocator,
//                    x);
// }

// void LocateInRange(LifeState *state, TargetLocator *targetLocator,
//                    LifeState *result, int minx, int maxx) {
//   return LocateInRange(state, targetLocator->onLocator,
//                        targetLocator->offLocator, result, minx, maxx);
// }

// void LocateTarget(LifeState *state, TargetLocator *targetLocator,
//                   LifeState *result) {
//   LocateInRange(state, targetLocator, result, state->min, state->max);
// }

// static TargetLocator *_glidersTarget[4];

// int RemoveAtX(LifeState *state, int x, int startGiderIdx) {
//   int removed = false;

//   for (int i = startGiderIdx; i < startGiderIdx + 2; i++) {
//     uint64_t gld = LocateAtX(state, _glidersTarget[i], x);

//     if (gld != 0) {
//       removed = true;
//       ClearAtX(state, _glidersTarget[i]->onLocator, x, gld);

//       for (int j = 0; j < 64; j++) {
//         if (gld % 2 == 1) {
//           // Append(state->emittedGliders, "(");
//           // Append(state->emittedGliders, i);
//           // Append(state->emittedGliders, ",");
//           // Append(state->emittedGliders, j);
//           // Append(state->emittedGliders, ",");
//           // Append(state->emittedGliders, state->gen);
//           // Append(state->emittedGliders, ",");
//           // Append(state->emittedGliders, x);
//           // Append(state->emittedGliders, ")");
//         }

//         gld = gld >> 1;

//         if (gld == 0)
//           break;
//       }
//     }
//   }

//   return removed;
// }

// void RemoveGliders(LifeState *state) {
//   int removed = false;

//   if (state->min <= 1)
//     if (RemoveAtX(state, 1, 0) == true)
//       removed = true;

//   if (state->max >= N - 2)
//     if (RemoveAtX(state, N - 2, 2) == true)
//       removed = true;

//   if (removed == true)
//     state->RecalculateMinMax();
// }

// void New() {
//   if (GlobalState == NULL) {
//     GlobalState = NewState();

//     // for (int i = 0; i < CAPTURE_COUNT; i++) {
//     //   Captures[i] = NewState();
//     // }

//     _glidersTarget[0] = NewTargetLocator("2o$obo$o!");
//     _glidersTarget[1] = NewTargetLocator("o$obo$2o!");
//     _glidersTarget[2] = NewTargetLocator("b2o$obo$2bo!", -2, 0);
//     _glidersTarget[3] = NewTargetLocator("2bo$obo$b2o!", -2, 0);

//   } else {
//     Clear(GlobalState);
//   }
// }
