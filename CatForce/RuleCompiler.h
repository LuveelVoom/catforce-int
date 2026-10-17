// RuleCompiler.h - writes a rule-specific stepping kernel as C++ source.
//
// The current rule (the global `table` set by SetRule) is turned into a
// straight-line, bit-parallel boolean circuit that computes the next state of
// all 64 cells of a column at once:
//   1. neighbour count bits (b0, t0, t1, t2) with an adder network,
//   2. the "totalistic" part (one value per centre state and neighbour count,
//      the majority of that group) as a small decision network over the count
//      bits, using the best of all variable orders,
//   3. the exceptions to it (the arrangements of a given count whose value
//      differs from the majority, e.g. 3q and 4i in B3/S23-q4i), each count's
//      set minimised as a sum of products over the 8 neighbours and the centre
//      using every other neighbour count as a "don't care", then factored.
// The result is XORed together. The kernel is a template so the stepper can run
// it on SIMD vectors of 4 columns (GCC/Clang vector extensions; AVX2 on x86,
// NEON on ARM). Correctness does not depend on any of the heuristics: every
// choice only changes the size of the circuit, never its truth table, and the
// generated header is checked against the table stepper when it is loaded.

#pragma once
#include <algorithm>
#include <functional>
#include <map>
#include <tuple>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace rulecompiler {

class Program {
public:
  std::vector<std::string> lines;
  std::map<std::string, std::string> memo;
  int counter = 0;

  static bool IsConst(const std::string &v) { return v == "ZERO" || v == "ONES"; }

  std::string Op(const std::string &kind, std::string a, std::string b = "") {
    if (kind == "and") {
      if (a == "ZERO" || b == "ZERO") return "ZERO";
      if (a == "ONES") return b;
      if (b == "ONES") return a;
      if (a == b) return a;
      if (b < a) std::swap(a, b);
    } else if (kind == "or") {
      if (a == "ONES" || b == "ONES") return "ONES";
      if (a == "ZERO") return b;
      if (b == "ZERO") return a;
      if (a == b) return a;
      if (b < a) std::swap(a, b);
    } else if (kind == "xor") {
      if (a == "ZERO") return b;
      if (b == "ZERO") return a;
      if (a == b) return "ZERO";
      if (a == "ONES") return Op("not", b);
      if (b == "ONES") return Op("not", a);
      if (b < a) std::swap(a, b);
    } else if (kind == "andn") { // a & ~b
      if (a == "ZERO" || b == "ONES") return "ZERO";
      if (b == "ZERO") return a;
      if (a == "ONES") return Op("not", b);
      if (a == b) return "ZERO";
    } else if (kind == "not") {
      if (a == "ZERO") return "ONES";
      if (a == "ONES") return "ZERO";
    }
    std::string key = kind + "|" + a + "|" + b;
    auto it = memo.find(key);
    if (it != memo.end()) return it->second;
    std::string v = "v" + std::to_string(counter++);
    std::string expr;
    if (kind == "and") expr = a + " & " + b;
    else if (kind == "or") expr = a + " | " + b;
    else if (kind == "xor") expr = a + " ^ " + b;
    else if (kind == "andn") expr = a + " & ~" + b;
    else expr = "~" + a;
    lines.push_back("  const T " + v + " = " + expr + ";");
    memo[key] = v;
    return v;
  }

  // x ? hi : lo
  std::string Mux(const std::string &x, const std::string &hi, const std::string &lo) {
    if (hi == lo) return hi;
    if (hi == "ONES" && lo == "ZERO") return x;
    if (hi == "ZERO" && lo == "ONES") return Op("not", x);
    if (lo == "ZERO") return Op("and", x, hi);
    if (hi == "ZERO") return Op("andn", lo, x);
    if (hi == "ONES") return Op("or", x, lo);
    if (lo == "ONES") return Op("or", Op("not", x), hi);
    return Op("or", Op("and", x, hi), Op("andn", lo, x));
  }
};

// ---- totalistic base: decision network over C, b0, t0, t1, t2 -------------
// value[a] for a 5-bit assignment (bit0 C, bit1 b0, bit2 t0, bit3 t1, bit4 t2):
// 0, 1, or -1 for "don't care" (neighbour counts above 8 cannot happen).
struct Node {
  int var;      // -1 for a constant
  int constVal; // for constants
  int lo, hi;   // child node ids
};

class Bdd {
public:
  std::vector<Node> nodes;
  std::map<std::tuple<int, int, int>, int> unique;
  int zero, one;
  Bdd() {
    nodes.push_back({-1, 0, -1, -1});
    nodes.push_back({-1, 1, -1, -1});
    zero = 0;
    one = 1;
  }
  int Make(int var, int lo, int hi) {
    if (lo == hi) return lo;
    auto key = std::make_tuple(var, lo, hi);
    auto it = unique.find(key);
    if (it != unique.end()) return it->second;
    nodes.push_back({var, 0, lo, hi});
    unique[key] = (int)nodes.size() - 1;
    return (int)nodes.size() - 1;
  }
  // Builds the diagram for `value` with the given variable order; don't-cares
  // take the value of the sibling branch when one branch is entirely don't-care.
  int Build(const int *value, const int *order, int depth, int fixedMask, int fixedVal) {
    if (depth == 5) {
      for (int a = 0; a < 32; a++)
        if ((a & fixedMask) == fixedVal && value[a] >= 0)
          return value[a] ? one : zero;
      return -2; // all don't-care
    }
    int v = order[depth];
    int lo = Build(value, order, depth + 1, fixedMask | (1 << v), fixedVal);
    int hi = Build(value, order, depth + 1, fixedMask | (1 << v), fixedVal | (1 << v));
    if (lo == -2) return hi;
    if (hi == -2) return lo;
    return Make(v, lo, hi);
  }
  int CountReachable(int root) {
    std::set<int> seen;
    std::vector<int> stack{root};
    while (!stack.empty()) {
      int n = stack.back();
      stack.pop_back();
      if (n < 2 || seen.count(n)) continue;
      seen.insert(n);
      stack.push_back(nodes[n].lo);
      stack.push_back(nodes[n].hi);
    }
    return (int)seen.size();
  }
};

// ---- two-level minimisation with don't-cares over 9 variables -------------
// Points are 9-bit: bits 0..7 = neighbour mask bits, bit 8 = centre.
struct Cube {
  int care, val;
};

static std::vector<Cube> Minimize(const std::vector<int> &onset, const std::vector<int> &offset) {
  std::vector<Cube> result;
  if (onset.empty()) return result;
  std::vector<char> isOff(512, 0);
  for (int p : offset) isOff[p] = 1;

  // all valid cubes containing at least one onset point
  std::set<std::pair<int, int>> valid;
  std::vector<Cube> cubes;
  for (int care = 0; care < 512; care++) {
    std::set<int> vals;
    for (int p : onset) vals.insert(p & care);
    for (int v : vals) {
      bool ok = true;
      for (int q : offset)
        if ((q & care) == v) { ok = false; break; }
      if (ok) {
        cubes.push_back({care, v});
        valid.insert({care, v});
      }
    }
  }
  // primes: cubes that cannot drop a literal and stay valid
  std::vector<Cube> primes;
  for (const Cube &c : cubes) {
    bool dominated = false;
    for (int b = 0; b < 9 && !dominated; b++)
      if ((c.care >> b) & 1)
        if (valid.count({c.care & ~(1 << b), c.val & ~(1 << b)})) dominated = true;
    if (!dominated) primes.push_back(c);
  }
  // greedy cover, cost = literals + 1
  std::set<int> uncovered(onset.begin(), onset.end());
  std::vector<std::vector<int>> covers(primes.size());
  for (size_t i = 0; i < primes.size(); i++)
    for (int p : onset)
      if ((p & primes[i].care) == primes[i].val) covers[i].push_back(p);
  std::vector<size_t> chosen;
  while (!uncovered.empty()) {
    double bestScore = -1;
    size_t best = 0;
    for (size_t i = 0; i < primes.size(); i++) {
      int gain = 0;
      for (int p : covers[i]) gain += uncovered.count(p);
      double score = gain / (double)(__builtin_popcount(primes[i].care) + 1);
      if (score > bestScore) { bestScore = score; best = i; }
    }
    chosen.push_back(best);
    for (int p : covers[best]) uncovered.erase(p);
  }
  // drop redundant cubes
  for (bool changed = true; changed;) {
    changed = false;
    for (size_t k = 0; k < chosen.size(); k++) {
      std::set<int> rest;
      for (size_t j = 0; j < chosen.size(); j++)
        if (j != k) rest.insert(covers[chosen[j]].begin(), covers[chosen[j]].end());
      bool redundant = true;
      for (int p : onset)
        if (!rest.count(p)) { redundant = false; break; }
      if (redundant) { chosen.erase(chosen.begin() + k); changed = true; break; }
    }
  }
  for (size_t i : chosen) result.push_back(primes[i]);
  return result;
}

static const char *const kVarName[9] = {"n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "C"};

// Greedy factoring of a sum of cubes (each a map var -> value).
static std::string Factor(Program &P, const std::vector<std::map<int, int>> &cubes) {
  if (cubes.empty()) return "ZERO";
  for (const auto &c : cubes)
    if (c.empty()) return "ONES";
  std::map<std::pair<int, int>, int> count;
  for (const auto &c : cubes)
    for (const auto &kv : c) count[kv]++;
  std::pair<int, int> best = count.begin()->first;
  for (const auto &kv : count)
    if (kv.second > count[best]) best = kv.first;
  std::vector<std::map<int, int>> with, without;
  for (const auto &c : cubes) {
    auto it = c.find(best.first);
    if (it != c.end() && it->second == best.second) {
      std::map<int, int> rest = c;
      rest.erase(best.first);
      with.push_back(rest);
    } else {
      without.push_back(c);
    }
  }
  std::string inner = Factor(P, with);
  std::string part = best.second ? P.Op("and", kVarName[best.first], inner) : P.Op("andn", inner, kVarName[best.first]);
  return P.Op("or", part, Factor(P, without));
}

// birth[m] / survive[m] indexed by the 8-bit neighbour mask (7 NW ... 0 SE).
// Returns the C++ text of `template <class T> static inline T name(T L, T C, T R)`;
// `ops` receives the number of operations after the neighbour count.
static std::string GenerateKernel(const bool *birth, const bool *survive, const std::string &name,
                                  const std::string &ruleText, int &ops) {
  Program P;
  int base[2][9];
  std::vector<int> exc[2][9];
  for (int c = 0; c < 2; c++)
    for (int k = 0; k <= 8; k++) {
      std::vector<int> on, off;
      for (int m = 0; m < 256; m++)
        if (__builtin_popcount(m) == k) ((c ? survive : birth)[m] ? on : off).push_back(m);
      base[c][k] = on.size() > off.size() ? 1 : 0;
      exc[c][k] = base[c][k] ? off : on;
    }

  // base function, best variable order
  int value[32];
  for (int a = 0; a < 32; a++) {
    int cnt = ((a >> 1) & 1) + 2 * ((a >> 2) & 1) + 4 * ((a >> 3) & 1) + 8 * ((a >> 4) & 1);
    value[a] = cnt <= 8 ? base[a & 1][cnt] : -1;
  }
  int order[5] = {0, 1, 2, 3, 4}, bestOrder[5] = {0, 1, 2, 3, 4};
  int bestSize = 1 << 30;
  do {
    Bdd b;
    int root = b.Build(value, order, 0, 0, 0);
    if (root == -2) root = b.zero;
    int size = b.CountReachable(root);
    if (size < bestSize) {
      bestSize = size;
      std::copy(order, order + 5, bestOrder);
    }
  } while (std::next_permutation(order, order + 5));
  Bdd bdd;
  int root = bdd.Build(value, bestOrder, 0, 0, 0);
  if (root == -2) root = bdd.zero;
  const char *baseVar[5] = {"C", "b0", "t0", "t1", "t2"};
  std::map<int, std::string> emitted;
  std::function<std::string(int)> emit = [&](int n) -> std::string {
    if (n == bdd.zero) return "ZERO";
    if (n == bdd.one) return "ONES";
    auto it = emitted.find(n);
    if (it != emitted.end()) return it->second;
    std::string r = P.Mux(baseVar[bdd.nodes[n].var], emit(bdd.nodes[n].hi), emit(bdd.nodes[n].lo));
    emitted[n] = r;
    return r;
  };
  std::string baseExpr = emit(root);

  // count == k (count = b0 + 2 t0 + 4 t1 + 8 t2, at most 8)
  auto eqk = [&](int k) -> std::string {
    if (k == 8) return "t2";
    // k < 8: the low three count bits must match k and t2 must be 0. When t2 is
    // set the count is exactly 8, so b0 = t0 = t1 = 0; any k with one of those
    // bits set therefore already excludes it, and only k = 0 needs the check.
    std::string e = "ONES";
    const char *bits[3] = {"b0", "t0", "t1"};
    for (int j = 0; j < 3; j++)
      if ((k >> j) & 1) e = P.Op("and", e, bits[j]);
    for (int j = 0; j < 3; j++)
      if (!((k >> j) & 1)) e = P.Op("andn", e, bits[j]);
    if (k == 0) e = P.Op("andn", e, "t2");
    return e;
  };

  std::string X = "ZERO";
  for (int k = 0; k <= 8; k++) {
    std::vector<int> on, off;
    for (int c = 0; c < 2; c++) {
      std::set<int> ex(exc[c][k].begin(), exc[c][k].end());
      for (int m = 0; m < 256; m++) {
        if (__builtin_popcount(m) != k) continue;
        (ex.count(m) ? on : off).push_back(m | (c << 8));
      }
    }
    if (on.empty()) continue;
    std::vector<Cube> cubes = Minimize(on, off);
    std::vector<std::map<int, int>> lits;
    for (const Cube &cb : cubes) {
      std::map<int, int> c;
      for (int i = 0; i < 9; i++)
        if ((cb.care >> i) & 1) c[i] = (cb.val >> i) & 1;
      lits.push_back(c);
    }
    std::string M = Factor(P, lits);
    X = P.Op("or", X, P.Op("and", M, eqk(k)));
  }
  std::string res = P.Op("xor", baseExpr, X);
  if (res == "ZERO") res = "T{}";
  if (res == "ONES") res = "~T{}";
  ops = (int)P.lines.size();

  std::stringstream ss;
  ss << "// rule " << ruleText << ": " << ops << " operations after the neighbour count\n";
  ss << "template <class T> static inline T " << name << "(T L, T C, T R) {\n";
  ss << "  const T n7 = (L << 1) | (L >> 63), n6 = (C << 1) | (C >> 63), n5 = (R << 1) | (R >> 63), n4 = L, n3 = R;\n";
  ss << "  const T n2 = (L >> 1) | (L << 63), n1 = (C >> 1) | (C << 63), n0 = (R >> 1) | (R << 63);\n";
  ss << "  const T a1 = n0 ^ n1, c1 = n0 & n1, a2 = n2 ^ n3, c2 = n2 & n3, a3 = n4 ^ n5, c3 = n4 & n5, a4 = n6 ^ n7, c4 = n6 & n7;\n";
  ss << "  const T x0 = a1 ^ a2, y0 = a1 & a2, x1 = a3 ^ a4, y1 = a3 & a4;\n";
  ss << "  const T b0 = x0 ^ x1, z = x0 & x1;\n";
  ss << "  const T p = c1 ^ c2, s1 = p ^ c3, k1 = (c1 & c2) | (p & c3);\n";
  ss << "  const T q = c4 ^ y0, s2 = q ^ y1, k2 = (c4 & y0) | (q & y1);\n";
  ss << "  const T r = s1 ^ s2, t0 = r ^ z, k3 = (s1 & s2) | (r & z);\n";
  ss << "  const T u = k1 ^ k2, t1 = u ^ k3, t2 = (k1 & k2) | (u & k3);\n";
  ss << "  (void)n0; (void)n1; (void)n2; (void)n3; (void)n4; (void)n5; (void)n6; (void)n7;\n";
  ss << "  (void)b0; (void)t0; (void)t1; (void)t2;\n";
  for (const std::string &l : P.lines) ss << l << "\n";
  ss << "  return " << res << ";\n}\n";
  return ss.str();
}

} // namespace rulecompiler
