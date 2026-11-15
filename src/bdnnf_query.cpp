/******************************************
Copyright (C) 2023 Authors of GANAK, see AUTHORS file

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

// bdnnf count / bdnnf sample: the exact (weighted, projected) model count of a
// BDNNF circuit in one streamed pass, one var bitset and one rational per node.
// Smooths on the fly, so it needs a decomposable, deterministic circuit but not
// a SMOOTH one. Weights ('c p weight', Arjun's semantics) come from --cnf; the
// projection is the file's "projection" header key (a --cnf must agree).
// --assume L1,L2,... conditions on literals of projection vars.
//
// sample draws N models, each with probability proportional to its weight
// (uniformly when unweighted) over the projection, top-down through the frame
// index; non-projection vars come from the path, the definitions, or at
// random. Models are printed as DIMACS lines over the projection (--full: all
// vars), to stdout or -o. With --cnf, each is checked against the clauses.

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>

#include <gmpxx.h>

#include "bdnnf.hpp"

namespace {

// gmpxx has no unsigned long long operators, and uint64_t is one on wasm32
mpz_class mpz_u64(uint64_t x) {
  mpz_class z = (unsigned long)(x >> 32);
  z <<= 32;
  z += (unsigned long)(x & 0xffffffffULL);
  return z;
}
uint64_t u64_mpz(const mpz_class& z) {
  const mpz_class hi = z >> 32;
  const mpz_class lo = z - (hi << 32);
  return ((uint64_t)hi.get_ui() << 32) | (uint64_t)lo.get_ui();
}

// Whether a literal list holds a literal and its negation, in linear time
// (arcs can carry hundreds of implied literals).
bool has_conflict(const std::vector<int>& ls) {
  if (ls.size() < 16) {
    for (size_t a = 0; a < ls.size(); a++)
      for (size_t b = a + 1; b < ls.size(); b++) if (ls[a] == -ls[b]) return true;
    return false;
  }
  static thread_local std::vector<uint32_t> stamp;
  static thread_local uint32_t cur = 0;
  if (++cur == 0) {
    std::fill(stamp.begin(), stamp.end(), 0);
    cur = 1;
  }
  for (int l : ls) {
    const size_t i = 2 * (size_t)std::abs(l) + (l < 0 ? 1 : 0);
    if (i + 2 > stamp.size()) stamp.resize(2 * (i + 2), 0);
    if (stamp[i ^ 1] == cur) return true;
    stamp[i] = cur;
  }
  return false;
}

struct QueryError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// throws, so that bdnnf stream survives a bad query
[[noreturn]] void die(const std::string& msg) { throw QueryError(msg); }

// cheap internal checks, on in every build
#define BDNNF_CHECK(c) \
  do { \
    if (!(c)) die(std::string("internal error: ") + #c + " (" + __FILE__ + ":" + std::to_string(__LINE__) + ")"); \
  } while (0)

// "a", "a/b", "-0.25", "1e-3"-free decimals.
mpq_class parse_num(const std::string& s) {
  mpq_class r;
  const auto dot = s.find('.');
  if (dot == std::string::npos) {
    if (r.set_str(s, 10) != 0) die("bad weight: " + s);
  } else {
    std::string digits = s.substr(0, dot) + s.substr(dot + 1);
    if (digits.empty() || digits == "-" || digits == "+") digits += "0";
    mpz_class num(digits, 10), den = 1;
    for (size_t i = dot + 1; i < s.size(); i++) den *= 10;
    r = mpq_class(num, den);
  }
  r.canonicalize();
  return r;
}

struct Weights {
  bool weighted = false;
  std::map<int, mpq_class> w;
  bool projected = false;
  std::vector<char> in_proj;  // by var
};

Weights read_cnf(const std::string& path, uint32_t nvars) {
  std::ifstream f(path);
  if (!f.good()) die("cannot open " + path);
  Weights W;
  std::vector<char> proj(nvars + 1, 0);
  std::string line;
  while (std::getline(f, line)) {
    std::istringstream ss(line);
    std::string a, b, c;
    ss >> a >> b;
    if (a != "c") continue;
    if (b == "ind") {
      W.projected = true;
      long v;
      while (ss >> v && v != 0) if (v > 0 && (uint64_t)v <= nvars) proj[v] = 1;
      continue;
    }
    if (b != "p") continue;
    ss >> c;
    if (c == "weight") {
      long l;
      std::string ws;
      ss >> l >> ws;
      const mpq_class x = parse_num(ws);
      W.weighted = true;
      if (!W.w.count((int)l) && !W.w.count((int)-l)) W.w[(int)-l] = 1 - x;
      W.w[(int)l] = x;
    } else if (c == "show") {
      W.projected = true;
      long v;
      while (ss >> v && v != 0) if (v > 0 && (uint64_t)v <= nvars) proj[v] = 1;
    }
  }
  if (!W.projected) proj.assign(nvars + 1, 1);
  W.in_proj = std::move(proj);
  return W;
}

std::vector<std::vector<int>> read_clauses(const std::string& path) {
  std::ifstream f(path);
  std::vector<std::vector<int>> cls;
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == 'c' || line[0] == 'p') continue;
    std::istringstream ss(line);
    std::vector<int> cl;
    long x;
    while (ss >> x && x != 0) cl.push_back((int)x);
    if (!cl.empty()) cls.push_back(std::move(cl));
  }
  return cls;
}

// Random access to nodes through the frame index, one decoded frame cached.
struct NodeAccess {
  bdnnf::Reader& r;
  size_t cur = SIZE_MAX;
  std::vector<uint8_t> raw;
  std::vector<size_t> starts;
  explicit NodeAccess(bdnnf::Reader& rd) : r(rd) {}
  void get(uint64_t id, bdnnf::Node& n) {
    if (r.is_cached()) {
      if (id < 2 || id - 2 >= r.cache.size()) throw bdnnf::Error("node id out of range");
      n = r.cache[id - 2];
      return;
    }
    const size_t fi = r.frame_of(id);
    if (fi != cur) {
      r.load_frame(fi, raw);
      starts.clear();
      bdnnf::Cursor c(raw.data(), raw.size());
      for (uint32_t i = 0; i < r.frames[fi].nnodes; i++) {
        starts.push_back((size_t)(c.p - raw.data()));
        bdnnf::decode_node(c, r.frames[fi].first_id + i, r.nvars, n);
      }
      cur = fi;
    }
    const size_t k = (size_t)(id - r.frames[fi].first_id);
    bdnnf::Cursor c(raw.data() + starts[k], raw.size() - starts[k]);
    bdnnf::decode_node(c, id, r.nvars, n);
  }
};

// a * (product of fs except fs[i]) for every i, without dividing.
void excl_products(const std::vector<mpq_class>& fs, const mpq_class& a, std::vector<mpq_class>& out) {
  const size_t k = fs.size();
  out.assign(k, a);
  mpq_class acc = 1;
  for (size_t i = 0; i < k; i++) { out[i] *= acc; acc *= fs[i]; }
  acc = 1;
  for (size_t i = k; i-- > 0;) { out[i] *= acc; acc *= fs[i]; }
}

// Weighted count of the models containing each projection literal (index
// 2v for v, 2v+1 for -v): the derivatives of the count, top-down in one
// backward scan. A literal's count is w(l) times the derivative by w(l).
template <class WFn>
std::vector<mpq_class> lit_counts(bdnnf::Reader& r, const std::vector<mpq_class>& val,
                                  const std::vector<char>& dead, const std::vector<uint64_t>& vs,
                                  size_t WD, const std::vector<char>& in_proj, WFn&& w) {
  const uint32_t nv = r.nvars;
  std::vector<mpq_class> cnt(2 * (size_t)nv + 2), dv(val.size());
  auto row = [&](uint64_t id) { return vs.data() + id * WD; };
  auto has = [](const uint64_t* s, uint32_t v) { return (s[v >> 6] >> (v & 63)) & 1; };
  auto lidx = [](int l) { return 2 * (size_t)std::abs(l) + (l < 0 ? 1 : 0); };
  auto both = [&](uint32_t v) -> mpq_class { return w((int)v) + w(-(int)v); };
  std::vector<mpq_class> fs, ex;
  std::vector<uint32_t> gap;
  // Each var of `gap` is a factor (w(v) + w(-v)) at fs[off + i].
  auto add_gap = [&](size_t off) {
    for (size_t i = 0; i < gap.size(); i++) {
      cnt[lidx((int)gap[i])] += ex[off + i] * w((int)gap[i]);
      cnt[lidx(-(int)gap[i])] += ex[off + i] * w(-(int)gap[i]);
    }
  };
  if (dead[r.root]) return cnt;
  fs.assign(1, val[r.root]);
  gap.clear();
  for (uint32_t v = 1; v <= nv; v++)
    if (in_proj[v] && !has(row(r.root), v)) { gap.push_back(v); fs.push_back(both(v)); }
  excl_products(fs, 1, ex);
  dv[r.root] = ex[0];
  add_gap(1);
  std::vector<uint64_t> av(WD);
  r.scan_backward([&](uint64_t id, const bdnnf::Node& nd) {
    if (dead[id] || dv[id] == 0) return;
    const mpq_class& g = dv[id];
    if (nd.kind == bdnnf::K_AND) {
      fs.clear();
      for (uint64_t k : nd.kids) fs.push_back(val[k]);
      gap.clear();
      for (uint32_t v : nd.free) if (in_proj[v]) { gap.push_back(v); fs.push_back(both(v)); }
      mpq_class u = 1;
      for (int l : nd.units) u *= w(l);
      excl_products(fs, g * u, ex);
      for (size_t i = 0; i < nd.kids.size(); i++) dv[nd.kids[i]] += ex[i];
      add_gap(nd.kids.size());
      for (int l : nd.units) if (in_proj[(uint32_t)std::abs(l)]) cnt[lidx(l)] += g * val[id];
      return;
    }
    const uint64_t* me = row(id);
    for (const auto& a : nd.arcs) {
      if (dead[a.child]) continue;
      if (has_conflict(a.lits)) continue;
      const uint64_t* cv = row(a.child);
      av.assign(cv, cv + WD);
      mpq_class u = 1;
      for (int l : a.lits) {
        u *= w(l);
        const uint32_t v = (uint32_t)std::abs(l);
        if (in_proj[v]) av[v >> 6] |= 1ULL << (v & 63);
      }
      fs.assign(1, val[a.child]);
      gap.clear();
      for (size_t q = 0; q < WD; q++)
        for (uint64_t m = me[q] & ~av[q]; m; m &= m - 1) {
          gap.push_back((uint32_t)(q * 64 + __builtin_ctzll(m)));
          fs.push_back(both(gap.back()));
        }
      excl_products(fs, g * u, ex);
      dv[a.child] += ex[0];
      add_gap(1);
      const mpq_class y = ex[0] * val[a.child];
      for (int l : a.lits) if (in_proj[(uint32_t)std::abs(l)]) cnt[lidx(l)] += y;
    }
  });
  return cnt;
}

// k-best combination of independent factors: each list sorted descending; a
// stage keeps, for each of its best k sums, the indices it came from.
struct KItem {
  mpq_class v;
  uint32_t a, b;
};
using KStage = std::vector<KItem>;

void kcombine(const std::vector<mpq_class>& A, const std::vector<mpq_class>& B, size_t k, KStage& out) {
  out.clear();
  for (uint32_t i = 0; i < A.size(); i++)
    for (uint32_t j = 0; j < B.size(); j++) out.push_back({A[i] + B[j], i, j});
  const size_t m = std::min(k, out.size());
  std::partial_sort(out.begin(), out.begin() + (long)m, out.end(), [](const KItem& x, const KItem& y) {
    return x.v != y.v ? x.v > y.v : (x.a != y.a ? x.a < y.a : x.b < y.b);
  });
  out.resize(m);
}

struct KFactor {
  uint64_t child = 0;            // a node, or 0 for a var
  std::vector<int> lits;         // a var: the literal of each entry
  std::vector<mpq_class> vals;
};

// The stages of base + factors[0] + factors[1] + ...; the last one is the list.
void kseq(const std::vector<KFactor>& fs, const mpq_class& base, size_t k, std::vector<KStage>& st) {
  st.assign(fs.size() + 1, {});
  st[0] = {{base, 0, 0}};
  std::vector<mpq_class> cur{base};
  for (size_t i = 0; i < fs.size(); i++) {
    kcombine(cur, fs[i].vals, k, st[i + 1]);
    cur.clear();
    for (const auto& it : st[i + 1]) cur.push_back(it.v);
  }
}

std::map<int, mpq_class> read_objective(const std::string& path, uint32_t nv) {
  std::ifstream f(path);
  if (!f.good()) die("cannot open " + path);
  std::map<int, mpq_class> o;
  std::string line;
  while (std::getline(f, line)) {
    std::istringstream ss(line);
    std::string a, b;
    if (!(ss >> a) || a[0] == 'c') continue;
    if (!(ss >> b)) die("objective line without a value: " + line);
    const long l = std::stol(a);
    if (l == 0 || (uint64_t)std::abs(l) > nv) die("bad objective literal: " + a);
    o[(int)l] = parse_num(b);
  }
  return o;
}

// The t-subsets of {0..p-1} in colex order: rank increases by one per next().
struct Combo {
  std::vector<uint32_t> c;
  uint32_t p;
  Combo(uint32_t t, uint32_t _p) : p(_p) { for (uint32_t i = 0; i < t; i++) c.push_back(i); }
  bool next() {
    const size_t t = c.size();
    for (size_t i = 0; i < t; i++) {
      const uint32_t lim = i + 1 < t ? c[i + 1] : p;
      if (c[i] + 1 < lim) {
        c[i]++;
        for (size_t j = 0; j < i; j++) c[j] = (uint32_t)j;
        return true;
      }
    }
    return false;
  }
};

// t-wise sampling, greedy as in YASA: take the first t-tuple of projection
// literals no sample covers, add further uncovered tuples while the circuit
// says the set stays satisfiable (64 sets per boolean pass), then complete
// it to a model by a random walk over satisfiable arcs. A tuple unsatisfiable
// on its own is skipped. Every satisfiable tuple ends up covered.
// For t = 2 (if the pair matrix fits in 512 MB) the set instead grows by the
// literals that cover most uncovered pairs with it, among those a backward
// pass shows to be in a model of it, and the walk picks arcs and free vars by
// the uncovered pairs they make.
void twise(bdnnf::Reader& r, const std::vector<char>& dead, const std::vector<char>& in_proj,
           const std::vector<int>& assume, uint32_t t, uint64_t seed, bool full,
           const std::vector<std::vector<int>>& cls, const std::vector<char>& lit_dead,
           std::ostream& out) {
  const uint32_t nv = r.nvars;
  std::vector<char> defined(nv + 1, 0);
  for (const auto& [v, rf] : r.defs.defs) defined[v] = 1;
  std::vector<uint32_t> P;
  for (uint32_t v = 1; v <= nv; v++) if (in_proj[v] && !defined[v]) P.push_back(v);
  const uint32_t p = (uint32_t)P.size();
  t = std::min(t, p);
  uint64_t ncombo = 1;
  for (uint32_t i = 0; i < t; i++) {
    ncombo = ncombo * (p - i) / (i + 1);
    if (ncombo > (1ULL << 32)) die("too many t-wise interactions");
  }
  const uint64_t total = ncombo << t;
  // infeasible: unsatisfiable on its own; feasible: known satisfiable on its own
  std::vector<uint64_t> infeasible((total >> 6) + 1), feasible((total >> 6) + 1);
  // smask[2i + neg]: the samples containing the literal of P[i]; a tuple is
  // covered when its literals' masks intersect
  std::vector<std::vector<uint64_t>> smask(2 * (size_t)p);
  // t = 2: cov[a] has bit b when literals a and b (2i + neg for P[i]) share a
  // sample, and the config grows by the possible literals that cover most
  // uncovered pairs with it
  const size_t np2 = 2 * (size_t)p, cw = (np2 + 63) / 64;
  const bool greedy2 = t == 2 && np2 * cw <= (64ULL << 20);
  std::vector<uint64_t> cov(greedy2 ? np2 * cw : 0);
  auto covered = [&](const Combo& c, uint32_t pol) {
    if (greedy2) {
      const size_t a = 2 * (size_t)c.c[0] + (pol & 1), b = 2 * (size_t)c.c[1] + ((pol >> 1) & 1);
      return ((cov[a * cw + (b >> 6)] >> (b & 63)) & 1) != 0;
    }
    const size_t words = smask.empty() ? 0 : smask[0].size();
    for (size_t w = 0; w < words; w++) {
      uint64_t m = ~0ULL;
      for (uint32_t i = 0; i < t && m; i++) m &= smask[2 * (size_t)c.c[i] + ((pol >> i) & 1)][w];
      if (m) return true;
    }
    return false;
  };
  auto bit = [](const std::vector<uint64_t>& b, uint64_t i) { return (b[i >> 6] >> (i & 63)) & 1; };
  auto setbit = [](std::vector<uint64_t>& b, uint64_t i) { b[i >> 6] |= 1ULL << (i & 63); };
  auto lidx = [](int l) { return 2 * (size_t)std::abs(l) + (l < 0 ? 1 : 0); };
  std::vector<uint64_t> kill(2 * (size_t)nv + 2), sm(dead.size());
  auto killed = [&](const std::vector<int>& ls) {
    uint64_t m = 0;
    for (int l : ls) if (in_proj[(uint32_t)std::abs(l)]) m |= kill[lidx(l)];
    return m;
  };
  auto conflicting = has_conflict;
  // bit j of the result: set j (plus --assume) is satisfiable
  // with kill set up for the sets in the lanes of `all`
  auto scan_pass = [&](uint64_t all) {
    for (int l : assume) kill[lidx(-l)] |= all;
    sm[0] = 0;
    sm[1] = all;
    r.scan([&](uint64_t id, const bdnnf::Node& nd) {
      if (dead[id]) { sm[id] = 0; return; }
      uint64_t m;
      if (nd.kind == bdnnf::K_AND) {
        m = all & ~killed(nd.units);
        for (uint64_t c : nd.kids) m &= sm[c];
      } else {
        m = 0;
        for (const auto& a : nd.arcs)
          if (!dead[a.child] && !conflicting(a.lits)) m |= sm[a.child] & ~killed(a.lits);
      }
      sm[id] = m;
    });
    uint64_t res = dead[r.root] ? 0 : sm[r.root];
    for (int l : assume) res &= ~kill[lidx(l)];
    return res;
  };
  auto sat_pass = [&](const std::vector<std::vector<int>>& sets) {
    const size_t k = sets.size();
    const uint64_t all = k == 64 ? ~0ULL : (1ULL << k) - 1;
    std::fill(kill.begin(), kill.end(), 0);
    for (size_t j = 0; j < k; j++) for (int l : sets[j]) kill[lidx(-l)] |= 1ULL << j;
    uint64_t res = scan_pass(all);
    for (size_t j = 0; j < k; j++)
      for (int l : sets[j]) if ((kill[lidx(l)] >> j) & 1) res &= ~(1ULL << j);
    return res;
  };
  // a tuple with a literal in no model is unsatisfiable, no pass needed
  auto has_dead = [&](const Combo& c, uint32_t pol) {
    for (uint32_t i = 0; i < t; i++)
      if (lit_dead[lidx((pol >> i) & 1 ? -(int)P[c.c[i]] : (int)P[c.c[i]])]) return true;
    return false;
  };
  auto tuple = [&](const Combo& c, uint32_t pol) {
    std::vector<int> ls;
    for (uint32_t i = 0; i < t; i++) ls.push_back((pol >> i) & 1 ? -(int)P[c.c[i]] : (int)P[c.c[i]]);
    return ls;
  };
  auto compatible = [](const std::vector<int>& a, const std::vector<int>& b) {
    for (int x : a) for (int y : b) if (x == -y) return false;
    return true;
  };
  // Tests the next up to 64 tuples of unknown status from (c0, rank0, pol0)
  // on their own in one pass, so the unsatisfiable ones cost one pass in all.
  auto prefilter = [&](const Combo& c0, uint64_t rank0, uint32_t pol0) {
    Combo c = c0;
    uint64_t rk = rank0;
    uint32_t pl = pol0;
    std::vector<std::vector<int>> sets;
    std::vector<uint64_t> ids;
    for (uint64_t scanned = 0; sets.size() < 64 && scanned < (1u << 16); scanned++) {
      const uint64_t ci = (rk << t) | pl;
      if (!bit(infeasible, ci) && !bit(feasible, ci) && !covered(c, pl)) {
        if (has_dead(c, pl)) setbit(infeasible, ci);
        else {
          sets.push_back(tuple(c, pl));
          ids.push_back(ci);
        }
      }
      if (++pl == (1U << t)) {
        pl = 0;
        rk++;
        if (!c.next()) break;
      }
    }
    if (sets.empty()) return;
    const uint64_t res = sat_pass(sets);
    for (size_t j = 0; j < ids.size(); j++) setbit((res >> j) & 1 ? feasible : infeasible, ids[j]);
  };
  std::vector<int32_t> pidx(nv + 1, -1);
  for (uint32_t i = 0; i < p; i++) pidx[P[i]] = (int32_t)i;
  std::vector<uint64_t> reach(dead.size()), posm(2 * (size_t)nv + 2), seen(nv + 1);
  auto plit = [&](size_t li) { return li & 1 ? -(int)P[li >> 1] : (int)P[li >> 1]; };
  auto pli = [&](int l) { return 2 * (size_t)pidx[std::abs(l)] + (l < 0 ? 1 : 0); };
  // After scan_pass gave the satisfiable lanes res: which literals occur in a
  // model of each, posm per mentioned literal, seen per var (an unseen var is free)
  auto back_pass = [&](uint64_t res) {
    std::fill(reach.begin(), reach.end(), 0);
    std::fill(posm.begin(), posm.end(), 0);
    std::fill(seen.begin(), seen.end(), 0);
    if (dead[r.root]) return res;
    reach[r.root] = res;
    r.scan_backward([&](uint64_t id, const bdnnf::Node& nd) {
      const uint64_t m = reach[id];
      if (!m) return;
      if (nd.kind == bdnnf::K_AND) {
        for (int l : nd.units) { posm[lidx(l)] |= m; seen[std::abs(l)] |= m; }
        for (uint64_t c : nd.kids) reach[c] |= m;
        return;
      }
      for (const auto& a : nd.arcs) {
        for (int l : a.lits) seen[std::abs(l)] |= m;
        if (dead[a.child] || conflicting(a.lits)) continue;
        const uint64_t am = m & sm[a.child] & ~killed(a.lits);
        if (!am) continue;
        for (int l : a.lits) posm[lidx(l)] |= am;
        reach[a.child] |= am;
      }
    });
    return res;
  };
  // lane 0 is A, lane i + 1 is A plus picks[0..i] (no two of them clash)
  auto poss_pass = [&](const std::vector<int>& A, const std::vector<size_t>& picks) {
    const uint64_t all = picks.size() == 63 ? ~0ULL : (2ULL << picks.size()) - 1;
    std::fill(kill.begin(), kill.end(), 0);
    for (int l : A) kill[lidx(-l)] = all;
    for (size_t i = 0; i < picks.size(); i++) kill[lidx(-plit(picks[i]))] = all & ~((2ULL << i) - 1);
    return back_pass(scan_pass(all));
  };
  // literal li of P is in a model of set j (after poss_pass)
  auto possible = [&](size_t li, uint32_t j) {
    const int l = plit(li);
    if ((kill[lidx(l)] >> j) & 1) return false;
    return ((posm[lidx(l)] >> j) & 1) || !((seen[std::abs(l)] >> j) & 1);
  };
  std::vector<int32_t> gain(np2), tg(np2);
  std::vector<char> inA(np2);
  std::vector<uint64_t> cs(cw);
  NodeAccess acc(r);
  bdnnf::Node nd;
  std::mt19937_64 rng(seed);
  std::vector<signed char> val(nv + 1);
  std::vector<uint64_t> stk;
  std::vector<size_t> ok_arcs;
  uint64_t nsamples = 0;
  // grows A greedily; false if A itself is unsatisfiable
  auto grow2 = [&](std::vector<int>& A) {
    std::fill(gain.begin(), gain.end(), 0);
    std::fill(inA.begin(), inA.end(), 0);
    auto add_gain = [&](std::vector<int32_t>& g, size_t a) {
      const uint64_t* row = &cov[a * cw];
      for (size_t w = 0; w < cw; w++) {
        uint64_t x = ~row[w];
        if (w + 1 == cw && (np2 & 63)) x &= (1ULL << (np2 & 63)) - 1;
        if (x == ~0ULL) { for (size_t b = w * 64; b < w * 64 + 64; b++) g[b]++; continue; }
        for (; x; x &= x - 1) g[w * 64 + (size_t)__builtin_ctzll(x)]++;
      }
    };
    auto addA = [&](int l) {
      const size_t li = pli(l);
      inA[li] = inA[li ^ 1] = 1;
      add_gain(gain, li);
    };
    for (int l : A) addA(l);
    std::vector<size_t> picks, cand;
    int lane = -1;  // lane of the last pass that holds A, -1 before the first
    while (true) {
      cand.clear();
      for (size_t b = 0; b < np2; b++) {
        if (inA[b] || gain[b] == 0) continue;
        if (lane >= 0 ? !possible(b, (uint32_t)lane) : lit_dead[lidx(plit(b))]) continue;
        cand.push_back(b);
        tg[b] = gain[b];
      }
      picks.clear();
      while (picks.size() < 63 && !cand.empty()) {
        size_t bi = 0;
        int32_t bg = 0;
        uint64_t ties = 0;
        for (size_t i = 0; i < cand.size(); i++) {
          const int32_t g = tg[cand[i]];
          if (g > bg) { bg = g; bi = i; ties = 1; }
          else if (g == bg && rng() % ++ties == 0) bi = i;
        }
        if (bg == 0) break;
        const size_t best = cand[bi];
        picks.push_back(best);
        const uint64_t* row = &cov[best * cw];
        size_t o = 0;
        for (size_t b : cand) {
          if (b == best || b == (best ^ 1)) continue;
          tg[b] += !((row[b >> 6] >> (b & 63)) & 1);
          cand[o++] = b;
        }
        cand.resize(o);
      }
      if (lane >= 0 && picks.empty()) {
        // a literal may be possible although its var is seen only with the
        // other sign (projected ORs need not be smooth): probe the best ones
        cand.clear();
        for (size_t b = 0; b < np2; b++)
          if (!inA[b] && gain[b] > 0 && !lit_dead[lidx(plit(b))] && !possible(b, (uint32_t)lane))
            cand.push_back(b);
        if (cand.empty()) break;
        const size_t m = std::min<size_t>(64, cand.size());
        std::partial_sort(cand.begin(), cand.begin() + (long)m, cand.end(),
                          [&](size_t x, size_t y) { return gain[x] > gain[y]; });
        std::vector<std::vector<int>> sets(m, A);
        for (size_t i = 0; i < m; i++) sets[i].push_back(plit(cand[i]));
        const uint64_t res = sat_pass(sets);
        if (!res) break;
        A.push_back(plit(cand[(size_t)__builtin_ctzll(res)]));
        addA(A.back());
        poss_pass(A, {});
        lane = 0;
        continue;
      }
      const uint64_t res = poss_pass(A, picks);
      if (!(res & 1)) {
        return false;
      }
      size_t k = 0;
      while (k < picks.size() && ((res >> (k + 1)) & 1)) k++;
      for (size_t i = 0; i < k; i++) {
        A.push_back(plit(picks[i]));
        addA(A.back());
      }
      lane = (int)k;
      if (picks.empty()) break;
    }
    return true;
  };
  // completes the satisfiable set A to a sample and writes it
  auto emit = [&](const std::vector<int>& A) {
    if (!(sat_pass({A}) & 1)) die("internal: t-wise set became unsatisfiable");
    std::fill(val.begin(), val.end(), 0);
    // t = 2: arcs and free vars are chosen by the uncovered pairs they make
    // with the literals set so far (cs), ties at random
    std::fill(cs.begin(), cs.end(), 0);
    auto setl = [&](int l) {
      if (val[std::abs(l)] == (l > 0 ? -1 : 1)) die("FAIL: t-wise walk has conflicting lits");
      val[std::abs(l)] = l > 0 ? 1 : -1;
      if (greedy2 && pidx[std::abs(l)] >= 0) cs[pli(l) >> 6] |= 1ULL << (pli(l) & 63);
    };
    auto score = [&](int l) {
      if (!greedy2 || pidx[std::abs(l)] < 0 || val[std::abs(l)]) return (uint64_t)0;
      const uint64_t* row = &cov[pli(l) * cw];
      uint64_t n = 0;
      for (size_t w = 0; w < cw; w++) n += (uint64_t)__builtin_popcountll(cs[w] & ~row[w]);
      return n;
    };
    if (greedy2) for (int l : A) setl(l);
    stk.assign(1, r.root);
    while (!stk.empty()) {
      const uint64_t x = stk.back();
      stk.pop_back();
      if (x <= 1) { if (x == 0) die("FAIL: t-wise walk reached FALSE"); continue; }
      acc.get(x, nd);
      if (nd.kind == bdnnf::K_AND) {
        for (int l : nd.units) setl(l);
        for (uint64_t k : nd.kids) stk.push_back(k);
        continue;
      }
      ok_arcs.clear();
      for (size_t ai = 0; ai < nd.arcs.size(); ai++) {
        const auto& a = nd.arcs[ai];
        if (!dead[a.child] && !conflicting(a.lits) && (sm[a.child] & 1) && !(killed(a.lits) & 1)) ok_arcs.push_back(ai);
      }
      if (ok_arcs.empty()) die("FAIL: t-wise walk found no satisfiable arc");
      size_t pick = ok_arcs[rng() % ok_arcs.size()];
      if (greedy2 && ok_arcs.size() > 1) {
        uint64_t bs = 0, ties = 0;
        for (size_t ai : ok_arcs) {
          uint64_t sc = 0;
          for (int l : nd.arcs[ai].lits) sc += score(l);
          if (sc > bs || ties == 0) { bs = sc; pick = ai; ties = 1; }
          else if (sc == bs && rng() % ++ties == 0) pick = ai;
        }
      }
      const auto& a = nd.arcs[pick];
      for (int l : a.lits) setl(l);
      stk.push_back(a.child);
    }
    for (uint32_t v = 1; v <= nv; v++) {
      if (val[v] || defined[v]) continue;
      if (!in_proj[v]) { val[v] = -1; continue; }
      if (kill[lidx((int)v)] & 1) val[v] = -1;
      else if (kill[lidx(-(int)v)] & 1) val[v] = 1;
      else if (greedy2) {
        const uint64_t sp = score((int)v), sn = score(-(int)v);
        setl(sp > sn || (sp == sn && (rng() & 1)) ? (int)v : -(int)v);
      } else val[v] = (rng() & 1) ? 1 : -1;
    }
    bdnnf::eval_defs(r.defs, nv, val);
    for (const auto& cl : cls) {
      bool ok = false;
      for (int l : cl) ok |= (val[std::abs(l)] == (l > 0 ? 1 : -1));
      if (!ok) die("FAIL: t-wise sample " + std::to_string(nsamples) + " violates a clause");
    }
    if (greedy2) {
      std::fill(cs.begin(), cs.end(), 0);
      for (size_t i = 0; i < p; i++) setl(val[P[i]] < 0 ? -(int)P[i] : (int)P[i]);
      for (size_t i = 0; i < p; i++) {
        uint64_t* row = &cov[(2 * i + (val[P[i]] < 0 ? 1 : 0)) * cw];
        for (size_t w = 0; w < cw; w++) row[w] |= cs[w];
      }
    } else for (size_t i = 0; i < p; i++) {
      auto& m = smask[2 * i + (val[P[i]] < 0 ? 1 : 0)];
      auto& o = smask[2 * i + (val[P[i]] < 0 ? 0 : 1)];
      m.resize((nsamples >> 6) + 1, 0);
      o.resize((nsamples >> 6) + 1, 0);
      m[nsamples >> 6] |= 1ULL << (nsamples & 63);
    }
    std::string line;
    for (uint32_t v = 1; v <= nv; v++) {
      if (!full && !in_proj[v]) continue;
      line += (val[v] < 0 ? "-" : "") + std::to_string(v) + " ";
    }
    out << line << "0\n";
    nsamples++;
  };
  Combo it(t, p);
  for (uint64_t rank = 0; rank < ncombo; rank++, it.next()) {
    for (uint32_t pol = 0; pol < (1U << t); pol++) {
      const uint64_t idx = (rank << t) | pol;
      if (bit(infeasible, idx) || covered(it, pol)) continue;
      if (has_dead(it, pol)) {
        setbit(infeasible, idx);
        continue;
      }
      if (!bit(feasible, idx)) prefilter(it, rank, pol);
      if (bit(infeasible, idx)) continue;
      std::vector<int> A = tuple(it, pol);
      if (greedy2) {
        if (!grow2(A)) setbit(infeasible, idx);
        else emit(A);
        if (!bit(infeasible, idx) && !covered(it, pol)) die("internal: t-wise sample misses its own tuple");
        continue;
      }
      // the cursor walks the later tuples; an unsatisfiable candidate stays so
      Combo cur = it;
      uint64_t crank = rank;
      uint32_t cpol = pol;
      auto advance = [&]() {
        if (++cpol < (1U << t)) return true;
        cpol = 0;
        crank++;
        return cur.next();
      };
      bool more = advance();
      for (int round = 0; round < 16; round++) {
        std::vector<std::vector<int>> sets{A};
        std::vector<std::tuple<Combo, uint64_t, uint32_t>> pos;
        while (more && sets.size() < 64) {
          const uint64_t ci = (crank << t) | cpol;
          if (!bit(infeasible, ci) && !has_dead(cur, cpol) && !covered(cur, cpol)) {
            const auto U = tuple(cur, cpol);
            if (compatible(A, U)) {
              std::vector<int> AU = A;
              for (int l : U) if (std::find(AU.begin(), AU.end(), l) == AU.end()) AU.push_back(l);
              sets.push_back(std::move(AU));
              pos.emplace_back(cur, crank, cpol);
            }
          }
          more = advance();
        }
        const uint64_t res = sat_pass(sets);
        if (!(res & 1)) {
          setbit(infeasible, idx);
          break;
        }
        size_t first = 1;
        while (first < sets.size() && !((res >> first) & 1)) first++;
        if (first == sets.size()) {
          if (!more) break;
          continue;
        }
        A = sets[first];
        std::tie(cur, crank, cpol) = pos[first - 1];
        more = advance();
      }
      if (bit(infeasible, idx)) continue;
      emit(A);
      if (!covered(it, pol)) die("internal: t-wise sample misses its own tuple");
    }
  }
  out.flush();
  std::cerr << "bdnnf twise: " << nsamples << " samples cover every satisfiable " << t
            << "-tuple of the " << p << " projection vars" << std::endl;
}

std::vector<std::vector<int>> read_queries(const std::string& path, const std::vector<char>& in_proj) {
  std::ifstream f(path);
  if (!f.good()) die("cannot open " + path);
  std::vector<std::vector<int>> qs;
  std::string line;
  while (std::getline(f, line)) {
    if (!line.empty() && line[0] == 'c') continue;
    std::istringstream ss(line);
    std::vector<int> q;
    long l;
    while (ss >> l && l != 0) {
      if ((size_t)std::abs(l) >= in_proj.size() || !in_proj[(size_t)std::abs(l)])
        die("query literal " + std::to_string(l) + " is not of a projection var");
      if (std::find(q.begin(), q.end(), (int)l) == q.end()) q.push_back((int)l);
    }
    qs.push_back(std::move(q));
  }
  return qs;
}

// Answers the queries (assumption sets) in chunks, one forward scan per chunk:
// --sat with one bit per query and node, 64 per scan; counts with one
// rational per query and node, as many per scan as fit kMaxVals.
// Count queries (assumption sets) as in ddnnife's marking: per node, a mask of
// the queries of the chunk that mention a var under it; only those get their
// own value there, the others read the unconditioned one (`base`). Num is
// mpz_class for unweighted counts, where free and smoothing vars are shifts.
template <class Num, class WFn>
void batch_counts(bdnnf::Reader& r, const std::vector<char>& dead, const std::vector<uint64_t>& vs,
                  size_t WD, const std::vector<char>& in_proj, WFn&& w, const std::vector<Num>& base,
                  const std::vector<std::vector<int>>& qs, std::vector<std::string>& out) {
  constexpr bool unw = std::is_same_v<Num, mpz_class>;
  const uint32_t nv = r.nvars;
  const size_t n = dead.size();
  auto row = [&](uint64_t id) { return vs.data() + id * WD; };
  auto has = [](const uint64_t* s, uint32_t v) { return (s[v >> 6] >> (v & 63)) & 1; };
  auto lidx = [](int l) { return 2 * (size_t)std::abs(l) + (l < 0 ? 1 : 0); };
  auto pv = [&](int l) { return in_proj[(uint32_t)std::abs(l)] != 0; };
  std::vector<mpq_class> wl(unw ? 0 : 2 * (size_t)nv + 2);
  if (!unw) for (uint32_t v = 1; v <= nv; v++) { wl[lidx((int)v)] = w((int)v); wl[lidx(-(int)v)] = w(-(int)v); }
  const mpq_class zero = 0;
  // a node stores values only for the queries touching it, so memory is the
  // marking work and every chunk can take 64
  constexpr size_t K = 64;
  std::vector<uint64_t> kill(2 * (size_t)nv + 2), qm(nv + 1), touch(n);
  std::vector<uint32_t> slot(n);
  std::vector<uint8_t> nconf(n, 0);  // per node: 0 unknown, 1 no conflicting arc, 2 some
  mpz_class tmpz;
  std::vector<Num> pool;
  std::vector<uint32_t> gap;
  std::vector<uint64_t> av(WD);
  Num total;  // the unconditioned count, for queries that touch nothing
  {
    Num m = 1;
    if (!dead[r.root])
      for (uint32_t v = 1; v <= nv; v++)
        if (in_proj[v] && !has(row(r.root), v)) {
          if constexpr (unw) m <<= 1; else m *= wl[lidx((int)v)] + wl[lidx(-(int)v)];
        }
    total = dead[r.root] ? Num(0) : Num(base[r.root] * m);
  }
  for (size_t q0 = 0; q0 < qs.size(); q0 += K) {
    const size_t k = std::min(K, qs.size() - q0);
    std::fill(kill.begin(), kill.end(), 0);
    std::fill(qm.begin(), qm.end(), 0);
    uint64_t conflict = 0;
    for (size_t j = 0; j < k; j++)
      for (int l : qs[q0 + j]) {
        kill[lidx(-l)] |= 1ULL << j;
        qm[(uint32_t)std::abs(l)] |= 1ULL << j;
      }
    for (size_t j = 0; j < k; j++)
      for (int l : qs[q0 + j]) if ((kill[lidx(l)] >> j) & 1) conflict |= 1ULL << j;
    std::vector<std::vector<uint32_t>> qv(k);  // each query's distinct vars
    for (size_t j = 0; j < k; j++)
      for (int l : qs[q0 + j]) {
        const uint32_t v = (uint32_t)std::abs(l);
        if (std::find(qv[j].begin(), qv[j].end(), v) == qv[j].end()) qv[j].push_back(v);
      }
    auto wj = [&](int l, size_t j) -> const mpq_class& {
      return (kill[lidx(l)] >> j) & 1 ? zero : wl[lidx(l)];
    };
    auto lit_mask = [&](const std::vector<int>& ls) {
      uint64_t m = 0;
      for (int l : ls) if (pv(l)) m |= qm[(uint32_t)std::abs(l)];
      return m;
    };
    size_t psz = 0;  // pool in use; its numbers are kept across chunks, so their storage is reused
    auto lane = [](uint64_t t, size_t j) { return (size_t)__builtin_popcountll(t & ((1ULL << j) - 1)); };
    auto get = [&](uint64_t c, size_t j) -> const Num& {
      return (touch[c] >> j) & 1 ? pool[slot[c] + lane(touch[c], j)] : base[c];
    };
    // unweighted: the factor of the free/gap vars `vs_` for query j is 2^shift
    auto gap_shift = [&](const std::vector<uint32_t>& vs_, auto&& in_vs, size_t j) {
      unsigned long sh = vs_.size();
      for (uint32_t v : qv[j])
        if (in_vs(v)) {
          BDNNF_CHECK(sh > 0);
          sh--;
        }
      return sh;
    };
    // the factor of the free/gap vars `vs_` (membership: in_vs) for query j
    auto gap_factor = [&](const std::vector<uint32_t>& vs_, auto&& in_vs, size_t j, Num& f) {
      if constexpr (unw) {
        unsigned long sh = vs_.size();
        for (uint32_t v : qv[j])
          if (in_vs(v)) {
            BDNNF_CHECK(sh > 0);
            sh--;
          }
        f = 1;
        f <<= sh;
      } else {
        f = 1;
        for (uint32_t v : vs_) f *= wj((int)v, j) + wj(-(int)v, j);
      }
    };
    touch[0] = touch[1] = 0;
    r.scan([&](uint64_t id, const bdnnf::Node& nd) {
      if (dead[id]) { touch[id] = 0; return; }
      uint64_t t;
      if (nd.kind == bdnnf::K_AND) {
        t = lit_mask(nd.units);
        for (uint64_t c : nd.kids) t |= touch[c];
        gap.clear();
        for (uint32_t v : nd.free) if (in_proj[v]) { gap.push_back(v); t |= qm[v]; }
      } else {
        t = 0;
        for (const auto& a : nd.arcs) if (!dead[a.child]) t |= touch[a.child] | lit_mask(a.lits);
      }
      touch[id] = t;
      if (!t) return;
      BDNNF_CHECK(pool.size() + 64 < UINT32_MAX);
      slot[id] = (uint32_t)psz;
      psz += (size_t)__builtin_popcountll(t);
      if (pool.size() < psz) pool.resize(psz);
      Num f;
      if (nd.kind == bdnnf::K_AND) {
        // once per node: the lanes a unit kills, and the units' weight
        uint64_t uk = 0;
        Num uw = 1;
        for (int l : nd.units) {
          if (!pv(l)) continue;
          uk |= kill[lidx(l)];
          if constexpr (!unw) uw *= wl[lidx(l)];
        }
        for (size_t j = 0; j < k; j++) {
          if (!((t >> j) & 1)) continue;
          Num& x = pool[slot[id] + lane(t, j)];
          if ((uk >> j) & 1) { x = 0; continue; }
          auto in_gap = [&](uint32_t v) { return std::find(gap.begin(), gap.end(), v) != gap.end(); };
          if constexpr (unw) {
            mpz_set_ui(x.get_mpz_t(), 1);
            for (uint64_t c : nd.kids) mpz_mul(x.get_mpz_t(), x.get_mpz_t(), get(c, j).get_mpz_t());
            mpz_mul_2exp(x.get_mpz_t(), x.get_mpz_t(), gap_shift(gap, in_gap, j));
          } else {
            x = 1;
            for (uint64_t c : nd.kids) x *= get(c, j);
            x *= uw;
            gap_factor(gap, in_gap, j, f);
            x *= f;
          }
        }
        return;
      }
      const uint64_t* me = row(id);
      for (size_t i = 0; i < (size_t)__builtin_popcountll(t); i++) pool[slot[id] + i] = 0;
      bool any_conflict = false;
      for (const auto& a : nd.arcs) {
        if (dead[a.child]) continue;
        if (nconf[id] != 1 && has_conflict(a.lits)) {
          any_conflict = true;
          continue;
        }
        const uint64_t* cv = row(a.child);
        av.assign(cv, cv + WD);
        for (int l : a.lits) if (pv(l)) { const uint32_t v = (uint32_t)std::abs(l); av[v >> 6] |= 1ULL << (v & 63); }
        gap.clear();
        for (size_t q = 0; q < WD; q++)
          for (uint64_t m = me[q] & ~av[q]; m; m &= m - 1) gap.push_back((uint32_t)(q * 64 + __builtin_ctzll(m)));
        // once per arc: the lanes its literals kill, and their weight
        uint64_t ak = 0;
        Num aw = 1;
        for (int l : a.lits) {
          if (!pv(l)) continue;
          ak |= kill[lidx(l)];
          if constexpr (!unw) aw *= wl[lidx(l)];
        }
        auto in_gap = [&](uint32_t v) { return has(me, v) && !has(av.data(), v); };
        for (size_t j = 0; j < k; j++) {
          if (!((t >> j) & 1) || ((ak >> j) & 1)) continue;
          Num& acc = pool[slot[id] + lane(t, j)];
          if constexpr (unw) {
            mpz_mul_2exp(tmpz.get_mpz_t(), get(a.child, j).get_mpz_t(), gap_shift(gap, in_gap, j));
            mpz_add(acc.get_mpz_t(), acc.get_mpz_t(), tmpz.get_mpz_t());
          } else {
            Num y = get(a.child, j);
            y *= aw;
            gap_factor(gap, in_gap, j, f);
            acc += y * f;
          }
        }
      }
      if (!nconf[id]) nconf[id] = any_conflict ? 2 : 1;
    });
    std::vector<uint32_t> missing;
    if (!dead[r.root])
      for (uint32_t v = 1; v <= nv; v++) if (in_proj[v] && !has(row(r.root), v)) missing.push_back(v);
    for (size_t j = 0; j < k; j++) {
      Num x;
      if (dead[r.root] || ((conflict >> j) & 1)) x = 0;
      else if (!((touch[r.root] >> j) & 1) && qs[q0 + j].empty()) x = total;
      else {
        Num f;
        gap_factor(missing, [&](uint32_t v) { return !has(row(r.root), v); }, j, f);
        x = get(r.root, j) * f;
      }
      if constexpr (!unw) x.canonicalize();
      out.push_back(x.get_str());
    }
  }
}

// SAT queries (assumption sets, plus --assume), 64 per forward scan with one
// bit per query and node.
void batch_sat(bdnnf::Reader& r, const std::vector<char>& dead, const std::vector<char>& in_proj,
               const std::vector<std::vector<int>>& qs, const std::vector<int>& assume,
               std::vector<std::string>& out) {
  const uint32_t nv = r.nvars;
  const size_t n = dead.size();
  auto conflicting = has_conflict;
  // kill[lit]: the queries of the chunk that assume -lit
  std::vector<uint64_t> kill(2 * (size_t)nv + 2), sm(n);
  auto lidx = [](int l) { return 2 * (size_t)std::abs(l) + (l < 0 ? 1 : 0); };
  for (size_t q0 = 0; q0 < qs.size(); q0 += 64) {
    const size_t k = std::min<size_t>(64, qs.size() - q0);
    const uint64_t all = k == 64 ? ~0ULL : (1ULL << k) - 1;
    std::fill(kill.begin(), kill.end(), 0);
    for (size_t j = 0; j < k; j++)
      for (int l : qs[q0 + j]) kill[lidx(-l)] |= 1ULL << j;
    for (int l : assume) kill[lidx(-l)] |= all;
    auto killed = [&](const std::vector<int>& ls) {
      uint64_t m = 0;
      for (int l : ls) if (in_proj[(uint32_t)std::abs(l)]) m |= kill[lidx(l)];
      return m;
    };
    sm[0] = 0;
    sm[1] = all;
    r.scan([&](uint64_t id, const bdnnf::Node& nd) {
      if (dead[id]) { sm[id] = 0; return; }
      uint64_t m;
      if (nd.kind == bdnnf::K_AND) {
        m = all & ~killed(nd.units);
        for (uint64_t c : nd.kids) m &= sm[c];
      } else {
        m = 0;
        for (const auto& a : nd.arcs)
          if (!dead[a.child] && !conflicting(a.lits)) m |= sm[a.child] & ~killed(a.lits);
      }
      sm[id] = m;
    });
    for (size_t j = 0; j < k; j++) {
      const bool ok = !dead[r.root] && ((sm[r.root] >> j) & 1);
      bool conflict = false;
      for (int l : qs[q0 + j]) conflict |= (kill[lidx(l)] >> j) & 1;
      for (int l : assume) conflict |= (kill[lidx(l)] >> j) & 1;
      out.push_back(ok && !conflict ? "SAT" : "UNSAT");
    }
  }
}


// Checks t-wise samples (DIMACS lines) against the circuit: every t-tuple of
// projection literals no sample covers must be unsatisfiable. Returns the
// number of satisfiable uncovered tuples, after printing a summary.
uint64_t twise_check(bdnnf::Reader& r, const std::vector<char>& dead, const std::vector<char>& in_proj,
                     const std::vector<int>& assume, uint32_t t, const std::string& samples_path,
                     const std::vector<char>& lit_dead) {
  const uint32_t nv = r.nvars;
  std::vector<char> defined(nv + 1, 0);
  for (const auto& [v, rf] : r.defs.defs) defined[v] = 1;
  std::vector<uint32_t> P;
  std::vector<int64_t> pos(nv + 1, -1);
  for (uint32_t v = 1; v <= nv; v++)
    if (in_proj[v] && !defined[v]) { pos[v] = (int64_t)P.size(); P.push_back(v); }
  const uint32_t p = (uint32_t)P.size();
  t = std::min(t, p);
  std::ifstream f(samples_path);
  if (!f.good()) die("cannot open " + samples_path);
  // smask[2i + neg]: the samples containing the literal of P[i]
  std::vector<std::vector<uint64_t>> smask(2 * (size_t)p);
  uint64_t ns = 0;
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == 'c') continue;
    std::istringstream ss(line);
    long l;
    for (auto& m : smask) m.resize((ns >> 6) + 1, 0);
    while (ss >> l && l != 0) {
      const uint64_t v = (uint64_t)std::abs(l);
      if (v > nv) die("sample literal out of range: " + std::to_string(l));
      if (pos[v] >= 0) smask[2 * (size_t)pos[v] + (l < 0 ? 1 : 0)][ns >> 6] |= 1ULL << (ns & 63);
    }
    ns++;
  }
  auto lit_of = [&](const Combo& c, uint32_t pol, uint32_t i) {
    return (pol >> i) & 1 ? -(int)P[c.c[i]] : (int)P[c.c[i]];
  };
  std::vector<std::vector<int>> qs;
  uint64_t bad = 0, unc = 0;
  std::vector<int> example;
  std::vector<std::string> ans;
  auto flush = [&]() {
    if (qs.empty()) return;
    ans.clear();
    batch_sat(r, dead, in_proj, qs, assume, ans);
    for (size_t j = 0; j < qs.size(); j++)
      if (ans[j] == "SAT") {
        if (!bad) example = qs[j];
        bad++;
      }
    qs.clear();
  };
  Combo c(t, p);
  bool more = p > 0;
  while (more) {
    for (uint32_t pol = 0; pol < (1U << t); pol++) {
      bool cov = false;
      for (size_t w = 0; w < smask[0].size() && !cov; w++) {
        uint64_t m = ~0ULL;
        for (uint32_t i = 0; i < t && m; i++) m &= smask[2 * (size_t)c.c[i] + ((pol >> i) & 1)][w];
        cov = m != 0;
      }
      if (cov) continue;
      unc++;
      std::vector<int> q;
      bool dl = false;
      for (uint32_t i = 0; i < t; i++) {
        q.push_back(lit_of(c, pol, i));
        dl |= lit_dead[2 * (size_t)std::abs(q.back()) + (q.back() < 0 ? 1 : 0)] != 0;
      }
      if (dl) continue;
      qs.push_back(std::move(q));
      if (qs.size() == 1024) flush();
    }
    more = c.next();
  }
  flush();
  std::cout << ns << " samples, " << unc << " uncovered " << t << "-tuples, " << bad
            << " of them satisfiable";
  if (bad) {
    std::cout << " (e.g.";
    for (int l : example) std::cout << " " << l;
    std::cout << ")";
  }
  std::cout << std::endl;
  return bad;
}

// argv[0] is the command. `pre`: an open (cached) reader of the file to use.
int run_query(int argc, char** argv, bdnnf::Reader* pre) {
  const std::string cmd = argv[0];
  const bool enumerating = cmd == "enum";
  const bool atomic = cmd == "atomic";
  const bool opt = cmd == "opt";
  const bool tw_check = cmd == "twise-check";
  const bool tw = cmd == "twise" || tw_check;
  const bool sampling = cmd == "sample" || enumerating || atomic;
  const bool lits = cmd == "lits" || cmd == "core";
  const bool core = cmd == "core";
  const bool is_batch = cmd == "batch";
  const std::string usage = tw_check
    ? "usage: bdnnf twise-check <in.bdnnf> <samples> [-t T] [--cnf <file.cnf>] [--assume L1,L2,...]"
    : tw
    ? "usage: bdnnf twise <in.bdnnf> [-t T] [--seed S] [--cnf <file.cnf>] [--assume L1,L2,...] [--full] [-o out]"
    : opt
    ? "usage: bdnnf opt <in.bdnnf> --obj <file> [-k K] [--min] [--full] [--cnf <file.cnf>] [--assume L1,L2,...]"
      "   (obj file: lines of 'literal value')"
    : atomic
    ? "usage: bdnnf atomic <in.bdnnf> [--cross] [--cnf <file.cnf>] [--assume L1,L2,...]"
    : enumerating
    ? "usage: bdnnf enum <in.bdnnf> [-n N] [--start I] [--cnf <file.cnf>] [--assume L1,L2,...] [--full] [-o out]"
    : sampling
    ? "usage: bdnnf sample <in.bdnnf> [-n N] [--seed S] [--cnf <file.cnf>] [--assume L1,L2,...] [--full] [-o out]"
    : is_batch
    ? "usage: bdnnf batch <in.bdnnf> <queries> [--sat] [--cnf <file.cnf>] [--cache GB]   (a query: a line of literals)"
    : "usage: bdnnf " + cmd + " <in.bdnnf> [--cnf <file.cnf>] [--assume L1,L2,...]";
  std::string path, cnf, out_path;
  uint64_t num = enumerating ? UINT64_MAX : atomic ? 64 : 1, seed = 1;
  mpz_class start = 0;
  bool full = false, sat = false, cross = false, minimize = false;
  double cache_gb = 5;  // --cache: load circuits up to this size into memory
  std::string obj_path;
  size_t kbest = 1;
  uint32_t tdeg = 2;
  std::string queries;
  std::vector<int> assume;
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    if (a == "--cnf" && i + 1 < argc) cnf = argv[++i];
    else if (sampling && (a == "-n" || a == "--num") && i + 1 < argc) num = std::stoull(argv[++i]);
    else if (sampling && !enumerating && a == "--seed" && i + 1 < argc) seed = std::stoull(argv[++i]);
    else if (enumerating && a == "--start" && i + 1 < argc) start = mpz_class(argv[++i], 10);
    else if (sampling && a == "--full") full = true;
    else if (sampling && a == "-o" && i + 1 < argc) out_path = argv[++i];
    else if (a == "--cache" && i + 1 < argc) cache_gb = std::stod(argv[++i]);
    else if (a == "--assume" && i + 1 < argc) {
      std::string t = argv[++i];
      for (char& ch : t) if (ch == ',') ch = ' ';
      std::istringstream ss(t);
      long l;
      while (ss >> l) if (l != 0) assume.push_back((int)l);
    }
    else if (is_batch && a == "--sat") sat = true;
    else if (atomic && a == "--cross") cross = true;
    else if (opt && a == "--obj" && i + 1 < argc) obj_path = argv[++i];
    else if (tw && a == "-t" && i + 1 < argc) tdeg = (uint32_t)std::stoul(argv[++i]);
    else if (tw && a == "--seed" && i + 1 < argc) seed = std::stoull(argv[++i]);
    else if (tw && a == "--full") full = true;
    else if (tw && a == "-o" && i + 1 < argc) out_path = argv[++i];
    else if (opt && a == "-k" && i + 1 < argc) kbest = std::stoull(argv[++i]);
    else if (opt && a == "--min") minimize = true;
    else if (opt && a == "--full") full = true;
    else if (path.empty() && (a == "-" || a.empty() || a[0] != '-')) path = a;
    else if ((is_batch || tw_check) && queries.empty() && (a.empty() || a[0] != '-')) queries = a;
    else die(usage);
  }
  if (path.empty() || ((is_batch || tw_check) && queries.empty()) || (opt && (obj_path.empty() || kbest == 0))) die(usage);
  try {
    bdnnf::Reader own;
    bdnnf::Reader& r = pre ? *pre : own;
    if (!pre) r.open(path);
    // queries making many passes keep the circuit in memory if it is small enough
    const bool can_cache = path != "-" && cache_gb > 0 &&
        (double)std::filesystem::file_size(path) <= cache_gb * (1ULL << 30);
    auto maybe_cache = [&]() {
      if (can_cache && !r.is_cached()) r.cache_all();
    };
    const uint32_t nv = r.nvars;
    std::vector<char> file_proj(nv + 1, 1);
    const bool file_projected = r.has(r.header_kv, "projection");
    const std::string hp = r.get(r.header_kv, "projection");
    if (file_projected) {
      file_proj.assign(nv + 1, 0);
      for (uint32_t v : bdnnf::parse_var_list(hp)) {
        if (v == 0 || v > nv) die("projection var out of range: " + std::to_string(v));
        file_proj[v] = 1;
      }
    }
    Weights W;
    if (!cnf.empty()) {
      W = read_cnf(cnf, nv);
      for (uint32_t v = 1; v <= nv; v++)
        if (W.in_proj[v] != file_proj[v])
          die("the CNF's projection differs from the circuit's (var " + std::to_string(v) + ")");
    } else {
      W.projected = file_projected;
      W.in_proj = file_proj;
    }
    for (const auto& [v, rf] : r.defs.defs)
      if (W.in_proj[v]) die("var " + std::to_string(v) + " is in the projection but defined: "
                            "its count is not a d-DNNF query");
    if (core || enumerating || atomic || opt || tw) {
      W.weighted = false;
      W.w.clear();
    }
    for (int l : assume) {
      if ((uint32_t)std::abs(l) > nv || !W.in_proj[(uint32_t)std::abs(l)])
        die("--assume " + std::to_string(l) + ": not a projection var");
      if (!W.weighted) W.w.clear();
      W.weighted = true;
      if (!W.w.count(l)) W.w[l] = 1;
      W.w[-l] = 0;
    }
    if (sampling && !enumerating)
      for (const auto& [l, x] : W.w)
        if (x < 0 && W.in_proj[(uint32_t)std::abs(l)])
          die("negative weight on literal " + std::to_string(l) + ": cannot sample");
    const mpq_class one = 1;
    auto w = [&](int l) -> mpq_class {
      const uint32_t v = (uint32_t)std::abs(l);
      if (!W.in_proj[v] || !W.weighted) return one;
      auto it = W.w.find(l);
      return it == W.w.end() ? one : it->second;
    };
    auto both = [&](uint32_t v) -> mpq_class { return w((int)v) + w(-(int)v); };

    const size_t n = (size_t)r.nnodes + 2;
    const size_t WD = (size_t)nv / 64 + 1;
    std::vector<mpq_class> val(n);
    std::vector<uint64_t> vs(n * WD, 0);
    std::vector<char> dead(n, 0);
    val[0] = 0;
    val[1] = 1;
    dead[0] = 1;
    auto row = [&](uint64_t id) { return vs.data() + id * WD; };
    auto set = [](uint64_t* s, uint32_t v) { s[v >> 6] |= 1ULL << (v & 63); };
    auto has = [](const uint64_t* s, uint32_t v) { return (s[v >> 6] >> (v & 63)) & 1; };
    auto conflicting = has_conflict;
    std::vector<uint64_t> av;
    // weight of an OR arc, smoothed over the vars of its node `me`
    auto arc_val = [&](const uint64_t* me, const bdnnf::Arc& a) -> mpq_class {
      if (dead[a.child] || conflicting(a.lits)) return 0;
      const uint64_t* cv = row(a.child);
      av.assign(cv, cv + WD);
      mpq_class y = val[a.child];
      for (int l : a.lits) {
        y *= w(l);
        if (W.in_proj[(uint32_t)std::abs(l)]) set(av.data(), (uint32_t)std::abs(l));
      }
      for (size_t q = 0; q < WD; q++)
        for (uint64_t m = me[q] & ~av[q]; m; m &= m - 1)
          y *= both((uint32_t)(q * 64 + __builtin_ctzll(m)));
      return y;
    };
    r.scan([&](uint64_t id, const bdnnf::Node& nd) {
      uint64_t* me = row(id);
      if (nd.kind == bdnnf::K_AND) {
        bool d = conflicting(nd.units);
        for (uint64_t k : nd.kids) d |= dead[k] != 0;
        if (d) { dead[id] = 1; val[id] = 0; return; }
        mpq_class x = 1;
        for (uint64_t k : nd.kids) {
          x *= val[k];
          const uint64_t* kv = row(k);
          for (size_t q = 0; q < WD; q++) me[q] |= kv[q];
        }
        for (int l : nd.units) {
          x *= w(l);
          if (W.in_proj[(uint32_t)std::abs(l)]) set(me, (uint32_t)std::abs(l));
        }
        for (uint32_t v : nd.free)
          if (W.in_proj[v]) { x *= both(v); set(me, v); }
        val[id] = std::move(x);
        return;
      }
      bool any = false;
      for (const auto& a : nd.arcs) {
        if (dead[a.child] || conflicting(a.lits)) continue;
        any = true;
        const uint64_t* cv = row(a.child);
        for (size_t q = 0; q < WD; q++) me[q] |= cv[q];
        for (int l : a.lits)
          if (W.in_proj[(uint32_t)std::abs(l)]) set(me, (uint32_t)std::abs(l));
      }
      if (!any) { dead[id] = 1; val[id] = 0; return; }
      mpq_class x = 0;
      for (const auto& a : nd.arcs) x += arc_val(me, a);
      val[id] = std::move(x);
    });
    mpq_class res = val[r.root];
    if (!dead[r.root])
      for (uint32_t v = 1; v <= nv; v++)
        if (W.in_proj[v] && !has(row(r.root), v)) res *= both(v);
    res.canonicalize();
    if (tw) {
      if (res == 0) die("no model");
      if (tdeg == 0) die("-t must be at least 1");
      std::ofstream of;
      if (!out_path.empty() && out_path != "-") {
        of.open(out_path);
        if (!of) die("cannot open " + out_path);
      }
      maybe_cache();
      const auto cnt = lit_counts(r, val, dead, vs, WD, W.in_proj, w);
      std::vector<char> lit_dead(2 * (size_t)nv + 2, 0);
      for (size_t i = 2; i < lit_dead.size(); i++) lit_dead[i] = cnt[i] == 0;
      if (tw_check) return twise_check(r, dead, W.in_proj, assume, tdeg, queries, lit_dead) ? EXIT_FAILURE : EXIT_SUCCESS;
      twise(r, dead, W.in_proj, assume, tdeg, seed, full,
            cnf.empty() ? std::vector<std::vector<int>>{} : read_clauses(cnf), lit_dead,
            of.is_open() ? of : std::cout);
      return EXIT_SUCCESS;
    }
    if (opt) {
      if (res == 0) die("no model");
      std::map<int, mpq_class> obj = read_objective(obj_path, nv);
      if (minimize) for (auto& [l, x] : obj) x = -x;
      auto ov = [&](int l) -> mpq_class {
        if (!W.in_proj[(uint32_t)std::abs(l)]) return 0;
        auto it = obj.find(l);
        return it == obj.end() ? mpq_class(0) : it->second;
      };
      auto var_factor = [&](uint32_t v) {
        KFactor f;
        for (int l : {(int)v, -(int)v}) if (w(l) != 0) f.lits.push_back(l);
        if (f.lits.size() == 2 && ov(f.lits[1]) > ov(f.lits[0])) std::swap(f.lits[0], f.lits[1]);
        for (int l : f.lits) f.vals.push_back(ov(l));
        return f;
      };
      std::vector<std::vector<mpq_class>> kb(n);
      auto list_of = [&](const KStage& st) {
        std::vector<mpq_class> l;
        for (const auto& it : st) l.push_back(it.v);
        return l;
      };
      std::vector<KFactor> fs;
      std::vector<KStage> st;
      auto and_factors = [&](const bdnnf::Node& nd, mpq_class& base) {
        fs.clear();
        base = 0;
        for (int l : nd.units) base += ov(l);
        for (uint64_t c : nd.kids) { KFactor f; f.child = c; f.vals = kb[c]; fs.push_back(std::move(f)); }
        for (uint32_t v : nd.free) if (W.in_proj[v]) fs.push_back(var_factor(v));
      };
      auto arc_factors = [&](const uint64_t* me, const bdnnf::Arc& a, mpq_class& base) {
        fs.clear();
        base = 0;
        for (int l : a.lits) base += ov(l);
        KFactor f;
        f.child = a.child;
        f.vals = kb[a.child];
        fs.push_back(std::move(f));
        const uint64_t* cv = row(a.child);
        av.assign(cv, cv + WD);
        for (int l : a.lits) {
          const uint32_t v = (uint32_t)std::abs(l);
          if (W.in_proj[v]) set(av.data(), v);
        }
        for (size_t q = 0; q < WD; q++)
          for (uint64_t m = me[q] & ~av[q]; m; m &= m - 1) fs.push_back(var_factor((uint32_t)(q * 64 + __builtin_ctzll(m))));
      };
      // the arcs' k-best lists merged: (value, arc, rank in the arc's list)
      auto or_merge = [&](uint64_t id, const bdnnf::Node& nd, KStage& out) {
        out.clear();
        mpq_class base;
        for (uint32_t ai = 0; ai < nd.arcs.size(); ai++) {
          if (arc_val(row(id), nd.arcs[ai]) == 0) continue;
          arc_factors(row(id), nd.arcs[ai], base);
          kseq(fs, base, kbest, st);
          for (uint32_t j = 0; j < st.back().size(); j++) out.push_back({st.back()[j].v, ai, j});
        }
        const size_t m = std::min(kbest, out.size());
        std::partial_sort(out.begin(), out.begin() + (long)m, out.end(), [](const KItem& x, const KItem& y) {
          return x.v != y.v ? x.v > y.v : (x.a != y.a ? x.a < y.a : x.b < y.b);
        });
        out.resize(m);
      };
      kb[1] = {0};
      KStage merged;
      r.scan([&](uint64_t id, const bdnnf::Node& nd) {
        if (dead[id] || val[id] == 0) return;
        if (nd.kind == bdnnf::K_AND) {
          mpq_class base;
          and_factors(nd, base);
          kseq(fs, base, kbest, st);
          kb[id] = list_of(st.back());
          return;
        }
        or_merge(id, nd, merged);
        kb[id] = list_of(merged);
      });
      std::vector<KFactor> root_fs(1);
      root_fs[0].child = r.root;
      root_fs[0].vals = kb[r.root];
      for (uint32_t v = 1; v <= nv; v++)
        if (W.in_proj[v] && !has(row(r.root), v)) root_fs.push_back(var_factor(v));
      std::vector<KStage> root_st;
      kseq(root_fs, 0, kbest, root_st);
      const auto cls = cnf.empty() ? std::vector<std::vector<int>>{} : read_clauses(cnf);
      NodeAccess acc(r);
      bdnnf::Node nd;
      std::vector<signed char> val2(nv + 1);
      std::vector<char> defined(nv + 1, 0);
      for (const auto& [v, rf] : r.defs.defs) defined[v] = 1;
      std::vector<std::pair<uint64_t, uint32_t>> stk;
      std::string line;
      auto setl = [&](int l) {
        if (val2[std::abs(l)] == (l > 0 ? -1 : 1)) die("FAIL: conflicting lits on the path");
        val2[std::abs(l)] = l > 0 ? 1 : -1;
      };
      // walk the stages back from rank j: each factor's rank, then recurse
      auto unwind = [&](const std::vector<KFactor>& f, const std::vector<KStage>& sts, uint32_t j) {
        for (size_t s = f.size(); s > 0; s--) {
          BDNNF_CHECK(j < sts[s].size());
          const KItem& it = sts[s][j];
          if (f[s - 1].child) stk.push_back({f[s - 1].child, it.b});
          else setl(f[s - 1].lits[it.b]);
          j = it.a;
        }
      };
      for (uint32_t rank = 0; rank < root_st.back().size(); rank++) {
        std::fill(val2.begin(), val2.end(), 0);
        stk.clear();
        unwind(root_fs, root_st, rank);
        while (!stk.empty()) {
          const auto [x, j] = stk.back();
          stk.pop_back();
          if (x <= 1) { if (x == 0) die("FAIL: reached FALSE"); continue; }
          acc.get(x, nd);
          mpq_class base;
          if (nd.kind == bdnnf::K_AND) {
            for (int l : nd.units) setl(l);
            and_factors(nd, base);
            kseq(fs, base, kbest, st);
            const auto f = fs;
            unwind(f, st, j);
            continue;
          }
          or_merge(x, nd, merged);
          BDNNF_CHECK(j < merged.size());
          const KItem pick = merged[j];
          const auto& a = nd.arcs[pick.a];
          for (int l : a.lits) setl(l);
          arc_factors(row(x), a, base);
          kseq(fs, base, kbest, st);
          const auto f = fs;
          unwind(f, st, pick.b);
        }
        for (uint32_t v = 1; v <= nv; v++) {
          if (val2[v] || defined[v]) continue;
          if (W.in_proj[v]) die("FAIL: optimum leaves var " + std::to_string(v) + " open");
          val2[v] = -1;
        }
        bdnnf::eval_defs(r.defs, nv, val2);
        for (const auto& cl : cls) {
          bool ok = false;
          for (int l : cl) ok |= (val2[std::abs(l)] == (l > 0 ? 1 : -1));
          if (!ok) die("FAIL: optimal model " + std::to_string(rank) + " violates a clause");
        }
        mpq_class value = root_st.back()[rank].v;
        if (minimize) value = -value;
        value.canonicalize();
        line = "c value " + value.get_str() + "\n";
        for (uint32_t v = 1; v <= nv; v++) {
          if (!full && !W.in_proj[v]) continue;
          line += (val2[v] < 0 ? "-" : "") + std::to_string(v) + " ";
        }
        std::cout << line << "0\n";
      }
      std::cout << std::flush;
      return EXIT_SUCCESS;
    }
    if (is_batch) {
      std::vector<std::string> ans;
      const auto qs = read_queries(queries, W.in_proj);
      if (qs.size() > 64) maybe_cache();
      if (sat) batch_sat(r, dead, W.in_proj, qs, assume, ans);
      else if (!W.weighted) {
        std::vector<mpz_class> base(val.size());
        for (size_t i = 0; i < val.size(); i++) base[i] = val[i].get_num();
        batch_counts<mpz_class>(r, dead, vs, WD, W.in_proj, w, base, qs, ans);
      } else batch_counts<mpq_class>(r, dead, vs, WD, W.in_proj, w, val, qs, ans);
      BDNNF_CHECK(ans.size() == qs.size());
      std::string o;
      for (const auto& a : ans) { o += a; o += '\n'; }
      std::cout << o << std::flush;
      return EXIT_SUCCESS;
    }
    if (lits) {
      const auto cnt = lit_counts(r, val, dead, vs, WD, W.in_proj, w);
      std::string line;
      if (core) {
        if (res == 0) die("no model");
        for (uint32_t v = 1; v <= nv; v++) {
          if (!W.in_proj[v]) continue;
          if (cnt[2 * v] == res) line += std::to_string(v) + " ";
          else if (cnt[2 * v + 1] == res) line += "-" + std::to_string(v) + " ";
        }
        std::cout << line << "0" << std::endl;
        return EXIT_SUCCESS;
      }
      for (uint32_t v = 1; v <= nv; v++) {
        if (!W.in_proj[v]) continue;
        mpq_class a = cnt[2 * v], b = cnt[2 * v + 1];
        a.canonicalize();
        b.canonicalize();
        line = std::to_string(v) + " " + a.get_str() + "\n-" + std::to_string(v) + " " + b.get_str() + "\n";
        std::cout << line;
      }
      std::cout << std::flush;
      return EXIT_SUCCESS;
    }
    if (!sampling) {
      std::cout << res.get_str() << std::endl;
      return EXIT_SUCCESS;
    }
    if (res == 0 && !enumerating) die(atomic ? "no model" : "no model of non-zero weight");
    if (enumerating) {
      const mpz_class left = start < res.get_num() ? mpz_class(res.get_num() - start) : mpz_class(0);
      if (left < mpz_u64(num)) num = u64_mpz(left);
    }

    const auto cls = cnf.empty() ? std::vector<std::vector<int>>{} : read_clauses(cnf);
    NodeAccess acc(r);
    gmp_randclass rnd(gmp_randinit_default);
    rnd.seed(seed);
    // Per node visited, decoded once; per OR node, its arc weights, and for
    // sampling scaled to integers with prefix sums.
    struct OrW {
      std::vector<mpq_class> ys;
      std::vector<mpz_class> cum;
    };
    maybe_cache();
    // uncached: the nodes visited, decoded, and the OR nodes' weights, each
    // bounded (cleared when full)
    constexpr size_t kMaxCached = 1 << 20;
    std::unordered_map<uint64_t, bdnnf::Node> ncache;
    std::unordered_map<uint64_t, OrW> orw;
    auto node_of = [&](uint64_t x) -> const bdnnf::Node& {
      if (r.is_cached()) return r.cache[x - 2];
      if (ncache.size() >= kMaxCached) ncache.clear();
      auto it = ncache.find(x);
      if (it != ncache.end()) return it->second;
      bdnnf::Node dec;
      acc.get(x, dec);
      return ncache.emplace(x, std::move(dec)).first->second;
    };
    // cumulative integer weights of ys, the same ratios
    auto cumulate = [](const std::vector<mpq_class>& ys, std::vector<mpz_class>& cum) {
      mpz_class den = 1;
      for (const auto& y : ys) mpz_lcm(den.get_mpz_t(), den.get_mpz_t(), y.get_den_mpz_t());
      cum.clear();
      mpz_class acc_ = 0;
      for (const auto& y : ys) {
        acc_ += y.get_num() * (den / y.get_den());
        cum.push_back(acc_);
      }
      if (cum.empty() || cum.back() <= 0) die("sampled into a node of weight 0");
    };
    // index i with probability (cum[i] - cum[i-1]) / cum.back()
    auto pick = [&](const std::vector<mpz_class>& cum) -> size_t {
      const mpz_class p = rnd.get_z_range(cum.back());
      const size_t i = (size_t)(std::upper_bound(cum.begin(), cum.end(), p) - cum.begin());
      BDNNF_CHECK(i < cum.size());
      return i;
    };
    auto or_weights = [&](uint64_t x, const bdnnf::Node& nd) -> const OrW& {
      if (orw.size() >= kMaxCached) orw.clear();
      auto it = orw.find(x);
      if (it != orw.end()) return it->second;
      OrW o;
      for (const auto& a : nd.arcs) o.ys.push_back(arc_val(row(x), a));
      if (!enumerating) cumulate(o.ys, o.cum);
      return orw.emplace(x, std::move(o)).first->second;
    };
    std::vector<std::vector<mpz_class>> var_cum(nv + 1);  // free projection var: w(-v), w(v)
    std::vector<char> defined(nv + 1, 0);
    for (const auto& [v, rf] : r.defs.defs) defined[v] = 1;
    std::vector<signed char> val2(nv + 1);
    std::vector<char> rnd_set(nv + 1);
    std::vector<uint64_t> stk, path_nodes;
    std::vector<mpz_class> istk;
    std::vector<uint64_t> cav(WD);
    std::vector<uint64_t> sig(atomic ? nv + 1 : 0);  // atomic: var values in the samples
    // Unranking: the models of a node are numbered in mixed radix over its
    // factors (children, then free vars, each var a digit of radix w(v)+w(-v)).
    auto digit = [](mpz_class& i, const mpz_class& radix) {
      mpz_class d;
      mpz_fdiv_qr(i.get_mpz_t(), d.get_mpz_t(), i.get_mpz_t(), radix.get_mpz_t());
      return d;
    };
    std::ofstream of;
    if (!out_path.empty() && out_path != "-") {
      of.open(out_path);
      if (!of) die("cannot open " + out_path);
    }
    std::ostream& out = of.is_open() ? of : std::cout;
    std::string line;
    for (uint64_t smp = 0; smp < num; smp++) {
      std::fill(val2.begin(), val2.end(), 0);
      std::fill(rnd_set.begin(), rnd_set.end(), 0);
      auto setl = [&](int l) {
        const signed char want = l > 0 ? 1 : -1;
        if (val2[std::abs(l)] == -want) die("FAIL: sampled path has conflicting lits");
        val2[std::abs(l)] = want;
      };
      // a free var of the unranked model: radix 2, or 1 when one literal has weight 0
      auto set_free = [&](uint32_t v, mpz_class& i) {
        const mpz_class rad = both(v).get_num();
        BDNNF_CHECK(rad == 1 || rad == 2);
        setl(rad == 1 ? (w((int)v) == 0 ? -(int)v : (int)v) : (digit(i, rad) == 1 ? (int)v : -(int)v));
      };
      stk.assign(1, r.root);
      istk.clear();
      if (enumerating) {
        mpz_class i = start + mpz_u64(smp);
        for (uint32_t v = 1; v <= nv; v++)
          if (W.in_proj[v] && !has(row(r.root), v)) set_free(v, i);
        istk.push_back(i);
      }
      path_nodes.clear();
      while (!stk.empty()) {
        const uint64_t x = stk.back();
        stk.pop_back();
        mpz_class i;
        if (enumerating) { i = istk.back(); istk.pop_back(); }
        path_nodes.push_back(x);
        if (x <= 1) {
          if (x == 0) die("FAIL: sampled into FALSE");
          if (enumerating) BDNNF_CHECK(i == 0);
          continue;
        }
        const bdnnf::Node& nd = node_of(x);
        if (nd.kind == bdnnf::K_AND) {
          for (int l : nd.units) setl(l);
          for (uint64_t k : nd.kids) {
            stk.push_back(k);
            if (enumerating) istk.push_back(digit(i, val[k].get_num()));
          }
          if (enumerating) for (uint32_t v : nd.free) if (W.in_proj[v]) set_free(v, i);
          continue;
        }
        const OrW& ow = or_weights(x, nd);
        const auto& ys = ow.ys;
        size_t ai = 0;
        if (!enumerating) ai = pick(ow.cum);
        else
          while (i >= ys[ai].get_num()) {
            i -= ys[ai].get_num();
            ai++;
            BDNNF_CHECK(ai < ys.size());
          }
        const auto& a = nd.arcs[ai];
        for (int l : a.lits) setl(l);
        stk.push_back(a.child);
        if (enumerating) {
          istk.push_back(digit(i, val[a.child].get_num()));
          const uint64_t* cv = row(a.child);
          cav.assign(cv, cv + WD);
          for (int l : a.lits) {
            const uint32_t v = (uint32_t)std::abs(l);
            if (W.in_proj[v]) set(cav.data(), v);
          }
          const uint64_t* me = row(x);
          for (size_t q = 0; q < WD; q++)
            for (uint64_t m = me[q] & ~cav[q]; m; m &= m - 1) set_free((uint32_t)(q * 64 + __builtin_ctzll(m)), i);
        }
      }
      // Projection vars off the path are free: each is drawn by its weight.
      for (uint32_t v = 1; v <= nv; v++) {
        if (val2[v] || defined[v]) continue;
        rnd_set[v] = 1;
        if (enumerating) {
          if (W.in_proj[v]) die("FAIL: model " + mpz_class(start + mpz_u64(smp)).get_str() + " leaves var " + std::to_string(v) + " open");
          val2[v] = -1;
        } else if (W.in_proj[v]) {
          if (var_cum[v].empty()) cumulate({w(-(int)v), w((int)v)}, var_cum[v]);
          val2[v] = pick(var_cum[v]) == 1 ? 1 : -1;
        } else val2[v] = rnd.get_z_bits(1) == 1 ? 1 : -1;
      }
      bdnnf::eval_defs(r.defs, nv, val2);
      for (const auto& cl : cls) {
        bool ok = false;
        for (int l : cl) ok |= (val2[std::abs(l)] == (l > 0 ? 1 : -1));
        if (ok) continue;
        std::ostringstream m;
        m << "FAIL: sample " << smp << " violates clause";
        for (int l : cl) m << " " << l << (rnd_set[std::abs(l)] ? "(not on path)" : "");
        m << "; nodes visited:";
        for (uint64_t x : path_nodes) m << " " << x;
        die(m.str());
      }
      if (atomic) {
        for (uint32_t v = 1; v <= nv; v++) if (val2[v] > 0) sig[v] |= 1ULL << smp;
        continue;
      }
      line.clear();
      for (uint32_t v = 1; v <= nv; v++) {
        if (!full && !W.in_proj[v]) continue;
        if (val2[v] < 0) line += '-';
        line += std::to_string(v);
        line += ' ';
      }
      line += "0\n";
      out << line;
    }
    if (atomic) {
      maybe_cache();
      const auto cnt = lit_counts(r, val, dead, vs, WD, W.in_proj, w);
      // Candidates: same literal count, same values in every sample. With
      // --cross each var is taken by its literal true in the first sample.
      std::map<std::pair<std::string, uint64_t>, std::vector<int>> groups;
      for (uint32_t v = 1; v <= nv; v++) {
        if (!W.in_proj[v]) continue;
        const int l = cross && !(sig[v] & 1) ? -(int)v : (int)v;
        const uint64_t sg = l > 0 ? sig[v] : ~sig[v];
        groups[{mpq_class(cnt[2 * v + (l < 0 ? 1 : 0)]).get_str(), sg}].push_back(l);
      }
      std::vector<std::vector<int>> cand, sets;
      for (auto& [k, g] : groups) if (g.size() >= 2) cand.push_back(std::move(g));
      // Confirm each member y against the first, r: {r, -y} and {-r, y} UNSAT.
      while (!cand.empty()) {
        std::vector<std::vector<int>> qs, next;
        for (const auto& g : cand)
          for (size_t i = 1; i < g.size(); i++) {
            qs.push_back({g[0], -g[i]});
            qs.push_back({-g[0], g[i]});
          }
        std::vector<std::string> ans;
        batch_sat(r, dead, W.in_proj, qs, assume, ans);
        size_t qi = 0;
        for (const auto& g : cand) {
          std::vector<int> eq{g[0]}, rest;
          for (size_t i = 1; i < g.size(); i++, qi += 2)
            (ans[qi] == "UNSAT" && ans[qi + 1] == "UNSAT" ? eq : rest).push_back(g[i]);
          if (eq.size() >= 2) sets.push_back(std::move(eq));
          if (rest.size() >= 2) next.push_back(std::move(rest));
        }
        cand = std::move(next);
      }
      for (auto& st : sets) {
        if (st[0] < 0) for (int& l : st) l = -l;
        std::sort(st.begin(), st.end(), [](int a, int b) { return std::abs(a) < std::abs(b); });
      }
      std::sort(sets.begin(), sets.end(), [](const auto& a, const auto& b) { return std::abs(a[0]) < std::abs(b[0]); });
      for (const auto& st : sets) {
        line.clear();
        for (int l : st) line += std::to_string(l) + " ";
        out << line << "0\n";
      }
    }
    out.flush();
    if (!out) die("write error");
    if (!cls.empty() && !atomic)
      std::cerr << "bdnnf " << cmd << ": " << num << " models satisfy the CNF" << std::endl;
  } catch (const bdnnf::Error& e) {
    die(path + ": " + e.what());
  }
  return EXIT_SUCCESS;
}

}  // namespace

int bdnnf_query(int argc, char** argv) {
  try {
    return run_query(argc, argv, nullptr);
  } catch (const QueryError& e) {
    std::cerr << "bdnnf: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}

// One query per stdin line, as on the command line without the file (and
// with the stream's --cnf unless the line gives one), on the circuit held in
// memory. Each answer ends with a "c done" line; errors give "c error <msg>".
int bdnnf_stream(int argc, char** argv) {
  std::string path, cnf;
  double cache_gb = 5;
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    if (a == "--cnf" && i + 1 < argc) cnf = argv[++i];
    else if (a == "--cache" && i + 1 < argc) cache_gb = std::stod(argv[++i]);
    else if (path.empty() && !a.empty() && a[0] != '-') path = a;
    else {
      std::cerr << "usage: bdnnf stream <in.bdnnf> [--cnf <file.cnf>] [--cache GB]" << std::endl;
      return EXIT_FAILURE;
    }
  }
  if (path.empty()) {
    std::cerr << "usage: bdnnf stream <in.bdnnf> [--cnf <file.cnf>] [--cache GB]" << std::endl;
    return EXIT_FAILURE;
  }
  bdnnf::Reader r;
  try {
    r.open(path);
    if (cache_gb > 0 && (double)std::filesystem::file_size(path) <= cache_gb * (1ULL << 30)) r.cache_all();
  } catch (const bdnnf::Error& e) {
    std::cerr << "bdnnf: " << path << ": " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  static const char* const cmds[] = {"count", "sample", "enum", "lits", "core", "batch", "atomic", "opt", "twise", "twise-check"};
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream ss(line);
    std::vector<std::string> tok;
    for (std::string t; ss >> t;) tok.push_back(t);
    if (tok.empty()) continue;
    if (tok[0] == "exit" || tok[0] == "quit") break;
    bool known = false;
    for (const char* c : cmds) known |= tok[0] == c;
    if (!known) {
      std::cout << "c error unknown query: " << tok[0] << "\nc done" << std::endl;
      continue;
    }
    std::vector<std::string> args{tok[0], path};
    bool has_cnf = false;
    for (size_t i = 1; i < tok.size(); i++) {
      has_cnf |= tok[i] == "--cnf";
      args.push_back(tok[i]);
    }
    if (!cnf.empty() && !has_cnf) { args.push_back("--cnf"); args.push_back(cnf); }
    std::vector<char*> av;
    for (auto& a : args) av.push_back(a.data());
    av.push_back(nullptr);
    try {
      run_query((int)args.size(), av.data(), &r);
    } catch (const QueryError& e) {
      std::cout << "c error " << e.what() << "\n";
    }
    std::cout << "c done" << std::endl;
  }
  return EXIT_SUCCESS;
}
