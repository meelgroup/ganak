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

// Fuzzer of our BuDDy fork, driven the way ganak's --buddy drives it: random
// CNFs conjoined clause by clause (bdd_clause, AND, bdd_appex early
// quantification) at local or global levels, node caps that force give-ups,
// tiny operator caches, var-count growth, init/done cycles. Small CNFs are
// checked against truth tables (function, canonicity, all counts, weighted
// mod p), wide ones against each other; BDDs kept across rounds are
// rechecked later, and after each round no node may leak.
// Usage: buddy-fuzz [--num N] [--seed S] [-v]   (no --num: forever)

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "bdd.h"

namespace {

constexpr uint64_t MOD = (1ULL << 61) - 1;
__extension__ typedef unsigned __int128 u128;
uint64_t mulm(uint64_t a, uint64_t b) {
  const u128 x = (u128)a * b;
  uint64_t r = (uint64_t)(x & MOD) + (uint64_t)(x >> 61);
  return r >= MOD ? r - MOD : r;
}
uint64_t addm(uint64_t a, uint64_t b) { return (a + b) % MOD; }

bool out_of_nodes = false;
uint64_t num_gc = 0;
void err_handler(int e) {
  if (e == BDD_NODENUM) { out_of_nodes = true; return; }
  fprintf(stderr, "buddy-fuzz: BuDDy error: %s\n", bdd_errstring(e));
  abort();
}
void gbc_handler(int pre, bddGbcStat*) { if (!pre) num_gc++; }

std::mt19937_64 rng;
uint64_t rnd(uint64_t n) { return std::uniform_int_distribution<uint64_t>(0, n-1)(rng); }
bool coin(uint32_t pct) { return rnd(100) < pct; }

struct Lit { uint32_t v; bool pos; };
struct Inst {
  uint32_t n = 0;
  uint32_t ni = 0;        // local vars < ni are indep, the rest projected away
  bool global = false;
  std::vector<int> level; // per local var; indep vars above the others
  std::vector<std::vector<Lit>> cls;
  std::vector<std::array<uint64_t, 2>> w;
};

struct Kept {
  bdd f;
  Inst in;
  std::vector<uint8_t> e;
  uint64_t cnt;
};

struct {
  uint64_t rounds = 0, wide = 0, epochs = 0, tight = 0, gave_up = 0, unsat = 0;
  uint64_t quant = 0, rechecks = 0, grown = 0, shannon = 0, retried = 0, siblings = 0, work_capped = 0;
} st;

uint64_t round_no = 0;
uint64_t seed = 1;
int varnum = 0;
int cap = 0;
bool verb = false;
std::vector<Kept> pool;

[[noreturn]] void fail(const Inst& in, const std::string& what) {
  fprintf(stderr, "buddy-fuzz: FAIL in round %llu (--seed %llu): %s\n",
      (unsigned long long)round_no, (unsigned long long)seed, what.c_str());
  fprintf(stderr, "c n %u indep %u %s, varnum %d, cap %d\nc levels:", in.n, in.ni,
      in.global ? "global" : "local", varnum, cap);
  for (uint32_t v = 0; v < in.n; v++) fprintf(stderr, " %u:%d", v+1, in.level[v]);
  fprintf(stderr, "\np cnf %u %zu\n", in.n, in.cls.size());
  for (const auto& c : in.cls) {
    for (const auto& l : c) fprintf(stderr, "%d ", l.pos ? (int)l.v+1 : -(int)l.v-1);
    fprintf(stderr, "0\n");
  }
  exit(1);
}

void add_random_clause(Inst& in, const uint32_t maxlen, const uint32_t window) {
  // by level, so a window is a set of vars close in the order
  std::vector<uint32_t> by_level(in.n);
  std::iota(by_level.begin(), by_level.end(), 0);
  std::sort(by_level.begin(), by_level.end(), [&](uint32_t a, uint32_t b) { return in.level[a] < in.level[b]; });
  const uint32_t len = std::min<uint32_t>(window, coin(5) || maxlen < 2 ? 1 : 2 + rnd(maxlen - 1));
  const uint32_t start = rnd(in.n - window + 1);
  std::vector<uint32_t> vs;
  while (vs.size() < len) {
    const uint32_t v = by_level[start + rnd(window)];
    if (std::find(vs.begin(), vs.end(), v) == vs.end()) vs.push_back(v);
  }
  std::vector<Lit> c;
  for (const uint32_t v : vs) c.push_back({v, coin(50)});
  in.cls.push_back(c);
  if (coin(5)) in.cls.push_back(c);
}

void reweight(Inst& in) {
  in.w.clear();
  for (uint32_t v = 0; v < in.n; v++) {
    std::array<uint64_t, 2> w;
    for (auto& x : w) x = coin(10) ? rnd(3) : rnd(MOD);
    in.w.push_back(w);
  }
}

Inst gen(const bool wide) {
  Inst in;
  in.n = wide ? 15 + rnd(34) : 1 + rnd(14);
  in.ni = coin(50) ? in.n : rnd(in.n+1);
  in.global = coin(40);
  std::vector<int> lv(in.n);
  if (in.global) {
    std::unordered_set<int> used;
    while (used.size() < in.n) used.insert((int)rnd(varnum));
    lv.assign(used.begin(), used.end());
    std::sort(lv.begin(), lv.end());
  } else std::iota(lv.begin(), lv.end(), 0);
  std::shuffle(lv.begin(), lv.begin()+in.ni, rng);
  std::shuffle(lv.begin()+in.ni, lv.end(), rng);
  in.level = lv;

  const uint32_t maxlen = std::min<uint32_t>(in.n, 2 + rnd(wide ? 3 : 4));
  const bool gates = in.n >= 3 && coin(30);
  uint32_t m = wide ? in.n/2 + rnd(2*in.n) : rnd(3*in.n + 1);
  if (gates) m /= 4;
  const uint32_t window = coin(50) ? in.n : std::min<uint32_t>(in.n, 3 + rnd(6));
  for (uint32_t i = 0; i < m; i++) add_random_clause(in, maxlen, window);
  // Tseitin definitions, as in count_fuzzer's circuit-like CNFs
  if (gates) for (uint32_t g = 0, num = 1 + rnd(in.n/2 + 1); g < num; g++) {
    uint32_t o = rnd(in.n), a = rnd(in.n), b = rnd(in.n);
    while (a == o) a = rnd(in.n);
    while (b == o || b == a) b = rnd(in.n);
    const bool sa = coin(50), sb = coin(50);
    switch (rnd(3)) {
      case 0: // o = a & b
        in.cls.push_back({{o, false}, {a, sa}});
        in.cls.push_back({{o, false}, {b, sb}});
        in.cls.push_back({{o, true}, {a, !sa}, {b, !sb}});
        break;
      case 1: // o = a | b
        in.cls.push_back({{o, true}, {a, !sa}});
        in.cls.push_back({{o, true}, {b, !sb}});
        in.cls.push_back({{o, false}, {a, sa}, {b, sb}});
        break;
      default: // o = a ^ b
        in.cls.push_back({{o, false}, {a, true}, {b, true}});
        in.cls.push_back({{o, false}, {a, false}, {b, false}});
        in.cls.push_back({{o, true}, {a, false}, {b, true}});
        in.cls.push_back({{o, true}, {a, true}, {b, false}});
    }
  }
  std::shuffle(in.cls.begin(), in.cls.end(), rng);
  reweight(in);
  return in;
}

// ganak's sibling comps share vars, levels and most clauses: their
// sub-BDDs and operator cache entries are reused, so check that reuse
Inst sibling(const Inst& prev) {
  Inst in = prev;
  in.cls.clear();
  for (const auto& c : prev.cls) if (!coin(25)) in.cls.push_back(c);
  for (uint32_t i = 0, num = rnd(4); i < num; i++) add_random_clause(in, std::min<uint32_t>(in.n, 3), in.n);
  if (coin(50)) reweight(in);
  return in;
}

std::vector<uint32_t> clause_order(const Inst& in) {
  std::vector<uint32_t> ord(in.cls.size());
  std::iota(ord.begin(), ord.end(), 0);
  const uint32_t how = rnd(3);
  if (how == 0) {
    // ganak's: clauses whose topmost var is deepest first
    std::vector<int> top(in.cls.size(), INT32_MAX);
    for (uint32_t c = 0; c < in.cls.size(); c++)
      for (const auto& l : in.cls[c]) top[c] = std::min(top[c], in.level[l.v]);
    std::sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b) {
      if (top[a] != top[b]) return top[a] > top[b];
      return a < b;
    });
  } else if (how == 1) std::shuffle(ord.begin(), ord.end(), rng);
  return ord;
}

// As ganak's buddy_count. False if it gave up.
bool build(const Inst& in, const std::vector<uint32_t>& ord, const bool quant, bdd& out) {
  bool retry = coin(30);
  std::vector<int> last;
  if (quant) {
    last.assign(in.n, -1);
    for (uint32_t p = 0; p < ord.size(); p++)
      for (const auto& l : in.cls[ord[p]]) last[l.v] = (int)p;
  }
  out_of_nodes = false;
  // ganak's --buddymaxwork: give up after so many op cache misses
  const bool work_cap = coin(30);
  if (work_cap) {
    bddCacheStat cs;
    bdd_cachestats(&cs);
    bdd_setopmisslimit(cs.opMiss + 1 + rnd(coin(50) ? 50 : 3000));
    st.work_capped++;
  }
  BDD f = 1;
  std::vector<std::pair<int, int>> lv;
  std::vector<int> levels, pos, qlevels;
  for (uint32_t p = 0; p < ord.size(); p++) {
    lv.clear();
    qlevels.clear();
    for (const auto& l : in.cls[ord[p]]) {
      lv.emplace_back(in.level[l.v], l.pos);
      if (quant && l.v >= in.ni && last[l.v] == (int)p) {
        qlevels.push_back(in.level[l.v]);
        last[l.v] = -1;
      }
    }
    std::sort(lv.begin(), lv.end());
    levels.clear();
    pos.clear();
    for (const auto& [l, s] : lv) { levels.push_back(l); pos.push_back(s); }
    const BDD cl = bdd_addref(bdd_clause(levels.data(), pos.data(), (int)lv.size()));
    const bool clause_ok = !out_of_nodes;
    auto op = [&]() {
      if (qlevels.empty()) return bdd_apply(f, cl, bddop_and);
      const bdd vs = bdd_makeset(qlevels.data(), (int)qlevels.size());
      return bdd_appex(f, cl, bddop_and, vs.id());
    };
    BDD next = op();
    // Gave up: the same op again, on the same live operands, with room. The
    // operator caches must not keep what was computed without nodes.
    if (out_of_nodes && clause_ok && retry) {
      bdd_clear_error();
      st.gave_up++;
      st.retried++;
      retry = false;
      if (work_cap) bdd_setopmisslimit(0);
      else {
        cap *= 2;
        bdd_setmaxnodenum(cap);
      }
      out_of_nodes = false;
      next = op();
    }
    bdd_addref(next);
    bdd_delref(cl);
    bdd_delref(f);
    f = next;
    if (out_of_nodes) break;
  }
  bdd_setopmisslimit(0);
  if (out_of_nodes) {
    bdd_delref(f);
    bdd_clear_error();
    st.gave_up++;
    return false;
  }
  out = bdd_fromroot(f);
  bdd_delref(f);
  return true;
}

bdd varset(const Inst& in, uint32_t from, uint32_t to) {
  std::vector<int> ls;
  for (uint32_t v = from; v < to; v++) ls.push_back(in.level[v]);
  return bdd_makeset(ls.data(), (int)ls.size());
}

// Reduced, ordered, and its support within [from, to) of the local vars
void check_shape(const Inst& in, const bdd& f, uint32_t from, uint32_t to, const char* name) {
  std::vector<int> loc(varnum+1, -1);
  for (uint32_t v = 0; v < in.n; v++) loc[in.level[v]] = (int)v;
  std::unordered_set<int> seen;
  std::vector<int> todo{f.id()};
  while (!todo.empty()) {
    const int x = todo.back();
    todo.pop_back();
    if (x < 2 || !seen.insert(x).second) continue;
    const int l = bdd_var(x);
    if (l < 0 || l >= varnum || loc[l] < 0 || (uint32_t)loc[l] < from || (uint32_t)loc[l] >= to)
      fail(in, std::string(name) + ": node at level " + std::to_string(l) + " outside its support");
    const int lo = bdd_low(x), hi = bdd_high(x);
    if (lo == hi) fail(in, std::string(name) + ": node with equal children");
    for (const int c : {lo, hi})
      if (c >= 2 && bdd_var(c) <= l) fail(in, std::string(name) + ": child not below its parent");
    todo.push_back(lo);
    todo.push_back(hi);
  }
}

bool eval(const bdd& f, const std::vector<int>& loc, uint64_t a) {
  int x = f.id();
  while (x >= 2) x = ((a >> loc[bdd_var(x)]) & 1) ? bdd_high(x) : bdd_low(x);
  return x == 1;
}

std::vector<int> local_of_level(const Inst& in) {
  std::vector<int> loc(varnum+1, -1);
  for (uint32_t v = 0; v < in.n; v++) loc[in.level[v]] = (int)v;
  return loc;
}

bool sat_by(const Inst& in, uint64_t a) {
  for (const auto& c : in.cls) {
    bool ok = false;
    for (const auto& l : c) if ((((a >> l.v) & 1) != 0) == l.pos) { ok = true; break; }
    if (!ok) return false;
  }
  return true;
}

std::vector<int> ranks(const Inst& in, std::vector<uint32_t>& var_of_rank) {
  std::vector<int> rank(varnum+1, (int)in.n);
  var_of_rank.resize(in.n);
  std::iota(var_of_rank.begin(), var_of_rank.end(), 0);
  std::sort(var_of_rank.begin(), var_of_rank.end(), [&](uint32_t a, uint32_t b) { return in.level[a] < in.level[b]; });
  for (uint32_t r = 0; r < in.n; r++) rank[in.level[var_of_rank[r]]] = (int)r;
  return rank;
}

// ganak's field-valued walk (buddy_field_count), mod p
uint64_t walk(const Inst& in, const bdd& f, const bool unit) {
  std::vector<uint32_t> vr;
  const std::vector<int> rank = ranks(in, vr);
  const uint32_t ni = in.ni;
  std::vector<uint64_t> sw(in.n);
  for (uint32_t r = 0; r < in.n; r++) sw[r] = unit ? 2 : addm(in.w[vr[r]][0], in.w[vr[r]][1]);
  auto skip = [&](uint32_t a, uint32_t b) {
    uint64_t p = 1;
    for (uint32_t r = a; r < b; r++) p = mulm(p, sw[r]);
    return p;
  };
  auto frank = [&](int x) { return x < 2 ? (int)in.n : rank[bdd_var(x)]; };
  std::unordered_map<int, uint64_t> memo;
  auto rec = [&](auto&& self, int x) -> uint64_t {
    if (x < 2) return (uint64_t)x;
    const int r = rank[bdd_var(x)];
    if (r >= (int)ni) return 1;
    if (auto it = memo.find(x); it != memo.end()) return it->second;
    uint64_t s = 0;
    for (const int sign : {0, 1}) {
      const int c = sign ? bdd_high(x) : bdd_low(x);
      if (c == 0) continue;
      uint64_t t = mulm(self(self, c), skip(r+1, std::min<uint32_t>(frank(c), ni)));
      if (!unit) t = mulm(t, in.w[vr[r]][sign]);
      s = addm(s, t);
    }
    memo[x] = s;
    return s;
  };
  const int root = f.id();
  if (root == 0) return 0;
  return mulm(skip(0, std::min<uint32_t>(frank(root), ni)), rec(rec, root));
}

uint64_t count_local(const Inst& in, const bdd& f) {
  const uint32_t proj_end = (in.ni == in.n) ? 63 : in.ni;
  const uint64_t c = bdd_satcount_i64(f, proj_end);
  return (proj_end == 63) ? c >> (63 - in.n) : c;
}

uint64_t count_rank(const Inst& in, const bdd& f) {
  std::vector<uint32_t> vr;
  std::vector<int> rank = ranks(in, vr);
  return bdd_satcount_i64_rank(f, rank.data(), (int)in.ni);
}

bool is_local(const Inst& in) {
  for (uint32_t v = 0; v < in.n; v++) if (in.level[v] >= (int)in.n) return false;
  return true;
}

void check_counts(const Inst& in, const bdd& f, const uint64_t cnt, const char* name) {
  const std::string nm(name);
  if (is_local(in)) {
    const uint64_t c = count_local(in, f);
    if (c != cnt) fail(in, nm + ": bdd_satcount_i64 " + std::to_string(c) + " != " + std::to_string(cnt));
  }
  const uint64_t r = count_rank(in, f);
  if (r != cnt) fail(in, nm + ": bdd_satcount_i64_rank " + std::to_string(r) + " != " + std::to_string(cnt));
  const uint64_t u = walk(in, f, true);
  if (u != cnt % MOD) fail(in, nm + ": unit-weight walk " + std::to_string(u) + " != " + std::to_string(cnt));
}

bdd shannon(const Inst& in, const std::vector<uint8_t>& e, const std::vector<uint32_t>& by_level,
    uint32_t depth, uint64_t a) {
  if (depth == in.n) return e[a] ? bddtrue : bddfalse;
  const uint32_t v = by_level[depth];
  const bdd lo = shannon(in, e, by_level, depth+1, a);
  const bdd hi = shannon(in, e, by_level, depth+1, a | (1ULL << v));
  return bdd_ite(bdd_ithvar(in.level[v]), hi, lo);
}

void recheck(const Kept& k) {
  check_shape(k.in, k.f, 0, k.in.n, "kept");
  const std::vector<int> loc = local_of_level(k.in);
  for (uint64_t a = 0; a < k.e.size(); a++)
    if (eval(k.f, loc, a) != (bool)k.e[a]) fail(k.in, "kept BDD changed its function");
  check_counts(k.in, k.f, k.cnt, "kept");
  st.rechecks++;
}

Inst prev_small, prev_wide;

void small_round() {
  const Inst in = (prev_small.n > 0 && coin(40)) ? (st.siblings++, sibling(prev_small)) : gen(false);
  prev_small = in;
  const uint32_t n = in.n, ni = in.ni;
  const uint64_t full = 1ULL << n;
  std::vector<uint8_t> t(full), p(1ULL << ni, 0), e(full);
  for (uint64_t a = 0; a < full; a++) {
    t[a] = sat_by(in, a);
    if (t[a]) p[a & ((1ULL << ni) - 1)] = 1;
  }
  uint64_t cnt = 0, wcnt = 0;
  for (uint64_t a = 0; a < p.size(); a++) {
    if (!p[a]) continue;
    cnt++;
    uint64_t w = 1;
    for (uint32_t v = 0; v < ni; v++) w = mulm(w, in.w[v][(a >> v) & 1]);
    wcnt = addm(wcnt, w);
  }

  const bool quant = ni < n && coin(70);
  if (quant) st.quant++;
  for (uint64_t a = 0; a < full; a++) e[a] = quant ? p[a & ((1ULL << ni) - 1)] : t[a];
  bdd f;
  if (!build(in, clause_order(in), quant, f)) return;
  if (f == bddfalse) st.unsat++;
  check_shape(in, f, 0, quant ? ni : n, "built");
  const std::vector<int> loc = local_of_level(in);
  for (uint64_t a = 0; a < full; a++)
    if (eval(f, loc, a) != (bool)e[a]) fail(in, "wrong function at assignment " + std::to_string(a));
  check_counts(in, f, cnt, "built");
  const uint64_t wc = walk(in, f, false);
  if (wc != wcnt) fail(in, "weighted walk " + std::to_string(wc) + " != " + std::to_string(wcnt));

  // the other way round: quantified if f is not, and the other one through bdd_exist
  bdd g;
  const bool gquant = !quant && ni < n;
  if (build(in, clause_order(in), gquant, g)) {
    if (ni < n) {
      out_of_nodes = false;
      const bdd vs = varset(in, ni, n);
      const bdd q = gquant ? g : f;
      const bdd ex = bdd_exist(gquant ? f : g, vs);
      if (out_of_nodes) { bdd_clear_error(); st.gave_up++; }
      else if (ex != q) fail(in, "early quantification gives another node than bdd_exist");
    } else if (f != g) fail(in, "another clause order: another node");
  }
  if (n <= 12) {
    std::vector<uint32_t> by_level(n);
    std::iota(by_level.begin(), by_level.end(), 0);
    std::sort(by_level.begin(), by_level.end(), [&](uint32_t a, uint32_t b) { return in.level[a] < in.level[b]; });
    out_of_nodes = false;
    const bdd s = shannon(in, e, by_level, 0, 0);
    if (out_of_nodes) { bdd_clear_error(); st.gave_up++; }
    else {
      st.shannon++;
      if (s != f) fail(in, "Shannon expansion of the truth table: another node");
    }
  }
  if (coin(20)) {
    if (pool.size() >= 6) pool.erase(pool.begin() + rnd(pool.size()));
    pool.push_back(Kept{f, in, e, cnt});
  }
}

void wide_round() {
  const Inst in = (prev_wide.n > 0 && coin(40)) ? (st.siblings++, sibling(prev_wide)) : gen(true);
  prev_wide = in;
  const uint32_t n = in.n, ni = in.ni;
  st.wide++;
  const bool quant = ni < n && coin(70);
  if (quant) st.quant++;
  bdd a, b;
  if (!build(in, clause_order(in), quant, a)) return;
  if (!build(in, clause_order(in), false, b)) return;
  check_shape(in, a, 0, quant ? ni : n, "wide");
  check_shape(in, b, 0, n, "wide unquantified");
  const std::vector<int> loc = local_of_level(in);
  for (int i = 0; i < 64; i++) {
    const uint64_t x = rng() & ((n == 64) ? ~0ULL : ((1ULL << n) - 1));
    if (eval(b, loc, x) != sat_by(in, x)) fail(in, "wide: wrong function at a sampled assignment");
  }
  if (a == bddfalse) st.unsat++;
  bdd ref = b;
  if (ni < n) {
    out_of_nodes = false;
    ref = bdd_exist(b, varset(in, ni, n));
    if (out_of_nodes) { bdd_clear_error(); st.gave_up++; return; }
  }
  if (quant && ref != a) fail(in, "wide: early quantification gives another node than bdd_exist");
  if (!quant && b != a) fail(in, "wide: another clause order gives another node");
  const uint64_t c = count_rank(in, ref);
  check_counts(in, a, c, "wide");
  if (walk(in, a, false) != walk(in, ref, false)) fail(in, "wide: weighted walks differ");
}

void check_leaks() {
  bdd_gbc();
  std::unordered_set<int> live;
  std::vector<int> todo;
  for (int l = 0; l < varnum; l++) {
    todo.push_back(bdd_ithvar(l).id());
    todo.push_back(bdd_nithvar(l).id());
  }
  for (const auto& k : pool) todo.push_back(k.f.id());
  while (!todo.empty()) {
    const int x = todo.back();
    todo.pop_back();
    if (x < 2 || !live.insert(x).second) continue;
    todo.push_back(bdd_low(x));
    todo.push_back(bdd_high(x));
  }
  const int used = bdd_getnodenum();
  if (used != (int)live.size() + 2) {
    fprintf(stderr, "buddy-fuzz: FAIL after round %llu (--seed %llu): %d nodes in use, %zu reachable\n",
        (unsigned long long)round_no, (unsigned long long)seed, used, live.size() + 2);
    exit(1);
  }
}

void start_epoch() {
  const bool tight = coin(30);
  cap = tight ? 300 + (int)rnd(1500) : 20000 + (int)rnd(200000);
  static const int caches[] = {1, 3, 17, 100, 2000, 10000};
  bdd_init(std::min(10000, cap/2), caches[rnd(6)]);
  bdd_error_hook(err_handler);
  bdd_gbc_hook(gbc_handler);
  varnum = 63;
  bdd_setvarnum(varnum);
  bdd_autoreorder(BDD_REORDER_NONE);
  bdd_setmaxnodenum(cap);
  st.epochs++;
  if (tight) st.tight++;
}

void end_epoch() {
  pool.clear();
  // its levels may be above the next epoch's varnum
  prev_small = prev_wide = Inst();
  check_leaks();
  bdd_done();
}

void grow_varnum() {
  const int v = varnum + 1 + (int)rnd(150);
  // as ganak: the new var nodes are permanent, the budget is for the BDDs
  cap += 2*(v+1) + 2;
  bdd_setmaxnodenum(cap);
  out_of_nodes = false;
  if (bdd_setvarnum(v) < 0 || out_of_nodes) {
    fprintf(stderr, "buddy-fuzz: FAIL in round %llu (--seed %llu): bdd_setvarnum(%d) failed\n",
        (unsigned long long)round_no, (unsigned long long)seed, v);
    exit(1);
  }
  varnum = v;
  st.grown++;
}

}

int main(int argc, char** argv) {
  uint64_t num = 0;
  bool forever = true;
  seed = std::random_device{}();
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    if (a == "--num" && i+1 < argc) { num = std::stoull(argv[++i]); forever = false; }
    else if (a == "--seed" && i+1 < argc) seed = std::stoull(argv[++i]);
    else if (a == "-v") verb = true;
    else {
      fprintf(stderr, "Usage: buddy-fuzz [--num N] [--seed S] [-v]\n");
      return 2;
    }
  }
  rng.seed(seed);
  printf("buddy-fuzz: seed %llu\n", (unsigned long long)seed);

  uint64_t epoch_left = 0;
  for (round_no = 0; forever || round_no < num; round_no++) {
    if (epoch_left == 0) {
      if (round_no > 0) end_epoch();
      start_epoch();
      epoch_left = 1 + rnd(200);
    }
    epoch_left--;
    if (coin(5)) grow_varnum();
    if (!pool.empty() && coin(30)) recheck(pool[rnd(pool.size())]);
    if (coin(25)) wide_round();
    else small_round();
    if (!pool.empty() && coin(10)) pool.erase(pool.begin() + rnd(pool.size()));
    check_leaks();
    st.rounds++;
    if ((verb && st.rounds % 1000 == 0) || (forever && st.rounds % 100000 == 0))
      printf("buddy-fuzz: %llu rounds\n", (unsigned long long)st.rounds);
  }
  end_epoch();
  printf("buddy-fuzz: OK, %llu rounds (%llu wide, %llu quantified, %llu unsat), %llu Shannon checks, "
      "%llu kept rechecks, %llu siblings; %llu epochs (%llu with a tight cap), %llu give-ups (%llu retried), "
      "%llu work-capped builds, %llu GCs, %llu var growths\n",
      (unsigned long long)st.rounds, (unsigned long long)st.wide, (unsigned long long)st.quant,
      (unsigned long long)st.unsat, (unsigned long long)st.shannon, (unsigned long long)st.rechecks,
      (unsigned long long)st.siblings, (unsigned long long)st.epochs, (unsigned long long)st.tight,
      (unsigned long long)st.gave_up, (unsigned long long)st.retried,
      (unsigned long long)st.work_capped, (unsigned long long)num_gc, (unsigned long long)st.grown);
  return 0;
}
