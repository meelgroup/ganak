/******************************************
Copyright (C) 2026 Authors of GANAK, see AUTHORS file

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

// Counts the models of a WHOLE input CNF that encodes a labeled exact cover,
// such as a Sudoku, Latin square or quasigroup completion, by branching on
// complete placements of a label instead of on single variables
// (by Christian Reitwiessner, PR #133).
//
// Accepted shape: every clause is a unit, a negative binary or an all-positive
// clause, and the CNF proves every positive clause exactly-one (all of its
// pairs are binaries, or a counting argument over disjoint binary cliques).
// One family of these constraints partitions the variables into cells; the
// remaining constraints connect the variables into labels with at most one
// variable per cell. A model then picks, for every label, a set of cells that
// meets each of that label's constraints exactly once (a placement), and the
// placements of all labels partition the cells. Anything else declines.
#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <numeric>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <gmpxx.h>

namespace GanakInt::placement {

constexpr uint32_t max_vars = 8192;
constexpr uint32_t max_cells = 256;
constexpr uint32_t max_labels = 64;
constexpr uint32_t max_scopes = 64; // constraints per label, one bit each

struct Structure {
  uint32_t cells = 0;
  // Per label: its constraints as lists of cells.
  std::vector<std::vector<std::vector<uint32_t>>> scopes;
  // Per label and cell: the label may take the cell / a unit forces it there.
  std::vector<std::vector<uint8_t>> allowed, required;
  bool contradictory = false;
  const char* certificate = "";
};

struct Stats {
  std::vector<uint64_t> domains;
  uint64_t enum_nodes = 0, join_nodes = 0, memo_hits = 0;
  double estimate = 0; // estimated join work, in placement comparisons
  bool fast_join = false;
};

namespace detail {
struct Union {
  std::vector<uint32_t> p;
  explicit Union(uint32_t n) : p(n) { std::iota(p.begin(), p.end(), 0u); }
  uint32_t find(uint32_t x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
  void join(uint32_t x, uint32_t y) { p[find(x)] = find(y); }
};

using Scope = std::vector<uint32_t>;

// The family is a set of exactly-one constraints over [0,n); `primary` is a
// candidate partition of the variables into cells.
inline std::optional<Structure> build(uint32_t n, const std::vector<Scope>& family,
    const std::vector<Scope>& primary, const std::vector<std::pair<uint32_t,bool>>& units) {
  if (primary.size() < 2 || primary.size() > max_cells) return {};
  std::vector<int32_t> owner(n, -1);
  for (uint32_t c = 0; c < primary.size(); c++) for (auto v : primary[c]) {
    if (owner[v] != -1) return {};
    owner[v] = c;
  }
  if (std::count(owner.begin(), owner.end(), -1)) return {};
  std::vector<Scope> rest;
  {
    std::vector<Scope> sorted_primary = primary;
    std::sort(sorted_primary.begin(), sorted_primary.end());
    for (const auto& g : family)
      if (!std::binary_search(sorted_primary.begin(), sorted_primary.end(), g)) rest.push_back(g);
  }
  Union labels(n);
  std::vector<uint8_t> constrained(n);
  for (const auto& g : rest) for (auto v : g) { labels.join(v, g.front()); constrained[v] = 1; }
  // A variable in no other constraint would let its label take any cell set
  if (std::count(constrained.begin(), constrained.end(), 0)) return {};
  std::vector<int32_t> id(n, -1), label(n);
  uint32_t k = 0;
  for (uint32_t v = 0; v < n; v++) {
    auto r = labels.find(v);
    if (id[r] == -1) id[r] = k++;
    label[v] = id[r];
  }
  if (k < 2 || k > max_labels) return {};
  Structure s;
  s.cells = primary.size();
  s.scopes.resize(k);
  s.allowed.assign(k, std::vector<uint8_t>(s.cells));
  s.required.assign(k, std::vector<uint8_t>(s.cells));
  for (uint32_t v = 0; v < n; v++) {
    auto& a = s.allowed[label[v]][owner[v]];
    if (a) return {}; // two variables of one label in one cell
    a = 1;
  }
  for (const auto& g : rest) {
    Scope cells;
    for (auto v : g) cells.push_back(owner[v]);
    std::sort(cells.begin(), cells.end());
    auto& sc = s.scopes[label[g.front()]];
    if (sc.size() >= max_scopes) return {};
    sc.push_back(std::move(cells));
  }
  for (auto [v, positive] : units) {
    const uint32_t l = label[v], c = owner[v];
    if (!positive) { s.allowed[l][c] = 0; continue; }
    s.required[l][c] = 1;
    for (uint32_t o = 0; o < k; o++) if (o != l) s.allowed[o][c] = 0;
  }
  for (uint32_t l = 0; l < k; l++) for (uint32_t c = 0; c < s.cells; c++)
    if (s.required[l][c] && !s.allowed[l][c]) s.contradictory = true;
  return s;
}

// Which constraints are the cells? In a Latin square rows, columns and cells
// play symmetric roles, and in a Sudoku rows also partition the variables, so
// search for a partition of the variables into constraints (an exact cover).
// Choosing a cell rules out every constraint meeting it, so those belong to
// the labels: they are merged into label components at once, and the search
// backtracks as soon as a component holds two variables of one chosen cell.
class CellSearch {
  uint32_t n;
  const std::vector<Scope>& family;
  std::vector<std::vector<uint32_t>> of_var;   // constraints containing a variable
  std::vector<int32_t> block;                  // chosen constraint covering a var, or -1
  std::vector<uint32_t> killed_by;             // 0: constraint still available
  std::vector<uint32_t> chosen;
  // Union-find with undo over variables; per root the chosen cells it meets
  std::vector<uint32_t> parent, sz;
  std::vector<std::vector<uint32_t>> cells_of;
  struct Undo { uint32_t child, root, old_cells; };
  std::vector<Undo> trail;
  std::vector<uint32_t> stamp;
  uint32_t now = 0;
  uint64_t nodes = 0;
  static constexpr uint64_t max_nodes = 200000;

  uint32_t find(uint32_t x) const { while (parent[x] != x) x = parent[x]; return x; }
  bool unite(uint32_t a, uint32_t b) {
    a = find(a); b = find(b);
    if (a == b) return true;
    if (sz[a] > sz[b]) std::swap(a, b);
    // Two variables of one cell may never share a label
    ++now;
    for (auto c : cells_of[b]) stamp[c] = now;
    bool ok = true;
    for (auto c : cells_of[a]) if (stamp[c] == now) ok = false;
    trail.push_back({a, b, (uint32_t)cells_of[b].size()});
    parent[a] = b; sz[b] += sz[a];
    cells_of[b].insert(cells_of[b].end(), cells_of[a].begin(), cells_of[a].end());
    return ok;
  }
  bool add_cell(uint32_t v, uint32_t c) {
    const uint32_t r = find(v);
    for (auto x : cells_of[r]) if (x == c) return false;
    trail.push_back({r, r, (uint32_t)cells_of[r].size()});
    cells_of[r].push_back(c);
    return true;
  }
  void undo_to(size_t mark) {
    while (trail.size() > mark) {
      auto u = trail.back(); trail.pop_back();
      if (u.child != u.root) { parent[u.child] = u.child; sz[u.root] -= sz[u.child]; }
      cells_of[u.root].resize(u.old_cells);
    }
  }
  // Choose constraint g as a cell; returns false on a label conflict
  bool choose(uint32_t g, uint32_t depth) {
    bool ok = true;
    for (auto v : family[g]) {
      block[v] = g;
      if (!add_cell(v, g)) ok = false;
    }
    for (auto v : family[g]) for (auto h : of_var[v]) {
      if (h == g || killed_by[h]) continue;
      killed_by[h] = depth;
      for (auto x : family[h]) if (!unite(x, family[h].front())) ok = false;
    }
    return ok;
  }
  void unchoose(uint32_t g, uint32_t depth) {
    for (auto v : family[g]) {
      block[v] = -1;
      for (auto h : of_var[v]) if (killed_by[h] == depth) killed_by[h] = 0;
    }
  }
  bool go(uint32_t depth) {
    if (++nodes > max_nodes) return false;
    // Most constrained uncovered variable
    int32_t best = -1;
    uint32_t best_live = ~0u;
    for (uint32_t v = 0; v < n; v++) if (block[v] == -1) {
      uint32_t live = 0;
      for (auto h : of_var[v]) live += !killed_by[h];
      if (live < best_live) { best_live = live; best = v; if (live <= 1) break; }
    }
    if (best == -1) return true;
    if (!best_live) return false;
    for (auto g : of_var[best]) {
      if (killed_by[g]) continue;
      killed_by[g] = depth; // the chosen constraint itself
      const size_t mark = trail.size();
      chosen.push_back(g);
      if (choose(g, depth) && go(depth+1)) return true;
      chosen.pop_back();
      unchoose(g, depth);
      killed_by[g] = 0;
      undo_to(mark);
      if (nodes > max_nodes) return false;
    }
    return false;
  }

public:
  CellSearch(uint32_t _n, const std::vector<Scope>& _family) : n(_n), family(_family),
    of_var(_n), block(_n, -1), killed_by(_family.size()), parent(_n), sz(_n, 1),
    cells_of(_n), stamp(_family.size()) {
    std::iota(parent.begin(), parent.end(), 0u);
    for (uint32_t g = 0; g < family.size(); g++) for (auto v : family[g]) of_var[v].push_back(g);
  }
  std::optional<std::vector<Scope>> run() {
    for (uint32_t v = 0; v < n; v++) if (of_var[v].empty()) return {};
    if (!go(1)) return {};
    std::vector<Scope> cells;
    for (auto g : chosen) cells.push_back(family[g]);
    return cells;
  }
};
} // namespace detail

template<class Clauses>
std::optional<Structure> decode(uint32_t n, const Clauses& clauses) {
  using detail::Scope;
  if (n < 4 || n > max_vars) return {};
  std::vector<Scope> pos;
  std::vector<uint64_t> bins;
  std::vector<std::pair<uint32_t,bool>> units;
  for (const auto& cl : clauses) {
    if (cl.empty()) return {};
    if (cl.size() == 1) { units.emplace_back(cl[0].var(), !cl[0].sign()); continue; }
    if (cl.size() == 2 && cl[0].sign() && cl[1].sign()) {
      uint32_t a = cl[0].var(), b = cl[1].var();
      if (a == b) return {};
      if (a > b) std::swap(a, b);
      bins.push_back(uint64_t(a)*n + b);
      continue;
    }
    Scope g;
    for (const auto& l : cl) { if (l.sign()) return {}; g.push_back(l.var()); }
    std::sort(g.begin(), g.end());
    if (std::adjacent_find(g.begin(), g.end()) != g.end()) return {};
    pos.push_back(std::move(g));
  }
  std::sort(pos.begin(), pos.end());
  pos.erase(std::unique(pos.begin(), pos.end()), pos.end());
  std::sort(bins.begin(), bins.end());
  bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
  if (pos.size() < 2) return {};

  // Pairs every positive clause needs for at-most-one. 64 vars per clause
  // keeps this list within a few million entries.
  std::vector<uint64_t> pairs;
  uint64_t num_pairs = 0;
  for (const auto& g : pos) {
    if (g.size() > 64) return {};
    num_pairs += g.size()*(g.size()-1)/2;
  }
  if (num_pairs > 8000000) return {};
  pairs.reserve(num_pairs);
  for (const auto& g : pos)
    for (size_t i = 0; i < g.size(); i++) for (size_t j = i+1; j < g.size(); j++)
      pairs.push_back(uint64_t(g[i])*n + g[j]);
  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());

  std::vector<Scope> family;
  std::vector<std::vector<Scope>> proposals; // tried before the search
  const char* certificate;
  if (pairs == bins) {
    certificate = "explicit";
    family = pos;
  } else {
    // Disjoint binary cliques (at most one per block) plus at-least-one
    // clauses in which every variable occurs equally often, d times. If m is
    // the number of positive clauses and m = d * #blocks then every model
    // sets exactly one variable per block, and every positive clause has
    // exactly one true variable: the d * #true incidences cover m clauses.
    certificate = "capacity";
    detail::Union u(n);
    for (auto e : bins) u.join(e / n, e % n);
    std::vector<uint64_t> size(n), edges(n);
    for (uint32_t v = 0; v < n; v++) size[u.find(v)]++;
    for (auto e : bins) edges[u.find(e / n)]++;
    std::vector<Scope> blocks;
    std::vector<int32_t> block_of(n, -1);
    for (uint32_t v = 0; v < n; v++) {
      const auto r = u.find(v);
      if (edges[r] != size[r]*(size[r]-1)/2) return {};
      if (block_of[r] == -1) { block_of[r] = blocks.size(); blocks.emplace_back(); }
      blocks[block_of[r]].push_back(v);
    }
    std::vector<uint64_t> degree(n);
    for (const auto& g : pos) for (auto v : g) degree[v]++;
    if (degree[0] == 0 || std::count(degree.begin(), degree.end(), degree[0]) != (int64_t)n) return {};
    if (pos.size() != degree[0] * blocks.size()) return {};
    family = pos;
    for (const auto& b : blocks) if (b.size() >= 2) family.push_back(b);
    std::sort(family.begin(), family.end());
    family.erase(std::unique(family.begin(), family.end()), family.end());
    proposals.push_back(std::move(blocks));
  }
  std::sort(family.begin(), family.end());
  for (const auto& p : proposals) {
    if (auto s = detail::build(n, family, p, units)) {
      s->certificate = certificate;
      return s;
    }
  }
  if (auto p = detail::CellSearch(n, family).run()) {
    if (auto s = detail::build(n, family, *p, units)) {
      s->certificate = certificate;
      return s;
    }
  }
  return {};
}

namespace detail {
template<uint32_t W>
struct Cells {
  std::array<uint64_t, W> w{};
  void add(uint32_t c) { w[c/64] |= uint64_t(1) << (c%64); }
  bool has(uint32_t c) const { return (w[c/64] >> (c%64)) & 1; }
  bool disjoint(const Cells& o) const {
    uint64_t x = 0;
    for (uint32_t i = 0; i < W; i++) x |= w[i] & o.w[i];
    return x == 0;
  }
  bool subset_of(const Cells& o) const {
    uint64_t x = 0;
    for (uint32_t i = 0; i < W; i++) x |= w[i] & ~o.w[i];
    return x == 0;
  }
  Cells operator|(const Cells& o) const { Cells r; for (uint32_t i = 0; i < W; i++) r.w[i] = w[i] | o.w[i]; return r; }
  Cells minus(const Cells& o) const { Cells r; for (uint32_t i = 0; i < W; i++) r.w[i] = w[i] & ~o.w[i]; return r; }
  bool operator==(const Cells& o) const { return w == o.w; }
};
template<uint32_t W>
struct CellsHash {
  size_t operator()(const Cells<W>& c) const {
    uint64_t h = 0x9e3779b97f4a7c15ULL;
    for (auto x : c.w) { h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); }
    return h;
  }
};

template<uint32_t W>
class Counter {
  using C = Cells<W>;
  const Structure& s;
  uint64_t budget;
  Stats& st;
  bool exhausted = false;
  C all;

  // Exact covers of one label's constraints.
  struct Enum {
    std::vector<C> scope_cells;
    std::vector<uint64_t> incidence; // per cell: bit i if the cell meets constraint i
    std::vector<uint32_t> required;
    std::vector<C>* out;
    uint64_t* budget;
    uint64_t* nodes;
    uint64_t max_out;
    bool fail = false;
    void go(uint64_t todo, C used) {
      if (fail) return;
      if (++*nodes, !(*budget)--) { fail = true; return; }
      for (auto c : required)
        if (!used.has(c) && (incidence[c] & ~todo)) return; // can no longer take c
      if (!todo) {
        if (out->size() >= max_out) { fail = true; return; }
        out->push_back(used);
        return;
      }
      int best = -1;
      uint32_t best_live = ~0u;
      for (uint64_t left = todo; left; left &= left-1) {
        const int g = std::countr_zero(left);
        uint32_t live = 0;
        for (uint32_t i = 0; i < W; i++) for (uint64_t m = scope_cells[g].w[i]; m; m &= m-1) {
          const uint32_t c = i*64 + std::countr_zero(m);
          live += (incidence[c] & ~todo) == 0;
        }
        if (live < best_live) { best_live = live; best = g; if (!live) return; }
      }
      for (uint32_t i = 0; i < W; i++) for (uint64_t m = scope_cells[best].w[i]; m; m &= m-1) {
        const uint32_t c = i*64 + std::countr_zero(m);
        if (incidence[c] & ~todo) continue;
        C next = used;
        next.add(c);
        go(todo & ~incidence[c], next);
        if (fail) return;
      }
    }
  };

  bool enumerate(uint32_t l, std::vector<C>& out, uint64_t max_out) {
    Enum e;
    const auto& sc = s.scopes[l];
    e.scope_cells.resize(sc.size());
    e.incidence.assign(s.cells, 0);
    for (uint32_t g = 0; g < sc.size(); g++) for (auto c : sc[g]) {
      if (!s.allowed[l][c]) continue;
      e.scope_cells[g].add(c);
      e.incidence[c] |= uint64_t(1) << g;
    }
    for (uint32_t c = 0; c < s.cells; c++) if (s.required[l][c]) {
      if (!e.incidence[c]) return out.clear(), true; // no constraint can take it
      e.required.push_back(c);
    }
    e.out = &out;
    e.budget = &budget;
    e.nodes = &st.enum_nodes;
    e.max_out = max_out;
    const uint64_t todo = sc.size() == 64 ? ~uint64_t(0) : (uint64_t(1) << sc.size()) - 1;
    e.go(todo, C{});
    return !e.fail;
  }

  // Join: choose placements label by label, always for the label with the
  // fewest placements left after removing those that overlap the cells taken
  // so far (forward checking). A branch dies as soon as some label has none.
  // The last label is forced to the uncovered cells.
  std::vector<std::vector<C>> base;                          // per label
  std::vector<std::unordered_set<C, CellsHash<W>>> member;   // per label, lazily
  std::vector<uint8_t> free_ok; // forced placement of this label always valid
  bool uniform_all = false;
  struct Live { uint32_t label; const std::vector<C>* dom; };
  std::vector<std::vector<std::vector<C>>> bufs; // per depth, per remaining label
  uint32_t probe_base = 0;
  static constexpr double estimate_slack = 4;

  bool forced_valid(uint32_t l, const C& occupied) {
    if (free_ok[l]) return true;
    if (member[l].empty()) member[l].insert(base[l].begin(), base[l].end());
    return member[l].count(all.minus(occupied));
  }

  // When every label has the same constraints, all of size k, and the forced
  // label may take any cell no other label is forced into, the uncovered cells
  // always form a valid placement for it: each constraint has k cells and the
  // other k-1 labels take exactly one each. The last four free labels are then
  // counted with bit-parallel compatibility masks over 64 placements at a time.
  struct Option { C cells; uint64_t compatible; };
  std::array<std::vector<Option>,3> rows;
  uint64_t finish(std::array<const std::vector<C>*,4> avail) {
    std::sort(avail.begin(), avail.end(), [](auto x, auto y){ return x->size() < y->size(); });
    const auto& last = *avail[3];
    uint64_t total = 0;
    for (size_t start = 0; start < last.size(); start += 64) {
      const size_t len = std::min<size_t>(64, last.size()-start);
      for (uint32_t i = 0; i < 3; i++) {
        rows[i].clear();
        if (budget < avail[i]->size()*len) { exhausted = true; return 0; }
        budget -= avail[i]->size()*len;
        for (const auto& p : *avail[i]) {
          uint64_t mask = 0;
          for (size_t j = 0; j < len; j++) mask |= uint64_t(p.disjoint(last[start+j])) << j;
          if (mask) rows[i].push_back({p, mask});
        }
      }
      for (const auto& x : rows[0]) for (const auto& y : rows[1]) {
        const uint64_t shared = x.compatible & y.compatible;
        if (!shared || !x.cells.disjoint(y.cells)) continue;
        if (budget < rows[2].size()+1) { exhausted = true; return 0; }
        budget -= rows[2].size()+1;
        st.join_nodes++;
        const C used = x.cells | y.cells;
        for (const auto& z : rows[2]) if (used.disjoint(z.cells))
          total += std::popcount(shared & z.compatible);
      }
    }
    return total;
  }

  // Which label to branch on, and whether the finish applies here.
  uint32_t pick(const std::vector<Live>& live) const {
    uint32_t best = 0;
    for (uint32_t i = 1; i < live.size(); i++)
      if (live[i].dom->size() < live[best].dom->size()) best = i;
    return best;
  }
  int32_t finish_forced(const std::vector<Live>& live) const {
    if (!uniform_all || live.size() != 5) return -1;
    int32_t f = -1;
    for (uint32_t i = 0; i < 5; i++) if (free_ok[live[i].label])
      if (f == -1 || live[i].dom->size() > live[f].dom->size()) f = i;
    return f;
  }
  // Filter the other labels' placements by p into bufs[depth]; false if one empties
  bool branch(uint32_t depth, const std::vector<Live>& live, uint32_t skip, const C& p,
      std::vector<Live>& next) {
    // Sized up front: children keep pointers into these buffers
    auto& b = bufs[depth];
    next.clear();
    for (uint32_t i = 0; i < live.size(); i++) if (i != skip) {
      const auto& d = *live[i].dom;
      if (budget < d.size()) { exhausted = true; return false; }
      budget -= d.size();
      b[i].clear();
      for (const auto& q : d) if (p.disjoint(q)) b[i].push_back(q);
      if (b[i].empty()) return false;
      next.push_back({live[i].label, &b[i]});
    }
    return true;
  }

  // Two labels can swap placements and leave the same cells to the same
  // remaining labels; the subtree below depends only on that pair, since the
  // filtered domains and so every later choice are functions of it.
  struct Key {
    C occupied; uint64_t labels;
    bool operator==(const Key& o) const { return labels == o.labels && occupied == o.occupied; }
  };
  struct KeyHash {
    size_t operator()(const Key& k) const { return CellsHash<W>()(k.occupied) ^ (k.labels * 0x9e3779b97f4a7c15ULL); }
  };
  std::unordered_map<Key, mpz_class, KeyHash> memo;
  static constexpr size_t max_memo = 1u << 20;

  mpz_class join(uint32_t depth, const C& occupied, const std::vector<Live>& live) {
    if (exhausted) return 0;
    st.join_nodes++;
    if (live.size() == 1) return forced_valid(live[0].label, occupied) ? 1 : 0;
    Key key{occupied, 0};
    for (const auto& l : live) key.labels |= uint64_t(1) << l.label;
    if (auto it = memo.find(key); it != memo.end()) { st.memo_hits++; return it->second; }
    mpz_class total = 0;
    if (const int32_t f = finish_forced(live); f >= 0) {
      std::array<const std::vector<C>*,4> four;
      for (uint32_t i = 0, j = 0; i < 5; i++) if ((int32_t)i != f) four[j++] = live[i].dom;
      static_assert(sizeof(unsigned long) == 8);
      if (const uint64_t t = finish(four)) total = mpz_class((unsigned long)t);
    } else {
      const uint32_t b = pick(live);
      std::vector<Live> next;
      for (const auto& p : *live[b].dom) {
        if (branch(depth, live, b, p, next)) total += join(depth+1, occupied | p, next);
        if (exhausted) return 0;
      }
    }
    if (!exhausted && memo.size() < max_memo) memo.emplace(key, total);
    return total;
  }

  // Knuth's estimator: follow random paths and weigh each node's work by the
  // product of the branching factors above it. Averaged over probes this is an
  // unbiased estimate of the join's total work, used to decline up front.
  double estimate(const std::vector<Live>& root, uint32_t probes) {
    uint64_t rng = 0x2545F4914F6CDD1DULL;
    auto next_rand = [&]{ rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; };
    const uint64_t saved = budget;
    double sum = 0;
    for (uint32_t pr = 0; pr < probes && !exhausted; pr++) {
      std::vector<Live> live = root, next;
      double weight = 1, work = 0;
      for (uint32_t depth = 0; live.size() > 1; depth++) {
        if (finish_forced(live) >= 0) {
          std::array<double,4> sz{};
          const int32_t f = finish_forced(live);
          for (uint32_t i = 0, j = 0; i < 5; i++) if ((int32_t)i != f) sz[j++] = live[i].dom->size();
          std::sort(sz.begin(), sz.end());
          work += weight * (sz[0]*sz[1] + (sz[0]+sz[1]+sz[2])*std::min(64.0, sz[3]));
          break;
        }
        const uint32_t b = pick(live);
        double others = 0;
        for (uint32_t i = 0; i < live.size(); i++) if (i != b) others += live[i].dom->size();
        const double fan = live[b].dom->size();
        work += weight * fan * others;
        weight *= fan;
        const C& p = (*live[b].dom)[next_rand() % live[b].dom->size()];
        // Probe paths use the second half of the buffers
        if (!branch(probe_base + depth, live, b, p, next)) break;
        live = next;
      }
      sum += work;
    }
    budget = saved;
    exhausted = false;
    return sum / probes;
  }

  bool uniform() const {
    const uint32_t k = s.scopes.size();
    if (k < 5) return false;
    auto norm = [](std::vector<std::vector<uint32_t>> v){ std::sort(v.begin(), v.end()); return v; };
    const auto first = norm(s.scopes[0]);
    if (first.empty()) return false;
    for (const auto& g : first) if (g.size() != k) return false;
    std::vector<uint8_t> in_scope(s.cells);
    for (const auto& g : first) for (auto c : g) in_scope[c] = 1;
    if (std::count(in_scope.begin(), in_scope.end(), 0)) return false;
    for (uint32_t l = 1; l < k; l++) if (norm(s.scopes[l]) != first) return false;
    return true;
  }

public:
  Counter(const Structure& _s, uint64_t _budget, Stats& _st) : s(_s), budget(_budget), st(_st) {}

  std::optional<mpz_class> count(uint64_t max_domain) {
    const uint32_t k = s.scopes.size();
    for (uint32_t c = 0; c < s.cells; c++) all.add(c);
    base.resize(k);
    for (uint32_t l = 0; l < k; l++) {
      if (!enumerate(l, base[l], max_domain)) return {};
      st.domains.push_back(base[l].size());
      if (base[l].empty()) return mpz_class(0);
    }
    member.resize(k);
    free_ok.assign(k, 0);
    uniform_all = uniform();
    if (uniform_all) for (uint32_t l = 0; l < k; l++) {
      // The forced label must be able to take every cell not forced elsewhere
      bool ok = true;
      for (uint32_t c = 0; c < s.cells; c++) {
        bool forced_elsewhere = false;
        for (uint32_t o = 0; o < k; o++) if (o != l && s.required[o][c]) forced_elsewhere = true;
        if (!forced_elsewhere && !s.allowed[l][c]) ok = false;
      }
      free_ok[l] = ok;
    }
    st.fast_join = uniform_all && std::count(free_ok.begin(), free_ok.end(), 1);
    std::vector<Live> root;
    for (uint32_t l = 0; l < k; l++) root.push_back({l, &base[l]});
    probe_base = k+1;
    bufs.assign(2*(k+1), std::vector<std::vector<C>>(k));
    // The estimate ignores the memo, so it only rules out hopeless joins;
    // the budget still bounds the actual work.
    st.estimate = estimate(root, 64);
    if (st.estimate > estimate_slack * (double)budget) return {};
    mpz_class total = join(0, C{}, root);
    if (exhausted) return {};
    return total;
  }
};
} // namespace detail

// Over budget: returns nullopt and the caller counts normally.
inline std::optional<mpz_class> count(const Structure& s, uint64_t budget, uint64_t max_domain, Stats& st) {
  if (s.contradictory) return mpz_class(0);
  if (s.cells <= 64) return detail::Counter<1>(s, budget, st).count(max_domain);
  if (s.cells <= 128) return detail::Counter<2>(s, budget, st).count(max_domain);
  return detail::Counter<4>(s, budget, st).count(max_domain);
}

} // namespace GanakInt::placement
