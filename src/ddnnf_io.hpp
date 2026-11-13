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

#pragma once

// In-memory circuit for the standalone d-DNNF tools, loaded from / saved to
// BDNNF (spec-bdnnf.md). The model is d4-shaped: node id -> 'f'|'t'|'a'|'o',
// arcs carry literals, ids 0/1 are FALSE/TRUE. An AND's units become one arc
// to TRUE, and each free var v becomes one shared synthetic node
// OR(v -> TRUE, -v -> TRUE) appended after the file's nodes; save() folds both
// back, so the graph passes of the tools need no BDNNF special cases.

#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "bdnnf.hpp"

namespace ddnnf_io {

struct Arc {
  int child;
  std::vector<int> lits;
};

struct Circuit {
  std::vector<char> type;
  std::vector<std::vector<Arc>> arcs;
  int root = 0;
  uint32_t nvars = 0;
  uint64_t props = 0;
  int nfile = 0;  // ids [2, nfile + 2) are nodes stored in the file
  bdnnf::KV footer_kv;
  bdnnf::Defs defs;  // passed through unchanged

  [[nodiscard]] bool is_free_or(int n) const {
    if (type[n] != 'o' || arcs[n].size() != 2) return false;
    const auto& a = arcs[n][0];
    const auto& b = arcs[n][1];
    return type[a.child] == 't' && type[b.child] == 't' && a.lits.size() == 1 &&
           b.lits.size() == 1 && a.lits[0] == -b.lits[0];
  }
};

inline bool load(const std::string& path, Circuit& c, std::string& err) {
  try {
    bdnnf::Reader r;
    r.open(path);
    if (r.nnodes + 2 > (uint64_t)INT32_MAX / 2) throw bdnnf::Error("circuit too large for this tool");
    c.nvars = r.nvars;
    c.nfile = (int)r.nnodes;
    const size_t n = (size_t)r.nnodes + 2;
    c.type.assign(n, 0);
    c.arcs.assign(n, {});
    c.type[0] = 'f';
    c.type[1] = 't';
    std::vector<int> free_or(r.nvars + 1, -1);
    auto free_node = [&](uint32_t v) {
      if (free_or[v] < 0) {
        free_or[v] = (int)c.type.size();
        c.type.push_back('o');
        c.arcs.push_back({Arc{1, {(int)v}}, Arc{1, {-(int)v}}});
      }
      return free_or[v];
    };
    r.scan([&](uint64_t id, const bdnnf::Node& nd) {
      std::vector<Arc> out;
      if (nd.kind == bdnnf::K_AND) {
        c.type[id] = 'a';
        for (uint64_t k : nd.kids) out.push_back(Arc{(int)k, {}});
        if (!nd.units.empty()) out.push_back(Arc{1, {nd.units.begin(), nd.units.end()}});
        for (uint32_t v : nd.free) out.push_back(Arc{free_node(v), {}});
      } else {
        c.type[id] = 'o';
        for (const auto& a : nd.arcs) out.push_back(Arc{(int)a.child, {a.lits.begin(), a.lits.end()}});
      }
      c.arcs[id] = std::move(out);
    });
    c.root = (int)r.root;
    c.props = r.props;
    c.footer_kv = r.footer_kv;
    c.defs = r.defs;
  } catch (const bdnnf::Error& e) {
    err = path + ": " + e.what();
    return false;
  }
  return true;
}

// Writes the part reachable from the root, children first, root last. The ids
// of `c` are not preserved.
inline bool save(const std::string& path, const Circuit& c, uint64_t props,
                 const std::string& producer, std::string& err,
                 const bdnnf::KV& footer_kv = {}) {
  const int nn = (int)c.type.size();
  // For an AND, the arcs to fold into free vars (one per var) and the children
  // to keep. Computed the same way during traversal and emission.
  auto and_parts = [&](int n, std::vector<int>& kids, std::vector<int>& units,
                       std::vector<int>& frees) {
    kids.clear();
    units.clear();
    frees.clear();
    for (const auto& a : c.arcs[n]) {
      units.insert(units.end(), a.lits.begin(), a.lits.end());
      if (c.is_free_or(a.child)) {
        const int v = std::abs(c.arcs[a.child][0].lits[0]);
        bool dup = false;
        for (int f : frees) dup |= (f == v);
        if (!dup) { frees.push_back(v); continue; }
      }
      if (c.type[a.child] == 't') continue;
      kids.push_back(a.child);
    }
  };
  std::vector<int64_t> newid(nn, -1);
  std::vector<int> kids, units, frees;
  try {
    bdnnf::Writer w;
    w.open(path, c.nvars, {{"producer", producer}});
    auto mapped = [&](int x) -> uint64_t {
      if (c.type[x] == 'f') return 0;
      if (c.type[x] == 't') return 1;
      return (uint64_t)newid[x];
    };
    auto emit = [&](int n) {
      if (c.type[n] == 'a') {
        and_parts(n, kids, units, frees);
        std::vector<uint64_t> mk;
        for (int k : kids) mk.push_back(mapped(k));
        newid[n] = (int64_t)w.add_and(mk, units, frees);
      } else {
        w.begin_or();
        for (const auto& a : c.arcs[n]) w.or_arc(a.lits, mapped(a.child));
        newid[n] = (int64_t)w.end_or();
      }
    };
    struct Frame { int n; std::vector<int> succ; size_t i; };
    std::vector<Frame> stack;
    std::vector<char> entered(nn, 0);
    auto push = [&](int n) {
      entered[n] = 1;
      Frame f{n, {}, 0};
      if (c.type[n] == 'a') {
        and_parts(n, kids, units, frees);
        f.succ = kids;
      } else if (c.type[n] == 'o') {
        for (const auto& a : c.arcs[n]) f.succ.push_back(a.child);
      }
      stack.push_back(std::move(f));
    };
    if (c.type[c.root] == 'a' || c.type[c.root] == 'o') push(c.root);
    while (!stack.empty()) {
      auto& top = stack.back();
      if (top.i < top.succ.size()) {
        const int ch = top.succ[top.i++];
        if (!entered[ch] && (c.type[ch] == 'a' || c.type[ch] == 'o')) push(ch);
        continue;
      }
      const int n = top.n;
      stack.pop_back();
      emit(n);
    }
    w.finish(mapped(c.root), props, footer_kv, &c.defs);
  } catch (const bdnnf::Error& e) {
    err = path + ": " + e.what();
    return false;
  }
  return true;
}

// d4 .nnf for legacy consumers: reachable nodes only, root = 1, ids
// contiguous in BFS order, declarations before arcs.
inline void write_d4(std::ostream& out, const Circuit& c) {
  std::vector<int> newid(c.type.size(), 0);
  std::vector<int> order;
  std::queue<int> q;
  newid[c.root] = 1;
  order.push_back(c.root);
  q.push(c.root);
  while (!q.empty()) {
    const int n = q.front();
    q.pop();
    for (const auto& a : c.arcs[n]) {
      if (newid[a.child] == 0) {
        newid[a.child] = (int)order.size() + 1;
        order.push_back(a.child);
        q.push(a.child);
      }
    }
  }
  for (int old : order) out << c.type[old] << " " << newid[old] << " 0\n";
  for (int old : order)
    for (const auto& a : c.arcs[old]) {
      out << newid[old] << " " << newid[a.child];
      for (int l : a.lits) out << " " << l;
      out << " 0\n";
    }
}

// Iterative post-order (children first) over the part reachable from `root`;
// fn(n) is called once per node. Recursion would blow the stack on big circuits.
template <class Fn>
void post_order(const Circuit& c, std::vector<char>& done, Fn&& fn) {
  std::vector<std::pair<int, int>> stack;
  if (done[c.root]) return;
  stack.emplace_back(c.root, 0);
  while (!stack.empty()) {
    auto& top = stack.back();
    const int n = top.first;
    if (done[n]) { stack.pop_back(); continue; }
    const int idx = top.second;
    if (idx >= (int)c.arcs[n].size()) {
      fn(n);
      done[n] = 1;
      stack.pop_back();
      continue;
    }
    const int ch = c.arcs[n][idx].child;
    top.second = idx + 1;
    if (!done[ch]) stack.emplace_back(ch, 0);
  }
}

struct VarSet {
  std::vector<uint64_t> w;
  void resize(size_t words) { w.assign(words, 0); }
  void set(int v) { w[v >> 6] |= 1ULL << (v & 63); }
  void or_with(const VarSet& o) {
    for (size_t i = 0; i < w.size(); i++) w[i] |= o.w[i];
  }
  [[nodiscard]] int first_intersect_var(const VarSet& o) const {
    for (size_t i = 0; i < w.size(); i++) {
      const uint64_t b = w[i] & o.w[i];
      if (b) return (int)(i * 64) + __builtin_ctzll(b);
    }
    return 0;
  }
  [[nodiscard]] bool operator==(const VarSet& o) const { return w == o.w; }
  [[nodiscard]] int first_diff_var(const VarSet& o) const {
    for (size_t i = 0; i < w.size(); i++) {
      const uint64_t b = w[i] ^ o.w[i];
      if (b) return (int)(i * 64) + __builtin_ctzll(b);
    }
    return 0;
  }
};

// Vars mentioned below each reachable node, arc literals included.
inline std::vector<VarSet> subtree_vars(const Circuit& c, std::vector<char>& reach) {
  int max_var = (int)c.nvars;
  for (const auto& kids : c.arcs)
    for (const auto& a : kids)
      for (int l : a.lits) max_var = std::max(max_var, std::abs(l));
  const size_t words = (size_t)(max_var / 64) + 1;
  std::vector<VarSet> sv(c.type.size());
  reach.assign(c.type.size(), 0);
  post_order(c, reach, [&](int n) {
    sv[n].resize(words);
    for (const auto& a : c.arcs[n]) {
      for (int l : a.lits) sv[n].set(std::abs(l));
      sv[n].or_with(sv[a.child]);
    }
  });
  return sv;
}

inline VarSet arc_vars(const std::vector<VarSet>& sv, const Arc& a) {
  VarSet v = sv[a.child];
  for (int l : a.lits) v.set(std::abs(l));
  return v;
}

// Every reachable AND's arcs (literals + subtree) are pairwise var-disjoint,
// and every arc's literals are var-disjoint from its child (an arc is the AND
// of its literals and its child).
inline std::pair<bool, std::string> check_decomposable(const Circuit& c) {
  std::vector<char> reach;
  const auto sv = subtree_vars(c, reach);
  for (int n = 0; n < (int)c.type.size(); n++) {
    if (!reach[n]) continue;
    for (const auto& a : c.arcs[n]) {
      VarSet lv;
      lv.resize(sv[n].w.size());
      for (int l : a.lits) {
        const int v = std::abs(l);
        if ((lv.w[v >> 6] >> (v & 63)) & 1) {
          return {false, "node " + std::to_string(n) + " has var " + std::to_string(v) +
                         " twice on one arc"};
        }
        lv.set(v);
      }
      if (const int v = lv.first_intersect_var(sv[a.child])) {
        std::ostringstream m;
        m << "node " << n << " arc literal on var " << v << " also occurs below its child "
          << a.child;
        return {false, m.str()};
      }
    }
    if (c.type[n] != 'a') continue;
    VarSet seen;
    seen.resize(sv[n].w.size());
    for (const auto& a : c.arcs[n]) {
      const VarSet cv = arc_vars(sv, a);
      if (const int v = seen.first_intersect_var(cv)) {
        std::ostringstream m;
        m << "AND " << n << " children share var " << v;
        return {false, m.str()};
      }
      seen.or_with(cv);
    }
  }
  return {true, ""};
}

// Every reachable OR's arcs mention the same vars, and the root mentions all
// of 1..nvars.
inline std::pair<bool, std::string> check_smooth(const Circuit& c) {
  std::vector<char> reach;
  const auto sv = subtree_vars(c, reach);
  for (int n = 0; n < (int)c.type.size(); n++) {
    if (!reach[n] || c.type[n] != 'o') continue;
    const VarSet* first = nullptr;
    VarSet f0;
    for (size_t i = 0; i < c.arcs[n].size(); i++) {
      if (c.type[c.arcs[n][i].child] == 'f') continue;
      const VarSet cv = arc_vars(sv, c.arcs[n][i]);
      if (!first) { f0 = cv; first = &f0; continue; }
      if (!(cv == *first)) {
        std::ostringstream m;
        m << "OR " << n << " arcs differ on var " << first->first_diff_var(cv);
        return {false, m.str()};
      }
    }
  }
  if (c.type[c.root] == 'f') return {true, ""};
  for (uint32_t v = 1; v <= c.nvars; v++)
    if (!((sv[c.root].w[v >> 6] >> (v & 63)) & 1)) {
      return {false, "root does not mention var " + std::to_string(v)};
    }
  return {true, ""};
}

// Syntactic determinism: any two arcs of a reachable OR carry complementary
// literals.
inline std::pair<bool, std::string> check_deterministic(const Circuit& c) {
  std::vector<char> reach(c.type.size(), 0);
  post_order(c, reach, [](int) {});
  for (int n = 0; n < (int)c.type.size(); n++) {
    if (!reach[n] || c.type[n] != 'o') continue;
    const auto& as = c.arcs[n];
    for (size_t i = 0; i < as.size(); i++)
      for (size_t j = i + 1; j < as.size(); j++) {
        bool clash = false;
        for (int l : as[i].lits)
          for (int m : as[j].lits) clash |= (l == -m);
        if (!clash) {
          std::ostringstream msg;
          msg << "OR " << n << " arcs " << i << " and " << j << " share no complementary literal";
          return {false, msg.str()};
        }
      }
  }
  return {true, ""};
}

}  // namespace ddnnf_io
