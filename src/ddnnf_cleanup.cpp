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

// ddnnf-cleanup: post-process a d4 .nnf circuit from `ganak --compile`.
//
// The streaming compiler writes a valid-but-loose file: root declared first (d4
// convention), but with the internal node ids and possibly unreachable ("dead")
// nodes (orphaned when an AND short-circuited to FALSE). Harmless to a root count,
// but a strict d4 reader wants a clean circuit. Cleanup steps:
//
//   1. STRICT-DECOMP REPAIR: when an AND node's children share a variable (a
//      known artifact of conflict-learned clauses bridging vars across the
//      analyzer's component split -- the analyzer only sees irreducible
//      clauses), scrub the shared var from one child's subtree. We pick the
//      child that DOESN'T force the var (the other child unconditionally
//      forces it, so the shared mentions are redundant in this AND context).
//      Cloning protects DAG-shared nodes; arcs with the wrong-polarity lit
//      become FALSE (dead under this AND anyway). Model count and model set
//      are preserved (each AND-of-siblings model is unchanged).
//   2. DEAD-NODE DROP + RENUMBER: BFS from the root, keep only reachable nodes,
//      renumber to root=1 contiguous (breadth-first).
//
// Final output is the classic two-section d4 file (all declarations, then all
// arcs). Pass --no-strict-decomp to skip step 1.
//
// Usage:  ddnnf-cleanup [--no-strict-decomp] <in.nnf> [out.nnf]
//                                                       ("-" or omitted => stdout)

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "ddnnf_io.hpp"

namespace {

using ddnnf_io::Arc;

[[noreturn]] void die(const std::string& msg) {
  std::cerr << "ddnnf-cleanup: " << msg << std::endl;
  exit(EXIT_FAILURE);
}

// Compute, for each reachable node, the set of vars appearing anywhere in its
// subtree (arc lits + child subtree vars). Bottom-up via post-order DFS.
void compute_subtree_vars(
    int root,
    const std::unordered_map<int, char>& type,
    const std::unordered_map<int, std::vector<Arc>>& arcs,
    std::unordered_map<int, std::unordered_set<int>>& out) {
  std::vector<int> stack = {root};
  std::vector<int> order;
  std::unordered_set<int> seen;
  while (!stack.empty()) {
    int n = stack.back(); stack.pop_back();
    if (!seen.insert(n).second) continue;
    order.push_back(n);
    auto it = arcs.find(n);
    if (it == arcs.end()) continue;
    for (const auto& a : it->second) stack.push_back(a.child);
  }
  // post-order processing: leaves first
  std::reverse(order.begin(), order.end());
  // Build child-->parents reverse map so we can process in dependency order via
  // topological-ish iteration. Simpler: just iterate twice (DAG depth is small in
  // practice) or use memoized recursion. We do memoized recursion to keep it simple.
  std::unordered_map<int, bool> in_progress;
  std::function<const std::unordered_set<int>&(int)> rec =
      [&](int n) -> const std::unordered_set<int>& {
    auto it = out.find(n);
    if (it != out.end()) return it->second;
    out[n] = {};   // placeholder against cycles (none expected, but safe)
    auto& vs = out[n];
    auto ait = arcs.find(n);
    if (ait != arcs.end()) {
      for (const auto& a : ait->second) {
        for (int l : a.lits) vs.insert(std::abs(l));
        const auto& cvs = rec(a.child);
        vs.insert(cvs.begin(), cvs.end());
      }
    }
    return vs;
  };
  rec(root);
  (void)type; (void)order; (void)in_progress;
}

// Compute, for each reachable node, the set of literals true in EVERY model of
// the subtree (the node's "forced" lits). nullopt = node is FALSE (no models).
// OR: intersection of (arc_lits + child_forced). AND: union (FALSE if contradictory).
void compute_forced(
    int root,
    const std::unordered_map<int, char>& type,
    const std::unordered_map<int, std::vector<Arc>>& arcs,
    std::unordered_map<int, std::unordered_set<int>>& forced,
    std::unordered_set<int>& is_false) {
  std::function<bool(int)> rec = [&](int n) -> bool {
    if (forced.count(n) || is_false.count(n)) return !is_false.count(n);
    char t = type.at(n);
    if (t == 't') { forced[n] = {}; return true; }
    if (t == 'f') { is_false.insert(n); return false; }
    auto it = arcs.find(n);
    const auto& kids = it->second;
    if (t == 'o') {
      std::vector<std::unordered_set<int>> sets;
      for (const auto& a : kids) {
        if (!rec(a.child)) continue;
        std::unordered_set<int> s(a.lits.begin(), a.lits.end());
        for (int l : forced[a.child]) s.insert(l);
        sets.push_back(std::move(s));
      }
      if (sets.empty()) { is_false.insert(n); return false; }
      auto r = sets[0];
      for (size_t i = 1; i < sets.size(); i++) {
        std::unordered_set<int> next;
        for (int l : r) if (sets[i].count(l)) next.insert(l);
        r = std::move(next);
      }
      forced[n] = std::move(r);
      return true;
    } else { // 'a'
      std::unordered_set<int> r;
      for (const auto& a : kids) {
        if (!rec(a.child)) { is_false.insert(n); return false; }
        for (int l : a.lits) r.insert(l);
        for (int l : forced[a.child]) r.insert(l);
      }
      for (int l : r) if (r.count(-l)) { is_false.insert(n); return false; }
      forced[n] = std::move(r);
      return true;
    }
  };
  rec(root);
}

// Count parents of each reachable node (so we know whether modifying it in place
// is safe -- nodes with multiple parents must be cloned before modification).
void compute_parents(
    int root,
    const std::unordered_map<int, std::vector<Arc>>& arcs,
    std::unordered_map<int, int>& parent_count) {
  std::unordered_set<int> seen;
  std::vector<int> stack = {root};
  parent_count[root] = 0;   // root has no parent
  while (!stack.empty()) {
    int n = stack.back(); stack.pop_back();
    if (!seen.insert(n).second) continue;
    auto it = arcs.find(n);
    if (it == arcs.end()) continue;
    for (const auto& a : it->second) {
      parent_count[a.child]++;
      stack.push_back(a.child);
    }
  }
}

// Scrub `var` from the subtree rooted at `nid`, knowing that an enclosing AND
// forces var to `polarity`. Returns the new root id (== nid if modified in
// place, else a fresh clone id). Caller passes max_id (mutable) for cloning.
//
// In each arc within the subtree:
//   - Drop the redundant (forced-polarity) literal from arc lits.
//   - If the OPPOSITE-polarity literal appears, the arc is dead under the
//     enclosing AND -- drop the arc entirely.
//   - Recurse into the child to scrub deeper occurrences.
// An OR whose arcs all get dropped becomes FALSE. Nodes with multiple parents
// (DAG sharing) get cloned before modification so other parents stay correct.
int scrub_var(
    int nid, int var, bool polarity,
    std::unordered_map<int, char>& type,
    std::unordered_map<int, std::vector<Arc>>& arcs,
    const std::unordered_map<int, int>& parent_count,
    std::unordered_map<long long, int>& memo,  // (nid, var, polarity) -> new id
    int& max_id) {
  long long key = ((long long)nid << 32) | (long long)(var * 2 + (polarity ? 1 : 0));
  auto mit = memo.find(key);
  if (mit != memo.end()) return mit->second;

  char t = type.at(nid);
  if (t == 't' || t == 'f') { memo[key] = nid; return nid; }

  int forced_lit = polarity ? var : -var;
  int opposite_lit = -forced_lit;

  std::vector<Arc> new_arcs;
  bool changed = false;
  for (const auto& a : arcs.at(nid)) {
    // arc dead under enclosing AND if it constrains var to opposite polarity
    if (std::find(a.lits.begin(), a.lits.end(), opposite_lit) != a.lits.end()) {
      changed = true;
      continue;
    }
    // strip the redundant forced literal
    Arc na{a.child, {}};
    for (int l : a.lits) if (l != forced_lit) na.lits.push_back(l);
    if (na.lits.size() != a.lits.size()) changed = true;
    // recurse into the child
    int new_c = scrub_var(a.child, var, polarity, type, arcs, parent_count, memo, max_id);
    if (new_c != na.child) changed = true;
    na.child = new_c;
    new_arcs.push_back(std::move(na));
  }

  if (!changed) { memo[key] = nid; return nid; }

  // Pick: modify in place (only if exactly one parent) or clone.
  auto pcit = parent_count.find(nid);
  int p_count = (pcit == parent_count.end()) ? 1 : pcit->second;
  int target = (p_count <= 1) ? nid : (++max_id);
  type[target] = t;
  if (t == 'o' && new_arcs.empty()) {
    type[target] = 'f';
    arcs[target] = {};
  } else {
    arcs[target] = std::move(new_arcs);
  }
  memo[key] = target;
  return target;
}

// Top-level driver: iterate until no AND violations remain.
int cleanup_decomp(
    int root,
    std::unordered_map<int, char>& type,
    std::unordered_map<int, std::vector<Arc>>& arcs) {
  int max_id = 0;
  for (const auto& p : type) max_id = std::max(max_id, p.first);
  int iterations = 0;
  while (true) {
    if (++iterations > 100) {
      std::cerr << "ddnnf-cleanup: WARN: strict-decomp cleanup did not converge in 100 iterations"
                << std::endl;
      break;
    }
    std::unordered_map<int, std::unordered_set<int>> subtree_vars;
    compute_subtree_vars(root, type, arcs, subtree_vars);
    std::unordered_map<int, std::unordered_set<int>> forced;
    std::unordered_set<int> is_false;
    compute_forced(root, type, arcs, forced, is_false);
    std::unordered_map<int, int> parent_count;
    compute_parents(root, arcs, parent_count);

    // Find first AND violation
    bool fixed = false;
    for (auto& tpair : type) {
      int nid = tpair.first;
      if (tpair.second != 'a') continue;
      auto ait = arcs.find(nid);
      if (ait == arcs.end()) continue;
      auto& kids = ait->second;
      // pairwise check for shared vars
      for (size_t i = 0; i < kids.size() && !fixed; i++) {
        for (size_t j = i + 1; j < kids.size() && !fixed; j++) {
          int ci = kids[i].child, cj = kids[j].child;
          auto svi = subtree_vars.find(ci);
          auto svj = subtree_vars.find(cj);
          if (svi == subtree_vars.end() || svj == subtree_vars.end()) continue;
          int shared_var = 0;
          for (int v : svi->second) {
            if (svj->second.count(v)) { shared_var = v; break; }
          }
          if (shared_var == 0) continue;
          // Decide which child forces the var
          auto fi = forced.find(ci);
          auto fj = forced.find(cj);
          int pi = 0, pj = 0;
          if (fi != forced.end()) {
            if (fi->second.count(shared_var))  pi = +1;
            if (fi->second.count(-shared_var)) pi = -1;
          }
          if (fj != forced.end()) {
            if (fj->second.count(shared_var))  pj = +1;
            if (fj->second.count(-shared_var)) pj = -1;
          }
          int forced_child = 0, scrub_child = 0;
          bool polarity = true;
          if (pi != 0) { forced_child = ci; scrub_child = cj; polarity = (pi > 0); }
          else if (pj != 0) { forced_child = cj; scrub_child = ci; polarity = (pj > 0); }
          else {
            std::cerr << "ddnnf-cleanup: WARN: AND " << nid << " shares var " << shared_var
                      << " but neither child forces it; leaving as-is" << std::endl;
            continue;
          }
          // Scrub
          std::unordered_map<long long, int> memo;
          int new_scrub = scrub_var(scrub_child, shared_var, polarity, type, arcs,
                                     parent_count, memo, max_id);
          if (new_scrub != scrub_child) {
            // Patch THIS AND's arc to point at the new (cloned) child
            for (auto& a : kids) if (a.child == scrub_child) a.child = new_scrub;
          }
          (void)forced_child;
          fixed = true;
        }
      }
      if (fixed) break;
    }
    if (!fixed) break;
  }
  return root;
}

} // namespace

int main(int argc, char** argv) {
  bool do_strict_decomp = true;
  int argi = 1;
  if (argi < argc && std::string(argv[argi]) == "--no-strict-decomp") {
    do_strict_decomp = false;
    argi++;
  }
  if (argc - argi < 1 || argc - argi > 2) {
    std::cerr << "Usage: " << argv[0]
              << " [--no-strict-decomp] <in.nnf> [out.nnf]"
                 "   (\"-\"/omitted out => stdout)" << std::endl;
    return EXIT_FAILURE;
  }
  const std::string in_path = argv[argi];
  const std::string out_path = (argc - argi == 2) ? argv[argi + 1] : "-";

  std::ifstream in(in_path);
  if (!in.good()) die("cannot open input file: " + in_path);

  std::unordered_map<int, char> type;             // node id -> 'f'|'t'|'a'|'o'
  std::unordered_map<int, std::vector<Arc>> arcs; // node id -> outgoing arcs
  std::vector<int> decl_order;                    // declaration order (unused; we renumber by BFS)
  std::string err;
  int root = ddnnf_io::parse_nnf(in, type, arcs, decl_order, err);
  if (root == -1) die(err);

  if (do_strict_decomp) root = cleanup_decomp(root, type, arcs);

  // BFS renumber from the root: root => 1, then children in arc order. Unvisited
  // (dead) nodes get no id and are dropped.
  std::unordered_map<int, int> newid;
  std::vector<int> order;                           // new-id order -> old id
  std::queue<int> q;
  newid[root] = 1;
  order.push_back(root);
  q.push(root);
  while (!q.empty()) {
    int nid = q.front();
    q.pop();
    auto it = arcs.find(nid);
    if (it == arcs.end()) continue;
    for (const auto& a : it->second) {
      if (newid.find(a.child) == newid.end()) {
        newid[a.child] = (int)order.size() + 1;
        order.push_back(a.child);
        q.push(a.child);
      }
    }
  }

  std::ofstream fout;
  std::ostream* out = &std::cout;
  if (out_path != "-") {
    fout.open(out_path);
    if (!fout.good()) die("cannot open output file: " + out_path);
    out = &fout;
  }

  // Section 1: declarations in new-id order.
  for (int oldid : order) *out << type[oldid] << " " << newid[oldid] << " 0\n";
  // Section 2: arcs in new-id order, with ids remapped.
  for (int oldid : order) {
    for (const auto& a : arcs[oldid]) {
      *out << newid[oldid] << " " << newid[a.child];
      for (int l : a.lits) *out << " " << l;
      *out << " 0\n";
    }
  }
  out->flush();

  std::cerr << "ddnnf-cleanup: declared=" << type.size()
            << " reachable=" << order.size()
            << " dropped=" << (type.size() - order.size()) << std::endl;
  return EXIT_SUCCESS;
}
