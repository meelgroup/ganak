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
//      are preserved (each AND-of-siblings model is unchanged). O(iters * N *
//      V) per pass with hash-set var/lit sets; can be slow + memory-hungry on
//      multi-million-node circuits, so pass --no-strict-decomp to skip it.
//   2. DEAD-NODE DROP + RENUMBER: BFS from the root, keep only reachable nodes,
//      renumber to root=1 contiguous (breadth-first).
//
// Final output is the classic two-section d4 file (all declarations, then all
// arcs).
//
// Usage:  ddnnf-cleanup [--no-strict-decomp] <in.nnf> [out.nnf]
//                                                  ("-" or omitted => stdout)

#include <algorithm>
#include <cstdint>
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

// Fixed-width bitset over [0, nvars]. Var v -> bit v. nvars is set globally per
// cleanup pass; the same word-count is used for every node's set so the inner
// loops are plain word-parallel ANDs/ORs (no hash-table overhead).
//
// For nvars=100 and 14M nodes, this is ~16 B per node (224 MB) instead of
// ~120 B+ for an unordered_set<int> (1.7 GB+), and intersection/union become
// 2-word loops instead of hashing.
struct VarSet {
  std::vector<uint64_t> w;
  void resize(int nvars_plus1, size_t words) { (void)nvars_plus1; w.assign(words, 0); }
  void set(int v) { w[v >> 6] |= 1ULL << (v & 63); }
  bool test(int v) const { return ((w[v >> 6] >> (v & 63)) & 1ULL) != 0; }
  void or_with(const VarSet& o) {
    const size_t n = w.size();
    for (size_t i = 0; i < n; i++) w[i] |= o.w[i];
  }
  bool intersects(const VarSet& o) const {
    const size_t n = w.size();
    for (size_t i = 0; i < n; i++) if (w[i] & o.w[i]) return true;
    return false;
  }
  // Return any var in (this & o), or 0 if none.
  int first_intersect_var(const VarSet& o) const {
    const size_t n = w.size();
    for (size_t i = 0; i < n; i++) {
      uint64_t b = w[i] & o.w[i];
      if (b) return (int)(i * 64) + __builtin_ctzll(b);
    }
    return 0;
  }
};

// Forced literals: track positive and negative bits per var separately. has(+v)
// = "subtree forces v=true", has(-v) = "subtree forces v=false". intersection
// (OR-node merge) and union (AND-node merge) are still word-parallel; an AND
// contradicts iff any var has both pos and neg bits set.
struct LitSet {
  std::vector<uint64_t> pos, neg;
  void resize(size_t words) { pos.assign(words, 0); neg.assign(words, 0); }
  void set(int l) {
    int v = (l > 0 ? l : -l);
    (l > 0 ? pos : neg)[v >> 6] |= 1ULL << (v & 63);
  }
  bool has(int l) const {
    int v = (l > 0 ? l : -l);
    const auto& bv = (l > 0 ? pos : neg);
    return ((bv[v >> 6] >> (v & 63)) & 1ULL) != 0;
  }
  void union_with(const LitSet& o) {
    const size_t n = pos.size();
    for (size_t i = 0; i < n; i++) { pos[i] |= o.pos[i]; neg[i] |= o.neg[i]; }
  }
  void intersect_with(const LitSet& o) {
    const size_t n = pos.size();
    for (size_t i = 0; i < n; i++) { pos[i] &= o.pos[i]; neg[i] &= o.neg[i]; }
  }
  bool contradicts() const {
    const size_t n = pos.size();
    for (size_t i = 0; i < n; i++) if (pos[i] & neg[i]) return true;
    return false;
  }
};

// Scan every arc literal to find the max variable; sets max_var and the
// per-node word count used for VarSet/LitSet allocation.
void measure_nvars(
    const std::unordered_map<int, std::vector<Arc>>& arcs,
    int& max_var, size_t& words) {
  max_var = 0;
  for (const auto& kv : arcs) {
    for (const auto& a : kv.second) {
      for (int l : a.lits) {
        int v = (l > 0 ? l : -l);
        if (v > max_var) max_var = v;
      }
    }
  }
  words = (size_t)(max_var / 64) + 1;
}

// Compute, for each reachable node, the set of vars appearing anywhere in its
// subtree (arc lits + child subtree vars). Bottom-up via memoized recursion.
void compute_subtree_vars(
    int root,
    const std::unordered_map<int, std::vector<Arc>>& arcs,
    size_t words,
    std::unordered_map<int, VarSet>& out) {
  std::function<const VarSet&(int)> rec = [&](int n) -> const VarSet& {
    auto it = out.find(n);
    if (it != out.end()) return it->second;
    VarSet& vs = out[n];
    vs.resize(0, words);
    auto ait = arcs.find(n);
    if (ait != arcs.end()) {
      for (const auto& a : ait->second) {
        for (int l : a.lits) vs.set(l > 0 ? l : -l);
        const VarSet& cvs = rec(a.child);
        vs.or_with(cvs);
      }
    }
    return vs;
  };
  rec(root);
}

// Compute, for each reachable node, the literals true in EVERY model of the
// subtree (the node's "forced" lits). is_false[nid] = true => node is FALSE
// (no models). OR: intersection of (arc_lits + child_forced) across arcs. AND:
// union; FALSE if union contradicts.
void compute_forced(
    int root,
    const std::unordered_map<int, char>& type,
    const std::unordered_map<int, std::vector<Arc>>& arcs,
    size_t words,
    std::unordered_map<int, LitSet>& forced,
    std::unordered_set<int>& is_false) {
  std::function<bool(int)> rec = [&](int n) -> bool {
    if (forced.count(n)) return true;
    if (is_false.count(n)) return false;
    char t = type.at(n);
    if (t == 't') { forced[n].resize(words); return true; }
    if (t == 'f') { is_false.insert(n); return false; }
    const auto& kids = arcs.at(n);
    if (t == 'o') {
      LitSet r;
      bool any = false;
      for (const auto& a : kids) {
        if (!rec(a.child)) continue;
        LitSet s = forced.at(a.child);     // copy (small, word-count words)
        for (int l : a.lits) s.set(l);
        if (!any) { r = std::move(s); any = true; }
        else      { r.intersect_with(s); }
      }
      if (!any) { is_false.insert(n); return false; }
      forced[n] = std::move(r);
      return true;
    } else { // 'a'
      LitSet r;
      r.resize(words);
      for (const auto& a : kids) {
        if (!rec(a.child)) { is_false.insert(n); forced.erase(n); return false; }
        const LitSet& cf = forced.at(a.child);
        r.union_with(cf);
        for (int l : a.lits) r.set(l);
      }
      if (r.contradicts()) { is_false.insert(n); return false; }
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

// Top-level driver: iterate, batching as many independent fixes per pass as we
// can, until no AND violations remain (or we run out of iterations).
int cleanup_decomp(
    int root,
    std::unordered_map<int, char>& type,
    std::unordered_map<int, std::vector<Arc>>& arcs) {
  int max_id = 0;
  for (const auto& p : type) max_id = std::max(max_id, p.first);
  constexpr int kMaxIters = 100;
  int iterations = 0;
  long long total_fixes = 0;
  int max_var; size_t words;
  measure_nvars(arcs, max_var, words);
  while (true) {
    if (++iterations > kMaxIters) {
      std::cerr << "ddnnf-cleanup: WARN: strict-decomp cleanup did not converge in "
                << kMaxIters << " iterations; output may still violate strict d-DNNF"
                << std::endl;
      break;
    }
    std::unordered_map<int, VarSet> subtree_vars;
    compute_subtree_vars(root, arcs, words, subtree_vars);
    std::unordered_map<int, LitSet> forced;
    std::unordered_set<int> is_false;
    compute_forced(root, type, arcs, words, forced, is_false);
    std::unordered_map<int, int> parent_count;
    compute_parents(root, arcs, parent_count);

    // One pass: collect every (AND, scrub_child, var, polarity) we can act on,
    // then apply them. Skipping ANDs whose modified arcs would interfere this
    // pass keeps things sane; remaining violations re-surface next iteration.
    int fixes = 0;
    std::unordered_set<int> touched_children;   // children we've already scrubbed this pass
    for (auto& tpair : type) {
      int nid = tpair.first;
      if (tpair.second != 'a') continue;
      auto ait = arcs.find(nid);
      if (ait == arcs.end()) continue;
      auto& kids = ait->second;
      for (size_t i = 0; i < kids.size(); i++) {
        for (size_t j = i + 1; j < kids.size(); j++) {
          int ci = kids[i].child, cj = kids[j].child;
          if (touched_children.count(ci) || touched_children.count(cj)) continue;
          auto svi = subtree_vars.find(ci);
          auto svj = subtree_vars.find(cj);
          if (svi == subtree_vars.end() || svj == subtree_vars.end()) continue;
          int shared_var = svi->second.first_intersect_var(svj->second);
          if (shared_var == 0) continue;
          auto fi = forced.find(ci);
          auto fj = forced.find(cj);
          int pi = 0, pj = 0;
          if (fi != forced.end()) {
            if (fi->second.has(shared_var))  pi = +1;
            if (fi->second.has(-shared_var)) pi = -1;
          }
          if (fj != forced.end()) {
            if (fj->second.has(shared_var))  pj = +1;
            if (fj->second.has(-shared_var)) pj = -1;
          }
          int scrub_child = 0;
          bool polarity = true;
          if (pi != 0)      { scrub_child = cj; polarity = (pi > 0); }
          else if (pj != 0) { scrub_child = ci; polarity = (pj > 0); }
          else {
            std::cerr << "ddnnf-cleanup: WARN: AND " << nid << " shares var " << shared_var
                      << " but neither child forces it; leaving as-is" << std::endl;
            continue;
          }
          std::unordered_map<long long, int> memo;
          int new_scrub = scrub_var(scrub_child, shared_var, polarity, type, arcs,
                                     parent_count, memo, max_id);
          if (new_scrub != scrub_child) {
            for (auto& a : kids) if (a.child == scrub_child) a.child = new_scrub;
          }
          touched_children.insert(scrub_child);
          touched_children.insert(new_scrub);
          fixes++;
        }
      }
    }
    total_fixes += fixes;
    if (fixes == 0) break;
  }
  std::cerr << "ddnnf-cleanup: strict-decomp iters=" << iterations
            << " fixes=" << total_fixes
            << " max_var=" << max_var << std::endl;
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
