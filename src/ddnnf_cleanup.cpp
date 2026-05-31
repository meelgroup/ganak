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
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <utility>
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
    const std::vector<std::vector<Arc>>& arcs,
    int& max_var, size_t& words) {
  max_var = 0;
  for (const auto& kids : arcs) {
    for (const auto& a : kids) {
      for (int l : a.lits) {
        int v = (l > 0 ? l : -l);
        if (v > max_var) max_var = v;
      }
    }
  }
  words = (size_t)(max_var / 64) + 1;
}

// Iterative post-order DFS from `root`: emits children before parents.
// `visited` must be sized to at least max_id+1 (cleared by the caller). At
// 14M-node circuits this matters -- a recursive version blows the 8 MB stack.
void post_order(
    int root,
    const std::vector<std::vector<Arc>>& arcs,
    std::vector<char>& visited,   // 0 = unseen, 1 = entered, 2 = emitted
    std::vector<int>& out) {
  // Stack of (node, next_child_index_to_visit). When we revisit a node and its
  // next index is past its arc count, we emit it.
  std::vector<std::pair<int, int>> stack;
  stack.reserve(64);
  if (visited[root]) return;
  visited[root] = 1;
  stack.emplace_back(root, 0);
  while (!stack.empty()) {
    auto& top = stack.back();
    int n = top.first;
    int idx = top.second;
    int sz = (n < (int)arcs.size()) ? (int)arcs[n].size() : 0;
    if (idx >= sz) {
      out.push_back(n);
      visited[n] = 2;
      stack.pop_back();
      continue;
    }
    int c = arcs[n][idx].child;
    top.second = idx + 1;
    if (!visited[c]) {
      visited[c] = 1;
      stack.emplace_back(c, 0);
    }
  }
}

// Compute, for each reachable node, the set of vars appearing anywhere in its
// subtree (arc lits + child subtree vars). Iterative bottom-up over post-order.
// `out` is sized to max_id+1; entry per reachable node.
void compute_subtree_vars(
    const std::vector<int>& post,
    const std::vector<std::vector<Arc>>& arcs,
    size_t words,
    std::vector<VarSet>& out) {
  for (int n : post) {
    VarSet& vs = out[n];
    vs.resize(0, words);
    if (n >= (int)arcs.size()) continue;
    for (const auto& a : arcs[n]) {
      for (int l : a.lits) vs.set(l > 0 ? l : -l);
      vs.or_with(out[a.child]);
    }
  }
}

// Compute, for each reachable node, the literals true in EVERY model of the
// subtree (the node's "forced" lits). is_false[nid] = 1 => node is FALSE (no
// models). OR: intersection of (arc_lits + child_forced) across non-FALSE arcs.
// AND: union; FALSE if union contradicts. Iterative bottom-up over post-order.
void compute_forced(
    const std::vector<int>& post,
    const std::vector<char>& type,             // node id -> 'f'|'t'|'a'|'o'|0
    const std::vector<std::vector<Arc>>& arcs,
    size_t words,
    std::vector<LitSet>& forced,
    std::vector<char>& is_false) {
  LitSet scratch;
  for (int n : post) {
    char t = type[n];
    if (t == 't') { forced[n].resize(words); continue; }
    if (t == 'f') { is_false[n] = 1; continue; }
    const auto& kids = arcs[n];
    if (t == 'o') {
      bool any = false;
      LitSet& r = forced[n];
      for (const auto& a : kids) {
        if (is_false[a.child]) continue;
        // scratch = forced[child] | arc_lits
        scratch = forced[a.child];
        for (int l : a.lits) scratch.set(l);
        if (!any) { r = std::move(scratch); any = true; scratch.resize(words); }
        else      { r.intersect_with(scratch); }
      }
      if (!any) { is_false[n] = 1; r.pos.clear(); r.neg.clear(); }
    } else { // 'a'
      LitSet& r = forced[n];
      r.resize(words);
      bool dead = false;
      for (const auto& a : kids) {
        if (is_false[a.child]) { dead = true; break; }
        r.union_with(forced[a.child]);
        for (int l : a.lits) r.set(l);
      }
      if (dead || r.contradicts()) { is_false[n] = 1; r.pos.clear(); r.neg.clear(); }
    }
  }
}

// Count parents of each reachable node (so we know whether modifying it in place
// is safe -- nodes with multiple parents must be cloned before modification).
// parent_count must be sized to max_id+1, zero-initialized.
void compute_parents(
    int root,
    const std::vector<std::vector<Arc>>& arcs,
    std::vector<int>& parent_count) {
  std::vector<char> seen(parent_count.size(), 0);
  std::vector<int> stack = {root};
  while (!stack.empty()) {
    int n = stack.back(); stack.pop_back();
    if (seen[n]) continue;
    seen[n] = 1;
    if (n >= (int)arcs.size()) continue;
    for (const auto& a : arcs[n]) {
      parent_count[a.child]++;
      stack.push_back(a.child);
    }
  }
}

// Top-down pass: for each reachable node, compute the lits forced on EVERY path
// from the root to it ("ambient context"). An ancestor's OR-arc lit, taken on
// every reaching path, pins a var's polarity even if neither subtree of an AND
// below it forces the var on its own.
//
// Per edge (p -> c): contribution = ambient[p] ∪ arc.lits. For a DAG-shared
// child, the ambient is the INTERSECTION across all incoming edges (lits common
// to every reaching path). Iterating reverse-post-order means every parent has
// already pushed before we touch the child, so the meet finishes in one pass.
//
// ambient[c] is left in its default (empty) state for unreachable c; the caller
// only reads it for nodes that appeared in `post`.
void compute_ambient_forced(
    const std::vector<int>& post,
    int root,
    const std::vector<std::vector<Arc>>& arcs,
    size_t words,
    std::vector<LitSet>& ambient,
    std::vector<char>& ambient_set) {
  ambient[root].resize(words);
  ambient_set[root] = 1;
  LitSet contrib;
  contrib.resize(words);
  for (auto it = post.rbegin(); it != post.rend(); ++it) {
    int n = *it;
    if (!ambient_set[n]) continue;  // unreachable from root
    if (n >= (int)arcs.size()) continue;
    for (const auto& a : arcs[n]) {
      int c = a.child;
      // contrib = ambient[n] ∪ arc lits  (reuse scratch capacity)
      contrib = ambient[n];
      for (int l : a.lits) contrib.set(l);
      if (!ambient_set[c]) {
        ambient[c] = contrib;
        ambient_set[c] = 1;
      } else {
        ambient[c].intersect_with(contrib);
      }
    }
  }
}

// Scrub `var` from the subtree rooted at `nid`, knowing that an enclosing AND
// forces var to `polarity`. Returns the new root id (== nid if modified in
// place, else a fresh clone id). Iterative post-order walk over the subtree --
// recursive at 14M nodes blew the 8 MB stack.
//
// In each arc within the subtree:
//   - Drop the redundant (forced-polarity) literal from arc lits.
//   - If the OPPOSITE-polarity literal appears, the arc is dead under the
//     enclosing AND -- drop the arc entirely.
// An OR whose arcs all get dropped becomes FALSE. Nodes with multiple parents
// (DAG sharing) get cloned before modification so other parents stay correct.
int scrub_var(
    int nid, int var, bool polarity,
    std::vector<char>& type,
    std::vector<std::vector<Arc>>& arcs,
    const std::vector<int>& parent_count,
    int& max_id) {
  auto ensure_node_slot = [&](int id) {
    if ((int)type.size() <= id) {
      type.resize((size_t)id + 1, 0);
      arcs.resize((size_t)id + 1);
    }
  };
  // 1. Iterative post-order over the subtree rooted at nid (gives us children
  //    before parents so each parent can look up its (possibly remapped) kids).
  std::vector<int> post;
  {
    // The subtree can repeat nodes that are shared up in the DAG, so dedupe
    // per-walk. The subtree is typically a small fraction of the full DAG, so
    // a hash map is cheaper than a max_id-sized vector for this scratch.
    std::unordered_map<int, char> visited;
    std::vector<std::pair<int, int>> stack;
    stack.reserve(64);
    visited[nid] = 1;
    stack.emplace_back(nid, 0);
    while (!stack.empty()) {
      auto& top = stack.back();
      int n = top.first;
      int idx = top.second;
      char t = type[n];
      int sz = (t == 't' || t == 'f') ? 0 : (int)arcs[n].size();
      if (idx >= sz) {
        post.push_back(n);
        stack.pop_back();
        continue;
      }
      int c = arcs[n][idx].child;
      top.second = idx + 1;
      if (!visited.count(c)) {
        visited[c] = 1;
        stack.emplace_back(c, 0);
      }
    }
  }

  const int forced_lit = polarity ? var : -var;
  const int opposite_lit = -forced_lit;
  std::unordered_map<int, int> remap;  // old id -> (possibly cloned) new id

  // 2. Process children-first, building new arc lists and deciding modify-in-
  //    place vs clone for each node touched. Nodes whose arcs+kids didn't
  //    change keep their identity (remap[n] = n).
  for (int n : post) {
    char t = type[n];
    if (t == 't' || t == 'f') { remap[n] = n; continue; }
    const auto& orig = arcs[n];
    std::vector<Arc> new_arcs;
    bool changed = false;
    for (const auto& a : orig) {
      if (std::find(a.lits.begin(), a.lits.end(), opposite_lit) != a.lits.end()) {
        changed = true;
        continue;
      }
      Arc na{a.child, {}};
      for (int l : a.lits) if (l != forced_lit) na.lits.push_back(l);
      if (na.lits.size() != a.lits.size()) changed = true;
      auto rit = remap.find(a.child);
      int new_c = (rit == remap.end()) ? a.child : rit->second;
      if (new_c != na.child) { na.child = new_c; changed = true; }
      new_arcs.push_back(std::move(na));
    }
    if (!changed) { remap[n] = n; continue; }
    int p_count = (n < (int)parent_count.size()) ? parent_count[n] : 1;
    int target = (p_count <= 1) ? n : (++max_id);
    ensure_node_slot(target);
    type[target] = t;
    if (t == 'o' && new_arcs.empty()) {
      type[target] = 'f';
      arcs[target].clear();
    } else {
      arcs[target] = std::move(new_arcs);
    }
    remap[n] = target;
  }
  auto rit = remap.find(nid);
  return (rit == remap.end()) ? nid : rit->second;
}

// Top-level driver: iterate, batching as many independent fixes per pass as we
// can, until no AND violations remain (or we run out of iterations). Per-node
// data (subtree_vars, forced, is_false, parent_count) lives in plain vectors
// indexed by node id, sized to the current max_id+1 -- way cheaper than
// unordered_map<int, X> at multi-million-node scale.
int cleanup_decomp(
    int root,
    std::vector<char>& type,
    std::vector<std::vector<Arc>>& arcs) {
  int max_id = (int)type.size() - 1;
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
    // Size scratch to current max_id (which grows as scrub_var clones nodes).
    const size_t sz = (size_t)max_id + 1;
    if (type.size() < sz) { type.resize(sz, 0); arcs.resize(sz); }
    std::vector<char> visited(sz, 0);
    std::vector<int> post;
    post.reserve(sz);
    post_order(root, arcs, visited, post);

    std::vector<VarSet> subtree_vars(sz);
    compute_subtree_vars(post, arcs, words, subtree_vars);

    std::vector<LitSet> forced(sz);
    std::vector<char> is_false(sz, 0);
    compute_forced(post, type, arcs, words, forced, is_false);

    std::vector<int> parent_count(sz, 0);
    compute_parents(root, arcs, parent_count);

    std::vector<LitSet> ambient(sz);
    std::vector<char> ambient_set(sz, 0);
    compute_ambient_forced(post, root, arcs, words, ambient, ambient_set);

    // One pass: collect every (AND, scrub_child, var, polarity) we can act on,
    // then apply them. Skipping ANDs whose modified arcs would interfere this
    // pass keeps things sane; remaining violations re-surface next iteration.
    int fixes = 0;
    std::vector<char> touched(sz, 0);
    for (int nid : post) {
      if (type[nid] != 'a') continue;
      // NB: do NOT hold a reference into `arcs` across scrub_var() -- it may
      // call arcs.resize() to clone a node, which reallocates the outer vector
      // and invalidates any reference/pointer into it. Index via arcs[nid].
      for (size_t i = 0; i < arcs[nid].size(); i++) {
        for (size_t j = i + 1; j < arcs[nid].size(); j++) {
          int ci = arcs[nid][i].child, cj = arcs[nid][j].child;
          if (touched[ci] || touched[cj]) continue;
          int shared_var = subtree_vars[ci].first_intersect_var(subtree_vars[cj]);
          if (shared_var == 0) continue;
          int pi = 0, pj = 0;
          if (!is_false[ci]) {
            if (forced[ci].has(shared_var))  pi = +1;
            if (forced[ci].has(-shared_var)) pi = -1;
          }
          if (!is_false[cj]) {
            if (forced[cj].has(shared_var))  pj = +1;
            if (forced[cj].has(-shared_var)) pj = -1;
          }
          int scrub_child = 0;
          bool polarity = true;
          if (pi != 0)      { scrub_child = cj; polarity = (pi > 0); }
          else if (pj != 0) { scrub_child = ci; polarity = (pj > 0); }
          else {
            // Neither child forces v in its own subtree. But if v is pinned by
            // the ambient context (an ancestor OR-arc lit on every path to
            // this AND), we can still scrub. Sound because every model passing
            // through `nid` agrees with the ambient polarity; the matching lit
            // is redundant and the opposite-lit arcs are dead under this AND.
            int pa = 0;
            if (ambient[nid].has(shared_var))  pa = +1;
            if (ambient[nid].has(-shared_var)) pa = -1;
            if (pa != 0) { scrub_child = ci; polarity = (pa > 0); }
            else {
              std::cerr << "ddnnf-cleanup: WARN: AND " << nid << " shares var " << shared_var
                        << " but neither child nor ambient context forces it; leaving as-is"
                        << std::endl;
              continue;
            }
          }
          int new_scrub = scrub_var(scrub_child, shared_var, polarity, type, arcs,
                                    parent_count, max_id);
          if (new_scrub != scrub_child) {
            for (auto& a : arcs[nid]) if (a.child == scrub_child) a.child = new_scrub;
          }
          touched[scrub_child] = 1;
          if (new_scrub < (int)sz) touched[new_scrub] = 1;
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

  std::vector<char> type;                          // node id -> 'f'|'t'|'a'|'o' (0 = no node)
  std::vector<std::vector<Arc>> arcs;              // node id -> outgoing arcs
  std::vector<int> decl_order;                     // declaration order (unused; we renumber by BFS)
  std::string err;
  int root = ddnnf_io::parse_nnf(in, type, arcs, decl_order, err);
  if (root == -1) die(err);

  if (do_strict_decomp) root = cleanup_decomp(root, type, arcs);

  // BFS renumber from the root: root => 1, then children in arc order. Unvisited
  // (dead) nodes get no id and are dropped. `newid` indexed by node id; 0 = not
  // yet reached (real new ids start at 1 so this is a safe sentinel).
  std::vector<int> newid(type.size(), 0);
  std::vector<int> order;                           // new-id order -> old id
  std::queue<int> q;
  newid[root] = 1;
  order.push_back(root);
  q.push(root);
  while (!q.empty()) {
    int nid = q.front();
    q.pop();
    if (nid >= (int)arcs.size()) continue;
    for (const auto& a : arcs[nid]) {
      if (newid[a.child] == 0) {
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

  // `type` includes 0-sentinels for unused ids; count the real declarations.
  size_t declared = 0;
  for (char c : type) if (c != 0) declared++;
  std::cerr << "ddnnf-cleanup: declared=" << declared
            << " reachable=" << order.size()
            << " dropped=" << (declared - order.size()) << std::endl;
  return EXIT_SUCCESS;
}
