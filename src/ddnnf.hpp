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

#include <vector>
#include <array>
#include <cstdint>
#include <string>
#include <fstream>

namespace GanakInt {

// Builds a (Decision-)d-DNNF circuit out of Ganak's search trace and writes it
// to disk in the d4 `.nnf` format. The mapping is the classic "DPLL with a
// trace = Decision-DNNF" (Huang & Darwiche, IJCAI 2005):
//   - a decision on a variable        -> OR node (its two arcs carry the
//                                         decision literal + implied literals)
//   - a component decomposition (AND)  -> AND node (arcs carry no literals)
//   - a free independent variable      -> OR(v -> T, -v -> T)   (factor 2)
//   - a SAT leaf / no remaining comps  -> TRUE leaf
//   - a conflict / UNSAT branch        -> FALSE leaf
//   - a component-cache hit            -> a shared edge to an existing node
//
// Counting the resulting circuit structurally (OR = sum of children, AND =
// product, TRUE = 1, FALSE = 0; literals on arcs ignored) reproduces exactly
// the value Ganak's internal count() returns, because every factor of two is
// made explicit as a free-variable node and every determined variable
// contributes a factor of one (it only appears as a literal on an arc).
class DDNNFCompiler {
public:
  enum NType : uint8_t { N_FALSE = 0, N_TRUE = 1, N_AND = 2, N_OR = 3 };
  struct Arc {
    int child;
    std::vector<int> lits; // DIMACS literals true on this arc
  };
  struct Node {
    NType type;
    std::vector<Arc> arcs;
  };

  DDNNFCompiler() {
    false_node = add_node(N_FALSE);
    true_node = add_node(N_TRUE);
  }

  int false_node = -1;
  int true_node = -1;
  int root = -1;
  int nvars = 0;

  std::vector<Node> nodes;

  int add_node(NType t) {
    nodes.push_back(Node{t, {}});
    return (int)nodes.size() - 1;
  }

  // AND of the given children (arcs carry no literals).
  // Simplifications: drop TRUE children; any FALSE child -> FALSE; empty -> TRUE;
  // single child -> that child (no spurious AND node).
  int mk_and(const std::vector<int>& comps) {
    std::vector<int> kids;
    kids.reserve(comps.size());
    for (int c : comps) {
      if (c == true_node) continue;
      if (c == false_node) return false_node;
      kids.push_back(c);
    }
    if (kids.empty()) return true_node;
    if (kids.size() == 1) return kids[0];
    int n = add_node(N_AND);
    for (int c : kids) nodes[n].arcs.push_back(Arc{c, {}});
    return n;
  }

  // OR of (literal-set -> child) arcs. Drop arcs to FALSE; empty -> FALSE.
  int mk_or(std::vector<Arc>&& arcs) {
    std::vector<Arc> kept;
    kept.reserve(arcs.size());
    for (auto& a : arcs) {
      if (a.child == false_node) continue;
      kept.push_back(std::move(a));
    }
    if (kept.empty()) return false_node;
    int n = add_node(N_OR);
    nodes[n].arcs = std::move(kept);
    return n;
  }

  // A free independent variable (unconstrained): contributes a factor of two.
  int mk_free_var(int dimacs_var) {
    int n = add_node(N_OR);
    nodes[n].arcs.push_back(Arc{true_node, {dimacs_var}});
    nodes[n].arcs.push_back(Arc{true_node, {-dimacs_var}});
    return n;
  }

  // Wrap a node under a single AND arc carrying `lits` (e.g. level-0 implied lits).
  int wrap_lits(int inner, const std::vector<int>& lits) {
    if (lits.empty()) return inner;
    if (inner == false_node) return false_node;
    int n = add_node(N_AND);
    nodes[n].arcs.push_back(Arc{inner, lits});
    return n;
  }

  // ---- per-search bookkeeping (not part of the emitted DAG) ----
  // children[lev][branch] = component node ids that are AND-ed into that branch
  std::vector<std::array<std::vector<int>, 2>> children;
  // left_lits[lev] = literals captured for the left branch of decision level lev
  std::vector<std::vector<int>> left_lits;
  // override_node[lev] >= 0: the SAT oracle solved this level; its node is a
  // witness leaf (set_override), bypassing the normal OR-from-children build.
  std::vector<int> override_node;

  void ensure_level(int lev) {
    if ((int)children.size() <= lev) {
      children.resize(lev + 1);
      left_lits.resize(lev + 1);
      override_node.resize(lev + 1, -1);
    }
  }
  void on_new_level(int lev) {
    ensure_level(lev);
    children[lev][0].clear();
    children[lev][1].clear();
    left_lits[lev].clear();
    override_node[lev] = -1;
  }
  void set_override(int lev, int node) {
    ensure_level(lev);
    override_node[lev] = node;
  }
  int take_override(int lev) {
    if (lev < 0 || (int)override_node.size() <= lev) return -1;
    int n = override_node[lev];
    override_node[lev] = -1;
    return n;
  }
  void add_child(int lev, bool branch, int node) {
    if (node == true_node) return;
    ensure_level(lev);
    children[lev][branch].push_back(node);
  }

  // Structural model count (debug cross-check only; may overflow on big counts).
  unsigned long long scount(int nid) const {
    std::vector<long long> memo(nodes.size(), -1);
    return scount_rec(nid, memo);
  }
  unsigned long long scount_rec(int nid, std::vector<long long>& memo) const {
    if (memo[nid] >= 0) return (unsigned long long)memo[nid];
    unsigned long long r;
    const Node& n = nodes[nid];
    if (n.type == N_FALSE) r = 0;
    else if (n.type == N_TRUE) r = 1;
    else if (n.type == N_OR) {
      r = 0;
      for (const auto& a : n.arcs) r += scount_rec(a.child, memo);
    } else {
      r = 1;
      for (const auto& a : n.arcs) r *= scount_rec(a.child, memo);
    }
    memo[nid] = (long long)r;
    return r;
  }

  size_t num_nodes() const { return nodes.size(); }
  size_t num_edges() const {
    size_t e = 0;
    for (const auto& n : nodes) e += n.arcs.size();
    return e;
  }

  // Write reachable-from-root nodes in the d4 .nnf format. The root is emitted
  // first, so it receives id 1.
  void write_d4(const std::string& fname) const {
    int r = (root < 0) ? false_node : root;
    std::vector<char> seen(nodes.size(), 0);
    std::vector<int> order;
    order.reserve(nodes.size());
    order.push_back(r);
    seen[r] = 1;
    for (size_t i = 0; i < order.size(); i++) {
      for (const auto& a : nodes[order[i]].arcs) {
        if (!seen[a.child]) { seen[a.child] = 1; order.push_back(a.child); }
      }
    }
    std::vector<int> outid(nodes.size(), 0);
    for (size_t i = 0; i < order.size(); i++) outid[order[i]] = (int)i + 1;

    std::ofstream o(fname);
    for (int n : order) {
      char c = nodes[n].type == N_FALSE ? 'f'
             : nodes[n].type == N_TRUE  ? 't'
             : nodes[n].type == N_AND   ? 'a'
                                        : 'o';
      o << c << " " << outid[n] << " 0\n";
    }
    for (int n : order) {
      for (const auto& a : nodes[n].arcs) {
        o << outid[n] << " " << outid[a.child];
        for (int l : a.lits) o << " " << l;
        o << " 0\n";
      }
    }
  }
};

}
