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

#include "common.hpp"

#include <vector>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <unistd.h>

namespace GanakInt {

// Builds a (Decision-)d-DNNF circuit from Ganak's search trace and STREAMS it to
// disk in d4 `.nnf` format ("DPLL trace = Decision-DNNF", Huang & Darwiche 2005):
//   decision          -> OR node (arcs carry decision + implied lits)
//   AND-decomposition -> AND node (no arc lits)
//   free var          -> OR(v->T, -v->T) (factor 2)
//   SAT leaf -> TRUE;  conflict -> FALSE;  cache hit -> shared edge
// Counting structurally (OR=sum, AND=product, T=1, F=0; arc lits ignored)
// reproduces count().
//
// STREAMING: the DAG is never fully in memory. Nodes are complete at creation and
// every arc points to a lower id (children exist first), so each node is emitted
// as it is built. Declarations stream to one temp file, arcs to another;
// write_d4() stitches them with the root declared first (d4 convention).
//
// The file may contain unreachable ("dead") nodes (mk_and short-circuits to FALSE,
// orphaning siblings); harmless to a root traversal, and `ddnnf-cleanup` removes
// them + renumbers root=1 for strict d4 consumers.
class DDNNFCompiler {
public:
  enum NType : uint8_t { N_FALSE = 0, N_TRUE = 1, N_AND = 2, N_OR = 3 };
  struct Arc {
    int child;
    std::vector<int> lits; // DIMACS literals true on this arc
  };

  explicit DDNNFCompiler(const std::string& target_fname) {
    open_temps(target_fname);
    false_node = emit(N_FALSE, {});
    true_node = emit(N_TRUE, {});
  }

  int false_node = -1;
  int true_node = -1;
  int root = -1;
  int nvars = 0;

  // AND of children (no arc lits): drop TRUE; any FALSE -> FALSE; empty -> TRUE;
  // single child -> that child.
  int mk_and(const std::vector<int>& comps) {
    std::vector<Arc> arcs;
    arcs.reserve(comps.size());
    for (int c : comps) {
      if (c == true_node) continue;
      if (c == false_node) return false_node;
      arcs.push_back(Arc{c, {}});
    }
    if (arcs.empty()) return true_node;
    if (arcs.size() == 1) return arcs[0].child;
    int id = emit(N_AND, arcs);
    // SLOW_DEBUG: strict d-DNNF requires the child subtrees to be pairwise
    // var-disjoint. Catches the violation at emission time (decision context
    // still live) instead of after-the-fact via ddnnf_verify.py.
    SLOW_DEBUG_DO(check_decomposable_at_and(id, arcs));
    return id;
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
    return emit(N_OR, kept);
  }

  // Free (unconstrained) variable: factor of two.
  int mk_free_var(int dimacs_var) {
    std::vector<Arc> arcs;
    arcs.push_back(Arc{true_node, {dimacs_var}});
    arcs.push_back(Arc{true_node, {-dimacs_var}});
    return emit(N_OR, arcs);
  }

  // Wrap a node under a single AND arc carrying `lits` (e.g. level-0 implied lits).
  int wrap_lits(int inner, const std::vector<int>& lits) {
    if (lits.empty()) return inner;
    if (inner == false_node) return false_node;
    std::vector<Arc> arcs;
    arcs.push_back(Arc{inner, lits});
    return emit(N_AND, arcs);
  }

  // ---- per-search bookkeeping (not part of the emitted DAG) ----
  // children[lev][branch] = component nodes AND-ed into that branch
  std::vector<std::array<std::vector<int>, 2>> children;
  // left_lits[lev] = left-branch literals of level lev
  std::vector<std::vector<int>> left_lits;
  // override_node[lev] >= 0: SAT oracle solved this level; node is a witness leaf
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

  size_t num_nodes() const { return n_nodes; }
  size_t num_edges() const { return n_edges; }

  // Stitch the temp files into the final d4 .nnf at `fname`, root declared first.
  // O(1) memory: two passes over the decls temp + a raw copy of the arcs temp.
  void write_d4(const std::string& fname) {
    decl_out.flush();
    arc_out.flush();
    decl_out.close();
    arc_out.close();
    int r = (root < 0) ? false_node : root;

    // Pass 1: find the root's type char.
    char root_char = (r == false_node) ? 'f' : (r == true_node) ? 't' : 'o';
    {
      std::ifstream d(decl_path);
      std::string line;
      while (std::getline(d, line)) {
        if (line.empty()) continue;
        int id = decl_line_id(line);
        if (id == r) { root_char = line[0]; break; }
      }
    }

    std::ofstream o(fname);
    o << root_char << " " << r << " 0\n";
    // Pass 2: copy every declaration except the root's.
    {
      std::ifstream d(decl_path);
      std::string line;
      while (std::getline(d, line)) {
        if (line.empty()) continue;
        if (decl_line_id(line) == r) continue;
        o << line << "\n";
      }
    }
    // Append all arc lines verbatim.
    {
      std::ifstream a(arc_path, std::ios::binary);
      o << a.rdbuf();
    }
    o.close();
    unlink(decl_path.c_str());
    unlink(arc_path.c_str());
  }

private:
  int next_id = 0;
  size_t n_nodes = 0;
  size_t n_edges = 0;
  std::ofstream decl_out;
  std::ofstream arc_out;
  std::string decl_path;
  std::string arc_path;

#ifdef SLOW_DEBUG
  // Per-node set of vars reachable in its subtree (arc lits + child subtree vars).
  // Built incrementally because every arc points to a lower id (children always
  // exist when we add an arc).
  std::vector<std::unordered_set<int>> subtree_vars;
  void record_subtree(int id, const std::vector<Arc>& arcs) {
    if ((int)subtree_vars.size() <= id) subtree_vars.resize(id + 1);
    auto& vs = subtree_vars[id];
    for (const auto& a : arcs) {
      for (int l : a.lits) vs.insert(std::abs(l));
      if (a.child >= 0 && a.child < (int)subtree_vars.size())
        for (int v : subtree_vars[a.child]) vs.insert(v);
    }
  }
  void check_decomposable_at_and(int id, const std::vector<Arc>& arcs) {
    std::vector<std::unordered_set<int>> kid_vars(arcs.size());
    for (size_t i = 0; i < arcs.size(); i++) {
      for (int l : arcs[i].lits) kid_vars[i].insert(std::abs(l));
      if (arcs[i].child >= 0 && arcs[i].child < (int)subtree_vars.size())
        for (int v : subtree_vars[arcs[i].child]) kid_vars[i].insert(v);
    }
    for (size_t i = 0; i < arcs.size(); i++) {
      for (size_t j = i + 1; j < arcs.size(); j++) {
        for (int v : kid_vars[i]) if (kid_vars[j].count(v)) {
          std::cerr << "*** STRICT-DECOMP VIOLATION at AND node " << id
                    << ": var " << v << " in both child " << i
                    << " (node " << arcs[i].child << ") and child " << j
                    << " (node " << arcs[j].child << ").\n    arcs:";
          for (const auto& a : arcs) {
            std::cerr << " (->n" << a.child << " lits:";
            for (int l : a.lits) std::cerr << " " << l;
            std::cerr << ")";
          }
          std::cerr << "\n";
          std::abort();
        }
      }
    }
  }
#endif

  static char type_char(NType t) {
    switch (t) {
      case N_FALSE: return 'f';
      case N_TRUE:  return 't';
      case N_AND:   return 'a';
      default:      return 'o';
    }
  }

  // Parse the node id (second token) out of a declaration line "<c> <id> 0".
  static int decl_line_id(const std::string& line) {
    size_t i = 0;
    while (i < line.size() && line[i] != ' ') i++;     // skip type char
    while (i < line.size() && line[i] == ' ') i++;      // skip spaces
    int id = 0;
    bool any = false;
    while (i < line.size() && line[i] >= '0' && line[i] <= '9') {
      id = id * 10 + (line[i] - '0'); i++; any = true;
    }
    return any ? id : -1;
  }

  void open_temps(const std::string& target_fname) {
    decl_path = mk_temp(target_fname + ".declXXXXXX");
    arc_path = mk_temp(target_fname + ".arcXXXXXX");
    decl_out.open(decl_path, std::ios::trunc);
    arc_out.open(arc_path, std::ios::trunc);
    if (!decl_out.good() || !arc_out.good()) {
      std::cerr << "ERROR: could not open d-DNNF temp files next to "
                << target_fname << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Reserve a unique temp path atomically (mkstemp), close the fd, return name.
  static std::string mk_temp(const std::string& tmpl) {
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    int fd = mkstemp(buf.data());
    if (fd < 0) {
      std::cerr << "ERROR: mkstemp failed for d-DNNF temp file " << tmpl << std::endl;
      exit(EXIT_FAILURE);
    }
    close(fd);
    return std::string(buf.data());
  }

  // Assign an id, stream the declaration + arc lines, and bump the counters.
  int emit(NType t, const std::vector<Arc>& arcs) {
    int id = next_id++;
    decl_out << type_char(t) << " " << id << " 0\n";
    for (const auto& a : arcs) {
      arc_out << id << " " << a.child;
      for (int l : a.lits) arc_out << " " << l;
      arc_out << " 0\n";
      n_edges++;
    }
    n_nodes++;
    // SLOW_DEBUG: maintain per-node subtree var sets so check_decomposable_at_and
    // sees correct child sets even when a child is a leaf (TRUE/FALSE) or OR.
    SLOW_DEBUG_DO(record_subtree(id, arcs));
    return id;
  }
};

}
