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
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <unistd.h>

namespace GanakInt {

// Builds a (Decision-)d-DNNF circuit out of Ganak's search trace and STREAMS it
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
// the value Ganak's internal count() returns.
//
// STREAMING. The whole DAG is never held in memory. Two invariants make this
// possible:
//   (1) a node is COMPLETE at creation -- every builder below sets all of a
//       node's arcs at the moment its id is assigned, and they are never
//       mutated afterwards; and
//   (2) every arc points to a LOWER id, because a node is only ever created
//       after all of its children already exist (cache hits, free vars,
//       SAT/override leaves and wrapped inner nodes are all pre-existing).
// So each node can be emitted the instant it is built, and a child is always
// already on disk before any parent references it. We stream declaration lines
// to one temp file and arc lines to another, then finalize() stitches them into
// the target file with the root declared first (the d4 convention: the
// first-declared node is the root).
//
// The emitted file may contain UNREACHABLE ("dead") nodes: mk_and short-circuits
// to FALSE when any child is FALSE, orphaning already-emitted non-false
// siblings. They are harmless to a root-rooted traversal (count/synthesis), and
// the separate `ddnnf-cleanup` tool removes them + renumbers root=1 contiguous
// for consumers that require a strict d4 file.
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
    std::vector<Arc> arcs;
    arcs.reserve(kids.size());
    for (int c : kids) arcs.push_back(Arc{c, {}});
    return emit(N_AND, arcs);
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

  // A free independent variable (unconstrained): contributes a factor of two.
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

  size_t num_nodes() const { return n_nodes; }
  size_t num_edges() const { return n_edges; }

  // Stitch the streamed temp files into the final d4 .nnf at `fname`, with the
  // root declared first. O(1) memory: two sequential passes over the (decls-only)
  // temp file plus a raw copy of the (arcs-only) temp file.
  void write_d4(const std::string& fname) {
    decl_out.flush();
    arc_out.flush();
    decl_out.close();
    arc_out.close();
    int r = (root < 0) ? false_node : root;

    // Pass 1: find the root's node type char (so we can declare it first).
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
    // Pass 2: copy every declaration except the root's (already written above).
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

  // Reserve a unique temp path atomically (O_CREAT|O_EXCL via mkstemp), close
  // the fd, and return the name for an ofstream to truncate+reuse.
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
    return id;
  }
};

}
