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

// bdnnf cleanup: canonicalise a BDNNF circuit from `ganak --compile`.
//
//   1. FLATTEN nested ANDs and ELIDE unary ANDs into their parent arcs.
//   2. DROP dead (unreachable) nodes: only the part reachable from the root is
//      written, root last.
//
// Output is flagged CLEAN, plus whichever of DECOMPOSABLE, DETERMINISTIC,
// SMOOTH it was checked to satisfy.
//
// Usage:  bdnnf cleanup <in.bdnnf> [out.bdnnf]      ("-" or omitted out => stdout)

#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "ddnnf_io.hpp"

namespace {

using ddnnf_io::Arc;

[[noreturn]] void die(const std::string& msg) {
  std::cerr << "bdnnf cleanup: " << msg << std::endl;
  exit(EXIT_FAILURE);
}

// FLATTEN NESTED ANDs (associativity): an AND reached from another AND through a
// *literal-free* arc is redundant -- splice its arcs straight into the parent.
// The dominant source is wrap_lits() component witnesses: a projected / SAT-
// solved sub-component is emitted as a single-arc AND (a TRUE child carrying the
// witness lits), then mk_and() conjoins it with its sibling components, giving
// AND( AND(TRUE; witness), sibling... ). Inlining rewrites the parent's arc to
// that wrapper into a direct arc carrying the witness lits, dropping a whole AND
// layer (the "AND -> AND -> literals" stack seen in bdnnf dot output).
//
// mk_and() always emits literal-free AND arcs, so this catches every such stack.
// Count, model set, per-path literals and strict decomposability are preserved:
// the inlined children keep their own arcs/lits and (being one decomposable AND
// child of the parent) were already var-disjoint from the parent's other
// children. A now-unreachable inner node is dropped by the BFS renumber.
//
// Literal-carrying AND arcs (e.g. the level-0 root wrap_lits) are left alone --
// moving their lits onto one spliced arc could place that var in a child that is
// no longer var-disjoint from its new siblings, breaking decomposability.
int flatten_ands(std::vector<char>& type, std::vector<std::vector<Arc>>& arcs) {
  int inlined = 0;
  bool changed = true;
  while (changed) {
    changed = false;
    for (int id = 0; id < (int)arcs.size(); id++) {
      if (id >= (int)type.size() || type[id] != 'a') continue;
      std::vector<Arc> out;
      out.reserve(arcs[id].size());
      for (auto& a : arcs[id]) {
        const bool inlineable =
            a.lits.empty() && a.child >= 0 && a.child != id &&
            a.child < (int)type.size() && type[a.child] == 'a' &&
            !arcs[a.child].empty();
        if (!inlineable) { out.push_back(std::move(a)); continue; }
        // Copy (not move): the inner AND may be DAG-shared by other parents.
        for (const auto& g : arcs[a.child]) out.push_back(g);
        inlined++;
        changed = true;
      }
      arcs[id] = std::move(out);
    }
  }
  return inlined;
}

// ELIDE UNARY ANDs: an AND with a single arc conjoins nothing -- it is just its
// child, plus the literals forced on that arc. wrap_lits() emits one for every
// projected / SAT-witnessed component (a lone witness arc into TRUE); when its
// parent is a decision OR (not another AND) flatten_ands above leaves it, so it
// survives as a degenerate one-input AND gate. Fold it into each parent arc:
// parent arc (A, Lp) becomes (G, Lp u La), where A's only arc is (G, La). The
// node, once unreferenced, is dropped by the BFS renumber. Count, model set,
// per-path literals and decomposability are preserved -- the lits stay together
// on the single arc that reaches G. The root is never elided: it has no parent
// to fold into and legitimately carries the top-level (level-0) literals.
int elide_unary_ands(int root, std::vector<char>& type,
                     std::vector<std::vector<Arc>>& arcs) {
  int elided = 0;
  bool changed = true;
  while (changed) {
    changed = false;
    for (int p = 0; p < (int)arcs.size(); p++) {
      if (p >= (int)type.size() || type[p] == 0) continue;
      for (auto& a : arcs[p]) {
        const int c = a.child;
        if (c == root || c == p || c < 0 || c >= (int)type.size()) continue;
        if (type[c] != 'a' || arcs[c].size() != 1) continue;
        // Fold the unary AND `c` into this arc (copy lits first: arcs[c] is
        // const here, and `c` may be DAG-shared by other parents).
        const Arc inner = arcs[c][0];
        a.child = inner.child;
        a.lits.insert(a.lits.end(), inner.lits.begin(), inner.lits.end());
        elided++;
        changed = true;
      }
    }
  }
  return elided;
}

// MERGE an AND's arcs into TRUE into one: their literals are all units of the
// AND (BDNNF stores them as one list), so an AND of cubes becomes one cube,
// which elide_unary_ands can then fold into its parent arc.
int merge_true_arcs(std::vector<char>& type, std::vector<std::vector<Arc>>& arcs) {
  int merged = 0;
  for (int id = 0; id < (int)arcs.size(); id++) {
    if (type[id] != 'a') continue;
    int first = -1;
    std::vector<Arc> out;
    for (auto& a : arcs[id]) {
      if (type[a.child] != 't') { out.push_back(std::move(a)); continue; }
      if (first < 0) { first = (int)out.size(); out.push_back(std::move(a)); continue; }
      out[first].lits.insert(out[first].lits.end(), a.lits.begin(), a.lits.end());
      merged++;
    }
    arcs[id] = std::move(out);
  }
  return merged;
}

} // namespace

int bdnnf_cleanup(int argc, char** argv) {
  int argi = 1;
  // Accepted for old scripts: there is no repair pass any more.
  if (argi < argc && std::string(argv[argi]) == "--no-strict-decomp") argi++;
  if (argc - argi < 1 || argc - argi > 2) {
    std::cerr << "Usage: bdnnf cleanup <in.bdnnf> [out.bdnnf]" << std::endl;
    return EXIT_FAILURE;
  }
  const std::string in_path = argv[argi];
  const std::string out_path = (argc - argi == 2) ? argv[argi + 1] : "-";

  ddnnf_io::Circuit c;
  std::string err;
  if (!ddnnf_io::load(in_path, c, err)) die(err);

  int flat = 0, unary = 0, f, u;
  do {
    f = flatten_ands(c.type, c.arcs);
    f += merge_true_arcs(c.type, c.arcs);
    u = elide_unary_ands(c.root, c.type, c.arcs);
    flat += f;
    unary += u;
  } while (f || u);
  std::cerr << "bdnnf cleanup: flattened-nested-and arcs=" << flat
            << " elided-unary-and=" << unary << std::endl;

  uint64_t props = bdnnf::P_CLEAN;
  if (ddnnf_io::check_decomposable(c).first) props |= bdnnf::P_DECOMPOSABLE;
  else std::cerr << "bdnnf cleanup: WARN: input is not decomposable" << std::endl;
  if ((c.props & bdnnf::P_DETERMINISTIC) && ddnnf_io::check_deterministic(c).first)
    props |= bdnnf::P_DETERMINISTIC;
  if (ddnnf_io::check_smooth(c).first) props |= bdnnf::P_SMOOTH;

  if (!ddnnf_io::save(out_path, c, props, "bdnnf cleanup", err, c.footer_kv)) die(err);
  std::cerr << "bdnnf cleanup: input nodes=" << c.nfile << " props=" << props << std::endl;
  return EXIT_SUCCESS;
}
