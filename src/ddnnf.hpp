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
#include "bdnnf.hpp"
#include "ddnnf_check.hpp"
#include "GitSHA1.hpp"

#include <algorithm>
#include <memory>
#include <vector>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <iostream>
#include <unordered_set>

namespace GanakInt {

// Builds a (Decision-)d-DNNF circuit from Ganak's search trace and streams it to
// disk as BDNNF (spec-bdnnf.md), Huang & Darwiche 2005 style:
//   decision          -> OR node (arcs carry decision + implied lits)
//   AND-decomposition -> AND node
//   free var          -> free var of the enclosing AND (factor 2)
//   SAT leaf -> TRUE;  conflict -> FALSE;  cache hit -> shared edge
// Every node is complete at creation and its children already exist, so it is
// written out immediately. The structural count reproduces count().
//
// Free vars are returned by mk_free_var() as the pseudo-id -v and only become
// part of a node in mk_and(); they never reach the cache or an OR arc.
class DDNNFCircuit {
public:
  struct Arc {
    int child;
    std::vector<int> lits; // DIMACS literals true on this arc
  };

  DDNNFCircuit(const std::string& fname, uint32_t _nvars, std::vector<int> _lit_map)
      : nvars(_nvars), lit_map(std::move(_lit_map)) {
    for (int l : lit_map) nvars = std::max<uint32_t>(nvars, (uint32_t)std::abs(l));
    guard([&] {
      w.open(fname, nvars, {{"producer", std::string("ganak ") +
          GANAK::get_version_tag() + " " + GANAK::get_version_sha1()}});
    });
  }

  static constexpr int false_node = 0;
  static constexpr int true_node = 1;
  int root = -1;
  uint32_t nvars;

  // AND of children: drop TRUE; any FALSE -> FALSE; empty -> TRUE;
  // single real child -> that child.
  int mk_and(const std::vector<int>& comps) {
    kids.clear();
    frees.clear();
    for (int c : comps) {
      if (c == true_node) continue;
      if (c == false_node) return false_node;
      if (c < 0) frees.push_back(-c);
      else kids.push_back(c);
    }
    if (kids.empty() && frees.empty()) return true_node;
    if (kids.size() == 1 && frees.empty()) return kids[0];
    const int id = emit_and(kids, no_lits, frees);
    return id;
  }

  // OR of (literal-set -> child) arcs. Drop arcs to FALSE; empty -> FALSE.
  int mk_or(std::vector<Arc>&& arcs) {
    size_t kept = 0;
    for (const auto& a : arcs) {
      assert(a.child >= 0);
      if (a.child != false_node) kept++;
    }
    if (kept == 0) return false_node;
    int id = -1;
    std::vector<Arc> written;
    guard([&] {
      w.begin_or();
      for (const auto& a : arcs) if (a.child != false_node) {
        w.or_arc(mapped(a.lits), (uint64_t)a.child);
        if (chk) written.push_back(Arc{a.child, mapped(a.lits)});
      }
      id = check_id(w.end_or());
    });
    if (chk) chk->on_or((uint64_t)id, written);
    return id;
  }

  static int mk_free_var(int dimacs_var) { return -dimacs_var; }

  // inner AND lits (e.g. level-0 implied lits, or a SAT witness on TRUE).
  int wrap_lits(int inner, const std::vector<int>& lits) {
    if (lits.empty()) return inner;
    if (inner == false_node) return false_node;
    kids.clear();
    if (inner != true_node) kids.push_back(inner);
    const int id = emit_and(kids, lits, no_lits);
    return id;
  }

  // ---- per-search bookkeeping (not part of the emitted DAG) ----
  // children[lev][branch] = component nodes AND-ed into that branch
  // TODO this should be somewhat more persistent for mem alloc efficiency
  std::vector<std::array<std::vector<int>, 2>> children;
  // left_lits[lev] = left-branch literals of level lev
  std::vector<std::vector<int>> left_lits;
  // override_node[lev] >= 0: SAT oracle solved this level; node is a witness leaf
  std::vector<int> override_node;
  // split_lits[lev][branch]: the branch's lits when it was (last) split into
  // comps, i.e. the assignment its count is based on; valid if has_split.
  std::vector<std::array<std::vector<int>, 2>> split_lits;
  std::vector<std::array<char, 2>> has_split;


  // TODO this should be insert, for mem perf
  void ensure_level(int lev) {
    if ((int)children.size() <= lev) {
      children.resize(lev + 1);
      left_lits.resize(lev + 1);
      override_node.resize(lev + 1, -1);
      split_lits.resize(lev + 1);
      has_split.resize(lev + 1, {0, 0});
    }
  }
  void on_new_level(int lev) {
    ensure_level(lev);
    children[lev][0].clear();
    children[lev][1].clear();
    left_lits[lev].clear();
    override_node[lev] = -1;
    has_split[lev] = {0, 0};
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

  [[nodiscard]] uint64_t num_nodes() const { return w.num_nodes(); }

  // --compilecheck: check every node as it is written (see ddnnf_check.hpp).
  void enable_check(bool expect_smooth) { chk = std::make_unique<DDNNFChecker>(nvars, expect_smooth); }
  [[nodiscard]] DDNNFChecker* checker() const { return chk.get(); }
  [[nodiscard]] uint64_t num_edges() const { return w.num_edges(); }

  void finish() {
    const int r = (root < 0) ? false_node : root;
    // Decision ORs split on x / -x. Arcs carry only their own component's lits of
    // their own level, and --compile does no chronological backtracking, so no
    // lit is recorded twice: decomposable (checked by --compilecheck and the fuzzers).
    guard([&] { w.finish((uint64_t)r, bdnnf::P_DETERMINISTIC | bdnnf::P_DECOMPOSABLE); });
  }

private:
  bdnnf::Writer w;
  std::unique_ptr<DDNNFChecker> chk;
  std::vector<int> lit_map;
  std::vector<int> scratch_lits;
  std::vector<int> kids;
  std::vector<int> frees;
  const std::vector<int> no_lits;

  template <class F>
  static void guard(F&& f) {
    try {
      f();
    } catch (const bdnnf::Error& e) {
      std::cerr << "ERROR: writing d-DNNF: " << e.what() << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  static int check_id(uint64_t id) {
    if (id > (uint64_t)INT32_MAX) {
      std::cerr << "ERROR: d-DNNF has more than 2^31 nodes" << std::endl;
      exit(EXIT_FAILURE);
    }
    return (int)id;
  }

  // Internal DIMACS literals -> the input's numbering.
  const std::vector<int>& mapped(const std::vector<int>& lits) {
    if (lit_map.empty()) return lits;
    scratch_lits.clear();
    for (int l : lits) {
      const int m = lit_map[std::abs(l) - 1];
      scratch_lits.push_back(l > 0 ? m : -m);
    }
    return scratch_lits;
  }

  int emit_and(const std::vector<int>& k, const std::vector<int>& units,
               const std::vector<int>& fr) {
    int id = -1;
    std::vector<int> mu = mapped(units);
    std::vector<int> mf = mapped(fr);
    for (int& v : mf) v = std::abs(v);
    guard([&] { id = check_id(w.add_and(k, mu, mf)); });
    if (chk) chk->on_and((uint64_t)id, k, mu, mf);
    return id;
  }


};

}
