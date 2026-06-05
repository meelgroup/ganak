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

#include "compiler.hpp"
#include "counter.hpp"
#include "ddnnf.hpp"

#include <vector>

namespace GanakInt {

// Trivial no-op compiler: every event is ignored. This is what the engine holds
// when `--compile` is off, which is why none of the hook call sites need guards.
struct NullCompiler final : public Compiler {};

// Real compiler: bridges the engine's search events to a streamed DDNNFCircuit.
// It is a friend of Counter (see counter.hpp) so it can read the search state
// (trail, decisions, components) it needs to label arcs and decompose nodes.
class DDNNFCompiler final : public Compiler {
public:
  DDNNFCompiler(Counter& _counter, const std::string& _fname)
      : counter(_counter), fname(_fname) {
    circuit = std::make_unique<DDNNFCircuit>(fname);
    circuit->nvars = counter.nVars();
  }

  bool active() const override { return true; }
  bool allows_cms_shortcut() const override { return false; }

  void new_decision_level() override { circuit->on_new_level(counter.dec_level()); }

  void save_left_lits() override {
    // d-DNNF: capture left-branch literals before the trail is undone.
    const int lev = counter.dec_level();
    circuit->ensure_level(lev);
    circuit->left_lits[lev] = cur_level_lits();
  }

  void build_level_node() override {
    // Build this level's OR node before the trail is undone (to read the right
    // branch's literals). The result is held until record_comp_node() /
    // attach_to_parent() consume it.
    const int lev = counter.dec_level();
    const auto right_lits = cur_level_lits();
    auto& top = counter.decisions.top();
    circuit->ensure_level(lev);
    // SAT-oracle leaf: level solved by the SAT solver, which recorded a witness for
    // the synthesized vars. Use that leaf instead of an OR over (stale) children.
    if (int ov = circuit->take_override(lev); ov >= 0) { compiled_node = ov; return; }
    int left  = top.is_zero(0) ? circuit->false_node : circuit->mk_and(circuit->children[lev][0]);
    int right = top.is_zero(1) ? circuit->false_node : circuit->mk_and(circuit->children[lev][1]);
    std::vector<DDNNFCircuit::Arc> arcs;
    if (left  != circuit->false_node) arcs.push_back(DDNNFCircuit::Arc{left,  circuit->left_lits[lev]});
    if (right != circuit->false_node) arcs.push_back(DDNNFCircuit::Arc{right, right_lits});
    VERBOSE_DEBUG_DO({
      std::cerr << "[OR-build] lev=" << lev << " super_comp_id=" << top.super_comp() << " super_vars={";
      const auto& sc = counter.comp_manager->get_super_comp(top);
      all_vars_in_comp(sc, vp) std::cerr << " " << *vp;
      std::cerr << " } left_lits={";
      for (int l : circuit->left_lits[lev]) std::cerr << " " << l;
      std::cerr << " } right_lits={";
      for (int l : right_lits) std::cerr << " " << l;
      std::cerr << " } left_kid=" << left << " right_kid=" << right << "\n";
    });
    compiled_node = circuit->mk_or(std::move(arcs));
  }

  void record_comp_node() override {
    counter.comp_manager->set_comp_node(counter.decisions.top().super_comp(), compiled_node);
  }

  void attach_to_parent() override {
    circuit->add_child(counter.dec_level() - 1,
        (counter.decisions.end() - 2)->is_right_branch(), compiled_node);
  }

  void branch_reset() override {
    // Mirror zero_out_branch_sol(): drop the current branch's children, as Ganak
    // discards its count to re-explore.
    const int lev = counter.dec_level();
    circuit->ensure_level(lev);
    circuit->children[lev][counter.decisions.top().is_right_branch()].clear();
  }

  void free_var(uint32_t v) override {
    // From the analyzer on a free var (factor two = OR(v, -v)).
    circuit->add_child(counter.dec_level(), counter.decisions.top().is_right_branch(),
        circuit->mk_free_var((int)v));
  }

  void cache_hit(int node) override {
    // A cached component must have its sub-DAG recorded (set_comp_node in
    // backtrack); a missing node would silently undercount, so fail loudly.
    release_assert(node >= 0); // d-DNNF cache hit with no recorded compile node
    circuit->add_child(counter.dec_level(), counter.decisions.top().is_right_branch(), node);
  }

  void sat_witness_capture(int sat_start_dec_level) override {
    // Record the SAT oracle's witness for this component's synthesized vars
    // (>= opt_indep_support_end) before backtracking.
    sat_witness.clear();
    all_vars_in_comp(counter.comp_manager->get_super_comp(counter.decisions.at(sat_start_dec_level)), it) {
      uint32_t const v = *it;
      if (v >= counter.opt_indep_support_end && counter.val(v) != X_TRI)
        sat_witness.push_back(Lit(v, counter.val(v) == T_TRI).to_visual_int());
    }
  }

  void sat_witness_apply() override {
    // This SAT level's node is the witness leaf.
    circuit->set_override(counter.dec_level(), circuit->wrap_lits(circuit->true_node, sat_witness));
  }

  void finalize_root() override {
    circuit->ensure_level(0);
    // Level 0 starts on the right branch (init_decision_stack), so top-level
    // components accumulate there.
    const bool b = counter.decisions.top().is_right_branch();
    int inner = circuit->mk_and(circuit->children[0][b]);
    circuit->root = circuit->wrap_lits(inner, level0_lits());
  }

  void finalize_and_write(bool unsat) override {
    if (unsat) {
      circuit->root = circuit->false_node;
      circuit->write_d4(fname);
      return;
    }
    if (circuit->root < 0) circuit->root = circuit->false_node;
    circuit->nvars = counter.nVars();
    circuit->write_d4(fname);
    const auto& conf = counter.get_conf();
    verb_print(0, "[compile] d-DNNF written to " << fname
        << " nodes: " << circuit->num_nodes() << " edges: " << circuit->num_edges());
  }

private:
  // Current top decision level's trail lits.
  std::vector<int> cur_level_lits() const {
    std::vector<int> r;
    for (auto it = counter.top_declevel_trail_begin(); it != counter.trail.end(); ++it)
      r.push_back(it->to_visual_int());
    return r;
  }
  // Literals forced at decision level 0.
  std::vector<int> level0_lits() const {
    std::vector<int> r;
    for (const auto& l : counter.trail)
      if (counter.var(l).decision_level == 0) r.push_back(l.to_visual_int());
    return r;
  }

  Counter& counter;
  std::string fname;
  std::unique_ptr<DDNNFCircuit> circuit;
  int compiled_node = -1;       // node built by build_level_node(), consumed downstream
  std::vector<int> sat_witness; // DIMACS lits of synthesized (Y) vars from the SAT oracle
};

std::unique_ptr<Compiler> make_null_compiler() {
  return std::make_unique<NullCompiler>();
}

std::unique_ptr<Compiler> make_ddnnf_compiler(Counter& counter, const std::string& fname) {
  return std::make_unique<DDNNFCompiler>(counter, fname);
}

}
