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

#include "clauseallocator.hpp"
#include "counter_config.hpp"
#include "statistics.hpp"
#include "comp_types/comp.hpp"
#include "comp_types/comp_archetype.hpp"

#include <map>
#include <gmpxx.h>
#include "containers.hpp"
#include "stack.hpp"

using std::map;

class ClauseAllocator;
class Counter;

// There is exactly ONE of this, inside CompManager, which is inside Solver
class CompAnalyzer {
public:
  CompAnalyzer(
        const LiteralIndexedVector<TriValue> & lit_values,
        const uint32_t& _indep_support_end,
        Counter* _solver);

  Lit const* get_idx_to_cl(uint32_t cl_id) const {
    return idx_to_cl_data.data() + idx_to_cl_map[cl_id];
  }

  const CompArchetype &current_archetype() const { return archetype; }

  void initialize(const LiteralIndexedVector<LitWatchList> & literals,
      const ClauseAllocator* alloc, const vector<ClauseOfs>& long_irred_cls);

  bool isUnseenAndSet(const uint32_t v) const {
    SLOW_DEBUG_DO(assert(v <= max_var));
    return archetype.var_unseen_in_sup_comp(v);
  }

  // manages the literal whenever it occurs in comp analysis
  // returns true iff the underlying variable was unseen before
  void manageSearchOccurrenceOf(const uint32_t v){
    if (archetype.var_unseen_in_sup_comp(v)) {
      comp_vars.push_back(v);
      archetype.set_var_seen(v);
    }
  }

  void setSeenAndStoreInSearchStack(const uint32_t v){
    assert(is_unknown(v));
    comp_vars.push_back(v);
    archetype.set_var_seen(v);
  }

  void setupAnalysisContext(StackLevel &top, const Comp & super_comp){
    archetype.re_initialize(top,super_comp);

    debug_print("Setting VAR/CL_SUP_COMP_UNSEEN in seen[] for vars&cls inside super_comp if unknown");
    for (auto vt = super_comp.vars_begin(); *vt != sentinel; vt++) {
      if (is_unknown(*vt)) {
        archetype.set_var_in_sup_comp_unseen(*vt);
        // TODO what is happening here....
        /* var_freq_scores[*vt] = 0; */
      }
    }

    for (auto it = super_comp.cls_begin(); *it != sentinel; it++)
      archetype.set_clause_in_sup_comp_unseen(*it);
  }

  bool explore_comp(const uint32_t v);

  // explore_comp has been called already
  // which set up search_stack, seen[] etc.
  inline Comp *make_comp_from_archetype(){ return archetype.make_comp(comp_vars.size()); }

  uint32_t get_max_clid() const { return max_clid; }
  uint32_t get_max_var() const { return max_var; }
  CompArchetype& get_archetype() { return archetype; }

private:
  // the id of the last clause
  // note that clause ID is the clause number,
  // different from the offset of the clause in the literal pool
  uint32_t max_clid = 0;
  uint32_t max_var = 0;


  // for every variable e have an array of
  // binarycls 0 ternary cls (consisting of: CLIDX LIT1 LIT2) 0 cls_idxs 0
  vector<uint32_t> unified_var_links_lists_pool;
  vector<uint32_t> variable_link_list_offsets; // offset into unified_var_links_lists_pool
                                                // indexed by variable.

  const CounterConfiguration& conf;
  const LiteralIndexedVector<TriValue> & values;
  const uint32_t& indep_support_end;
  CompArchetype  archetype;
  Counter* solver = nullptr;

  // Quick lookup of cl based on ID
  vector<Lit> idx_to_cl_data; //packed clauses separated by NOT_A_LIT, idx_to_cl_map indexes in
  vector<uint32_t> idx_to_cl_map; //ID goes in, offset of id_to_cl_data comes out. Ends with NOT_A_LIT

  // Used to figure out which vars are in a component
  // used in  record_comp
  // its size is the number of variables in the component
  vector<uint32_t> comp_vars;


  bool is_false(const Lit lit) const { return values[lit] == F_TRI; }
  bool is_true(const Lit lit) const { return values[lit] == T_TRI; }
  bool is_unknown(const Lit lit) const { return values[lit] == X_TRI; }
  bool is_unknown(const uint32_t v) const { return values[Lit(v, true)] == X_TRI; }
  uint32_t const* begin_cls_of_var(const uint32_t v) const {
    assert(v > 0);
    return &unified_var_links_lists_pool[variable_link_list_offsets[v]];
  }
  void bump_var_occs(const uint32_t v);

  // stores all information about the comp of var
  // in variables_seen_, clauses_seen_ and
  // comp_search_stack
  // we have an isolated variable iff
  // after execution comp_search_stack.size()==1
  void record_comp(const uint32_t var);

  void get_cl(vector<uint32_t> &tmp, const Clause& cl, const Lit & omit_lit) {
    tmp.clear();
    for (const auto&l: cl) {
      if (l.var() != omit_lit.var()) tmp.push_back(l.raw());
    }
  }

  // This is called from record_comp, i.e. during figuring out what
  // belongs to a component. It's called on every long clause.
  void search_clause([[maybe_unused]] uint32_t vt, ClauseIndex cl_id, Lit const* pstart_cls){
    const auto it_v_end = comp_vars.end();
    bool all_lits_unkn = true;
    for (auto it_l = pstart_cls; *it_l != SENTINEL_LIT; it_l++) {
      assert(it_l->var() <= max_var);
      if (!archetype.var_nil(it_l->var())) manageSearchOccurrenceOf(it_l->var());
      else {
        assert(!is_unknown(*it_l));
        all_lits_unkn = false;
        if (is_false(*it_l)) continue;

        //accidentally entered a satisfied clause: undo the search process
        while (comp_vars.end() != it_v_end) {
          assert(comp_vars.back() <= max_var);
          archetype.set_var_in_sup_comp_unseen(comp_vars.back()); //unsets it from being seen
          comp_vars.pop_back();
        }
        archetype.clear_cl(cl_id);
        break;
      }
    }

    if (!archetype.clause_nil(cl_id)) archetype.set_clause_seen(cl_id,all_lits_unkn);
  }
};
