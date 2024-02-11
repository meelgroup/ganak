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
#include "comp_types/base_packed_comp.hpp"
#include "comp_types/comp_archetype.hpp"

#include <vector>
#include <cmath>
#include <map>
#include <gmpxx.h>
#include "containers.hpp"
#include "stack.hpp"

using std::map;

class ClauseAllocator;
class Counter;

// There is exactly ONE of this, inside ComponentManager, which is inside Solver
class ComponentAnalyzer {
public:
  ComponentAnalyzer(
        const LiteralIndexedVector<TriValue> & lit_values,
        const uint32_t& _indep_support_end,
        Counter* _solver);

  const map<uint32_t, vector<Lit>>& get_idx_to_cl() const {
    return idx_to_cl;
  }

  double score_of(VariableIndex v) const { return var_frequency_scores_[v]/act_inc; }
  void un_bump_score(VariableIndex v) {
    var_frequency_scores_[v] -= act_inc;
  }
  inline void bump_score(VariableIndex v) {
    var_frequency_scores_[v] += act_inc;
    if (var_frequency_scores_[v] > 1e100) {
      for(auto& f: var_frequency_scores_) f *= 1e-90;
      act_inc *= 1e-90;
    }
    if ((conf.decide & 2) == 0) act_inc *= 1.0/0.98;
  }
  const ComponentArchetype &current_archetype() const { return archetype_; }

  void initialize(LiteralIndexedVector<LitWatchList> & literals,
      const ClauseAllocator* alloc, const vector<ClauseOfs>& long_irred_cls);

  bool isUnseenAndSet(const VariableIndex v) const {
    SLOW_DEBUG_DO(assert(v <= max_variable_id_));
    return archetype_.var_unseen_in_sup_comp(v);
  }

  // manages the literal whenever it occurs in comp analysis
  // returns true iff the underlying variable was unseen before
  bool manageSearchOccurrenceOf(const uint32_t v){
    if (archetype_.var_unseen_in_sup_comp(v)) {
      search_stack_.push_back(v);
      archetype_.setVar_seen(v);
      return true;
    }
    return false;
  }

  bool manageSearchOccurrenceAndScoreOf(Lit lit){
  if (is_unknown(lit)) bump_score(lit.var());
    return manageSearchOccurrenceOf(lit.var());
  }

  void setSeenAndStoreInSearchStack(const VariableIndex v){
    assert(is_unknown(v));
    search_stack_.push_back(v);
    archetype_.setVar_seen(v);
  }

  void setupAnalysisContext(StackLevel &top, const Component & super_comp){
    archetype_.reInitialize(top,super_comp);

    debug_print("Setting VAR/CL_SUP_COMP_UNSEEN in seen[] for vars&cls inside super_comp if unknown");
    for (auto vt = super_comp.varsBegin(); *vt != sentinel; vt++) {
      if (is_unknown(*vt)) {
        archetype_.setVar_in_sup_comp_unseen(*vt);
        // TODO what is happening here....
        /* var_frequency_scores_[*vt] = 0; */
      }
    }

    for (auto itCl = super_comp.clsBegin(); *itCl != sentinel; itCl++)
      archetype_.setClause_in_sup_comp_unseen(*itCl);
  }

  bool exploreRemainingCompOf(const VariableIndex v);

  // exploreRemainingCompOf has been called already
  // which set up search_stack, seen[] etc.
  inline Component *makeComponentFromArcheType(){
    return archetype_.makeComponentFromState(search_stack_.size());
  }

  uint32_t max_clause_id() const {
     return max_clause_id_;
  }
  uint32_t max_variable_id() const {
    return max_variable_id_;
  }

  ComponentArchetype &getArchetype() {
    return archetype_;
  }

private:
  // the id of the last clause
  // note that clause ID is the clause number,
  // different from the offset of the clause in the literal pool
  uint32_t max_clause_id_ = 0;
  uint32_t max_variable_id_ = 0;


  // for every variable e have an array of
  // binarycls 0 ternary cls (consisting of: CLIDX LIT1 LIT2) 0 cls_idxs 0
  vector<uint32_t> unified_variable_links_lists_pool_;
  vector<uint32_t> variable_link_list_offsets_; // offset into unified_variable_links_lists_pool_
                                                // indexed by variable.

  const CounterConfiguration& conf;
  const LiteralIndexedVector<TriValue> & lit_values_;
  const uint32_t& indep_support_end;
  vector<double> var_frequency_scores_;
  double act_inc = 1.0;
  ComponentArchetype  archetype_;
  Counter* solver = nullptr;
  map<uint32_t, vector<Lit>> idx_to_cl;
  vector<VariableIndex> search_stack_; // Used to figure out which vars are in a component
                                       // used in  recordComponentOf
                                       // its size is the number of variables in the component

  bool is_false(const Lit lit) const {
    return lit_values_[lit] == F_TRI;
  }

  bool is_true(const Lit lit) const {
    return lit_values_[lit] == T_TRI;
  }
  bool is_unknown(const Lit lit) const {
      return lit_values_[lit] == X_TRI;
  }

  bool is_unknown(const VariableIndex v) const {
    return lit_values_[Lit(v, true)] == X_TRI;
  }

  uint32_t const* begin_cls_of_var(const VariableIndex v) const {
    assert(v > 0);
    return &unified_variable_links_lists_pool_[variable_link_list_offsets_[v]];
  }

  // stores all information about the comp of var
  // in variables_seen_, clauses_seen_ and
  // comp_search_stack
  // we have an isolated variable iff
  // after execution comp_search_stack.size()==1
  void recordComponentOf(const VariableIndex var);

  void getClause(vector<uint32_t> &tmp, const Clause& cl, const Lit & omitLit) {
    tmp.clear();
    for (const auto&l: cl) {
      if (l.var() != omitLit.var()) tmp.push_back(l.raw());
    }
  }

  // This is called from recordComponentOf, i.e. during figuring out what
  // belongs to a component. It's called on every long clause.
  void searchClause(VariableIndex vt, ClauseIndex clID, Lit const* pstart_cls){
    const auto itVEnd = search_stack_.end();
    bool all_lits_set = true;
    for (auto itL = pstart_cls; *itL != SENTINEL_LIT; itL++) {
      assert(itL->var() <= max_variable_id_);
      if(!archetype_.var_nil(itL->var()))
        manageSearchOccurrenceAndScoreOf(*itL); // sets var to be seen
      else {
        assert(!is_unknown(*itL));
        all_lits_set = false;
        if (is_false(*itL)) continue;

        //accidentally entered a satisfied clause: undo the search process
        while (search_stack_.end() != itVEnd) {
          assert(search_stack_.back() <= max_variable_id_);
          archetype_.setVar_in_sup_comp_unseen(search_stack_.back()); //unsets it from being seen
          search_stack_.pop_back();
        }
        archetype_.setClause_nil(clID);
        while(*itL != SENTINEL_LIT)
          if(is_unknown(*(--itL))) un_bump_score(itL->var());
        break;
      }
    }

    if (!archetype_.clause_nil(clID)){
      bump_score(vt);
      archetype_.setClause_seen(clID,all_lits_set);
    }
    if ((conf.decide & 2)) act_inc *= 1.0/0.98;
  }
};
