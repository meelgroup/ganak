/*
 * alt_comp_analyzer.h
 *
 *  Created on: Mar 5, 2013
 *      Author: mthurley
 */

#pragma once

#include "statistics.h"
#include "comp_types/comp.h"
#include "comp_types/base_packed_comp.h"
#include "comp_types/comp_archetype.h"



#include <vector>
#include <cmath>
#include <gmpxx.h>
#include "containers.h"
#include "stack.h"
#include <set>

using std::set;

// There is exactly ONE of this, inside ComponentManager, which is inside Solver
class ComponentAnalyzer {
public:
	ComponentAnalyzer(
        const LiteralIndexedVector<TriValue> & lit_values,
        const set <uint32_t> & indep_support) :
        lit_values_(lit_values),
        indep_support_(indep_support)
  {}

  uint32_t& scoreOf(VariableIndex v) {
    return var_frequency_scores_[v];
  }

  const uint32_t& scoreOf(VariableIndex v) const {
    return var_frequency_scores_[v];
  }

  ComponentArchetype &current_archetype(){
    return archetype_;
  }

  void initialize(LiteralIndexedVector<LitWatchList> & literals,
      vector<Lit> &lit_pool);

  bool isUnseenAndActive(const VariableIndex v) const {
    assert(v <= max_variable_id_);
    return archetype_.var_unseen_in_sup_comp(v);
  }

  // manages the literal whenever it occurs in comp analysis
  // returns true iff the underlying variable was unseen before
  bool manageSearchOccurrenceOf(const Lit lit){
    if (archetype_.var_unseen_in_sup_comp(lit.var())) {
      search_stack_.push_back(lit.var());
      archetype_.setVar_seen(lit.var());
      return true;
    }
    return false;
  }

  bool manageSearchOccurrenceAndScoreOf(Lit lit){
    var_frequency_scores_[lit.var()]+= isUnknown(lit);
    return manageSearchOccurrenceOf(lit);
  }

  void setSeenAndStoreInSearchStack(const VariableIndex v){
    assert(isUnknown(v));
    search_stack_.push_back(v);
    archetype_.setVar_seen(v);
  }

  void setupAnalysisContext(StackLevel &top, const Component & super_comp){
    archetype_.reInitialize(top,super_comp);

    print_debug("Setting VAR/CL_SUP_COMP_UNSEEN in seen[] for vars&cls inside super_comp if unknown");
    for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++) {
      if (isUnknown(*vt)) {
        archetype_.setVar_in_sup_comp_unseen(*vt);
        var_frequency_scores_[*vt] = 0;
      }
    }

    for (auto itCl = super_comp.clsBegin(); *itCl != clsSENTINEL; itCl++)
      archetype_.setClause_in_sup_comp_unseen(*itCl);
  }

  // returns true, iff the comp found is non-trivial
  bool exploreRemainingCompOf(const VariableIndex v) {
    assert(archetype_.var_unseen_in_sup_comp(v));
    recordComponentOf(v); // finds the comp that "v" is in

    // comp only contains one variable
    if (search_stack_.size() == 1) {
      if (indep_support_.count(v) == 0) {
        archetype_.stack_level().includeSolution(1);
      } else {
        archetype_.stack_level().includeSolution(2);
      }
      archetype_.setVar_in_other_comp(v);
      return false;
    }
    return true;
  }

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

  const LiteralIndexedVector<TriValue> & lit_values_;
  const set <uint32_t> & indep_support_;
  vector<uint32_t> var_frequency_scores_;
  ComponentArchetype  archetype_;
  vector<VariableIndex> search_stack_; // Used to figure out which vars are in a component
                                       // used in  recordComponentOf
                                       // its size is the number of variables in the component

  bool isFalse(const Lit lit) const {
    return lit_values_[lit] == F_TRI;
  }

  bool isTrue(const Lit lit) const {
    return lit_values_[lit] == T_TRI;
  }
  bool isUnknown(const Lit lit) const {
      return lit_values_[lit] == X_TRI;
  }

  bool isUnknown(const VariableIndex v) const {
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

  // Gets a full clause until SENTINEL_LIT, except for the omitLit
  void getClause(
    vector<uint32_t> &tmp,
    vector<Lit>::iterator & it_start_of_cl,
    Lit & omitLit)
  {
    tmp.clear();
    for (auto it_lit = it_start_of_cl; *it_lit != SENTINEL_LIT; it_lit++) {
      if (it_lit->var() != omitLit.var())
        tmp.push_back(it_lit->raw());
    }
  }

  // This is called from recordComponentOf, i.e. during figuring out what
  // belongs to a component. It's called on every long clause.
  void searchClause(VariableIndex vt, ClauseIndex clID, Lit const* pstart_cls){
    const auto itVEnd = search_stack_.end();
    bool all_lits_active = true;
    for (auto itL = pstart_cls; *itL != SENTINEL_LIT; itL++) {
      assert(itL->var() <= max_variable_id_);
      if(!archetype_.var_nil(itL->var()))
        manageSearchOccurrenceAndScoreOf(*itL); // sets var to be seen
      else {
        assert(!isUnknown(*itL));
        all_lits_active = false;
        if (isFalse(*itL)) continue;

        //accidentally entered a satisfied clause: undo the search process
        while (search_stack_.end() != itVEnd) {
          assert(search_stack_.back() <= max_variable_id_);
          archetype_.setVar_in_sup_comp_unseen(search_stack_.back()); //unsets it from being seen
          search_stack_.pop_back();
        }
        archetype_.setClause_nil(clID);
        while(*itL != SENTINEL_LIT)
          if(isUnknown(*(--itL))) var_frequency_scores_[itL->var()]--;
        break;
      }
    }

    if (!archetype_.clause_nil(clID)){
      var_frequency_scores_[vt]++;
      archetype_.setClause_seen(clID,all_lits_active);
    }
  }

};
