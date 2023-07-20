/*
 * comp_archetype.h
 *
 *  Created on: Feb 9, 2013
 *      Author: mthurley
 */

#ifndef COMPONENT_ARCHETYPE_H_
#define COMPONENT_ARCHETYPE_H_

#include <cstring>
#include <algorithm>
#include <iostream>

#include "../primitive_types.hpp"
#include "common.hpp"
#include "comp.hpp"
#include "cacheable_comp.hpp"

// State values for variables found during comp
// analysis (CA)
typedef uint8_t CA_SearchState;
#define   CA_VAR_IN_SUP_COMP_UNSEEN  1
#define   CA_VAR_SEEN 2
#define   CA_VAR_IN_OTHER_COMP  4

// 1 + 2 + 4 == 7
#define   CA_VAR_MASK  7

#define   CA_CL_IN_SUP_COMP_UNSEEN  8
#define   CA_CL_SEEN 16
#define   CA_CL_IN_OTHER_COMP  32
#define   CA_CL_ALL_LITS_SET  64

// 64+32+16+8 == 120
#define   CA_CL_MASK  120

class StackLevel;

// There is exactly ONE of this. Inside ComponentAnalyzer, which is inside ComponentManager, which is inside Solver
class ComponentArchetype {
public:
  ComponentArchetype() { }
  ~ComponentArchetype() { delete[] seen_; }
  ComponentArchetype(StackLevel &stack_level, const Component &super_comp) :
      p_super_comp_(&super_comp), p_stack_level_(&stack_level) {
  }

  void reInitialize(StackLevel &stack_level, const Component &super_comp) {
    debug_print("Reinitializing seen[] to all-zero in ComponentArchetype");
    p_super_comp_ = &super_comp;
    p_stack_level_ = &stack_level;
    clearArrays();
    num_bin_cls = 0;
    current_comp_for_caching_.reserveSpace(super_comp.nVars(),super_comp.numLongClauses());
  }

  const Component &super_comp() {
    return *p_super_comp_;
  }

  StackLevel & stack_level() {
    return *p_stack_level_;
  }

  void setVar_in_sup_comp_unseen(const VariableIndex v) {
    seen_[v] = CA_VAR_IN_SUP_COMP_UNSEEN | (seen_[v] & CA_CL_MASK);
  }

  void setClause_in_sup_comp_unseen(const ClauseIndex cl) {
    seen_[cl] = CA_CL_IN_SUP_COMP_UNSEEN | (seen_[cl] & CA_VAR_MASK);
  }

  void setVar_nil(const VariableIndex v) {
    seen_[v] &= CA_CL_MASK;
  }

  void setClause_nil(const ClauseIndex cl) {
    seen_[cl] &= CA_VAR_MASK;
  }

  void setVar_seen(const VariableIndex v) {
    seen_[v] = CA_VAR_SEEN | (seen_[v] & CA_CL_MASK);
  }

  void setClause_seen(const ClauseIndex cl) {
    setClause_nil(cl);
    seen_[cl] = CA_CL_SEEN | (seen_[cl] & CA_VAR_MASK);
  }

  void setClause_seen(const ClauseIndex cl, const bool all_lits_act) {
      setClause_nil(cl);
      seen_[cl] = CA_CL_SEEN | (all_lits_act?CA_CL_ALL_LITS_SET:0) | (seen_[cl] & CA_VAR_MASK);
    }

  void setVar_in_other_comp(const VariableIndex v) {
    seen_[v] = CA_VAR_IN_OTHER_COMP | (seen_[v] & CA_CL_MASK);
  }

  void setClause_in_other_comp(const ClauseIndex cl) {
    seen_[cl] = CA_CL_IN_OTHER_COMP | (seen_[cl] & CA_VAR_MASK);
  }

  bool var_seen(const VariableIndex v) const {
    return seen_[v] & CA_VAR_SEEN;
  }

  bool clause_seen(const ClauseIndex cl) const {
    return seen_[cl] & CA_CL_SEEN;
  }

  bool clause_all_lits_set(const ClauseIndex cl) const {
    return seen_[cl] & CA_CL_ALL_LITS_SET;
  }

  bool var_nil(const VariableIndex v) const {
    return (seen_[v] & CA_VAR_MASK) == 0;
  }

  bool clause_nil(const ClauseIndex cl) const {
    return (seen_[cl] & CA_CL_MASK) == 0;
  }

  bool var_unseen_in_sup_comp(VariableIndex v) const {
    return seen_[v] & CA_VAR_IN_SUP_COMP_UNSEEN;
  }

  bool clause_unseen_in_sup_comp(ClauseIndex cl) const {
    return seen_[cl] & CA_CL_IN_SUP_COMP_UNSEEN;
  }

  bool var_seen_in_peer_comp(VariableIndex v) const {
    return seen_[v] & CA_VAR_IN_OTHER_COMP;
  }

  bool clause_seen_in_peer_comp(ClauseIndex cl) const {
    return seen_[cl] & CA_CL_IN_OTHER_COMP;
  }

  void initSeen(uint32_t max_variable_id, uint32_t max_clause_id) {
    uint32_t seen_size = std::max(max_variable_id,max_clause_id)  + 1;
    debug_print("Creating new seen[] of size: " << seen_size << " and zeroing it.");
    seen_ = new CA_SearchState[seen_size];
    seen_byte_size_ = sizeof(CA_SearchState) * (seen_size);
    clearArrays();
  }

  void clearArrays() {
    memset(seen_, 0, seen_byte_size_);
  }

  // At this point exploreRemainingCompOf has been called already which
  // set up search_stack_, seen[] etc. so this is now quite easy.
  Component *makeComponentFromState(const uint32_t search_stack_size) {
    debug_print(COLREDBG << __PRETTY_FUNCTION__ << " start.");
    Component *p_new_comp = new Component();
    p_new_comp->reserveSpace(search_stack_size, super_comp().numLongClauses());
    current_comp_for_caching_.clear();

    // Fill variables in new comp
    for (auto v_it = super_comp().varsBegin(); *v_it != varsSENTINEL;  v_it++)
      if (var_seen(*v_it)) { //we have to put a var into our comp
        p_new_comp->addVar(*v_it);
        current_comp_for_caching_.addVar(*v_it);
        setVar_in_other_comp(*v_it);
      }
    p_new_comp->closeVariableData();
    current_comp_for_caching_.closeVariableData();

    // Fill clauses in new comp
    for (auto it_cl = super_comp().clsBegin(); *it_cl != clsSENTINEL; it_cl++)
      if (clause_seen(*it_cl)) {
        p_new_comp->addCl(*it_cl);
        if (!clause_all_lits_set(*it_cl)) current_comp_for_caching_.addCl(*it_cl);
        setClause_in_other_comp(*it_cl);
      }
    p_new_comp->closeClauseData();
    current_comp_for_caching_.closeClauseData();
    p_new_comp->setNumBinCls(num_bin_cls);

    debug_print(COLREDBG << __PRETTY_FUNCTION__ << " finish." <<
        " New comp vars: " << p_new_comp->nVars() <<
        " long cls:" << p_new_comp->numLongClauses());
    return p_new_comp;
  }

  Component current_comp_for_caching_;
  uint32_t num_bin_cls = 0;

private:
  Component const* p_super_comp_;
  StackLevel *p_stack_level_;
  CA_SearchState* seen_ = nullptr; // all variables and all clause IDXs can be indexed here
  uint32_t seen_byte_size_ = 0;
};

#endif /* COMPONENT_ARCHETYPE_H_ */
