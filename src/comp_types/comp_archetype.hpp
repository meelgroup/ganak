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

using std::cout;
using std::endl;

// State values for variables found during comp
// analysis (CA)
using CA_SearchState = uint8_t;
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

// There is exactly ONE of this. Inside CompAnalyzer, which is inside CompManager, which is inside Solver
class CompArchetype {
public:
  CompArchetype() = default;
  ~CompArchetype() { delete[] seen_; }
  CompArchetype(StackLevel &stack_level, const Comp &super_comp) :
      p_super_comp_(&super_comp), p_stack_level_(&stack_level) {
  }

  // called every time we want to deal with a new component
  void re_initialize(StackLevel &stack_level, const Comp &super_comp) {
    debug_print("Reinitializing seen[] to all-zero in CompArchetype");
    p_super_comp_ = &super_comp;
    p_stack_level_ = &stack_level;
    clear_arrays();
    current_comp_for_caching_.reserve_space(super_comp.nVars(),super_comp.num_long_cls());
    BUDDY_DO(num_bin_cls = 0);
  }

  const Comp &super_comp() {
    return *p_super_comp_;
  }

  StackLevel & stack_level() {
    return *p_stack_level_;
  }

  void set_var_in_sup_comp_unseen(const uint32_t v) {
    seen_[v] = CA_VAR_IN_SUP_COMP_UNSEEN | (seen_[v] & CA_CL_MASK);
  }

  void setClause_in_sup_comp_unseen(const ClauseIndex cl) {
    seen_[cl] = CA_CL_IN_SUP_COMP_UNSEEN | (seen_[cl] & CA_VAR_MASK);
  }

  void set_var_nil(const uint32_t v) {
    seen_[v] &= CA_CL_MASK;
  }

  void setClause_nil(const ClauseIndex cl) {
    seen_[cl] &= CA_VAR_MASK;
  }

  void set_var_seen(const uint32_t v) {
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

  void set_var_in_other_comp(const uint32_t v) {
    seen_[v] = CA_VAR_IN_OTHER_COMP | (seen_[v] & CA_CL_MASK);
  }

  void setClause_in_other_comp(const ClauseIndex cl) {
    seen_[cl] = CA_CL_IN_OTHER_COMP | (seen_[cl] & CA_VAR_MASK);
  }

  bool var_seen(const uint32_t v) const {
    return seen_[v] & CA_VAR_SEEN;
  }

  bool clause_seen(const ClauseIndex cl) const {
    return seen_[cl] & CA_CL_SEEN;
  }

  bool clause_all_lits_set(const ClauseIndex cl) const {
    return seen_[cl] & CA_CL_ALL_LITS_SET;
  }

  bool var_nil(const uint32_t v) const {
    return (seen_[v] & CA_VAR_MASK) == 0;
  }

  bool clause_nil(const ClauseIndex cl) const {
    return (seen_[cl] & CA_CL_MASK) == 0;
  }

  bool var_unseen_in_sup_comp(uint32_t v) const {
    return seen_[v] & CA_VAR_IN_SUP_COMP_UNSEEN;
  }

  bool clause_unseen_in_sup_comp(ClauseIndex cl) const {
    return seen_[cl] & CA_CL_IN_SUP_COMP_UNSEEN;
  }

  bool var_seen_in_peer_comp(uint32_t v) const {
    return seen_[v] & CA_VAR_IN_OTHER_COMP;
  }

  bool clause_seen_in_peer_comp(ClauseIndex cl) const {
    return seen_[cl] & CA_CL_IN_OTHER_COMP;
  }

  // called exactly once during lifetime of counter
  void init_seen(uint32_t max_var_id, uint32_t max_cl_id) {
    uint32_t seen_size = std::max(max_var_id,max_cl_id)  + 1;
    debug_print("Creating new seen[] of size: " << seen_size << " and zeroing it.");
    seen_ = new CA_SearchState[seen_size];
    seen_byte_size_ = sizeof(CA_SearchState) * (seen_size);
    clear_arrays();
  }

  void clear_arrays() { memset(seen_, 0, seen_byte_size_); }

  // At this point exploreRemainingCompOf has been called already which
  // set up search_stack_, seen[] etc. so this is now quite easy.
  Comp *makeCompFromState(const uint32_t search_stack_size) {
    debug_print(COLREDBG << __PRETTY_FUNCTION__ << " start.");
    Comp *p_new_comp = new Comp();
    p_new_comp->reserve_space(search_stack_size, super_comp().num_long_cls());
    current_comp_for_caching_.clear();

    // Fill variables in new comp
    for (auto v_it = super_comp().vars_begin(); *v_it != sentinel;  v_it++)
      if (var_seen(*v_it)) { //we have to put a var into our comp
        p_new_comp->add_var(*v_it);
        current_comp_for_caching_.add_var(*v_it);
        set_var_in_other_comp(*v_it);
      }
    p_new_comp->close_vars_data();
    current_comp_for_caching_.close_vars_data();

    // Fill clauses in new comp
    for (auto it_cl = super_comp().cls_begin(); *it_cl != sentinel; it_cl++)
      if (clause_seen(*it_cl)) {
        p_new_comp->add_cl(*it_cl);
        if (!clause_all_lits_set(*it_cl)) current_comp_for_caching_.add_cl(*it_cl);
        setClause_in_other_comp(*it_cl);
      }
    p_new_comp->close_cls_data();
    current_comp_for_caching_.close_cls_data();
    BUDDY_DO(p_new_comp->setNumBinCls(num_bin_cls/2));

    debug_print(COLREDBG << __PRETTY_FUNCTION__ << " finish." <<
        " New comp vars: " << p_new_comp->nVars() <<
        " long cls:" << p_new_comp->numLongClauses());
    return p_new_comp;
  }

  Comp current_comp_for_caching_;
#ifdef BUDDY_ENABLED
  uint32_t num_bin_cls = 0;
#endif

private:
  Comp const* p_super_comp_;
  StackLevel *p_stack_level_;
  CA_SearchState* seen_ = nullptr; // all variables and all clause IDXs can be indexed here
  uint32_t seen_byte_size_ = 0;
};

#endif /* COMPONENT_ARCHETYPE_H_ */
