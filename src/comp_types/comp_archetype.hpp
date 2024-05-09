/*
 * comp_archetype.h
 *
 *  Created on: Feb 9, 2013
 *      Author: mthurley
 */

#pragma once

#include <cstring>
#include <algorithm>
#include <iostream>

#include "../primitive_types.hpp"
#include "../common.hpp"
#include "comp.hpp"
#include "cacheable_comp.hpp"

using std::cout;
using std::endl;

// State values for variables found during comp analysis (CA)
#define   CA_VAR_IN_SUP_COMP_UNVISITED  1
#define   CA_VAR_VISITED 2
#define   CA_VAR_IN_PEER_COMP  4
// 1 + 2 + 4 == 7
#define   CA_VAR_MASK  7

#define   CA_CL_IN_SUP_COMP_UNVISITED  8
#define   CA_CL_VISITED 16
#define   CA_CL_IN_PEER_COMP  32
#define   CA_CL_ALL_LITS_UNK  64
// 64+32+16+8 == 120
#define   CA_CL_MASK  120

class StackLevel;

// There is exactly ONE of this. Inside CompAnalyzer, which is inside CompManager, which is inside Solver
class CompArchetype {
public:
  CompArchetype() = default;
  ~CompArchetype() { delete[] data; }
  CompArchetype(StackLevel &stack_level, const Comp &super_comp) :
      super_comp_ptr(&super_comp), stack_lvl_ptr(&stack_level) {
  }

  // called every time we want to deal with a new component
  void re_initialize(StackLevel &stack_level, const Comp &super_comp) {
    debug_print("Reinitializing data to all-zero in CompArchetype");
    super_comp_ptr = &super_comp;
    stack_lvl_ptr = &stack_level;
    clear_data();
    curr_comp.reserve_space(super_comp.nVars(),super_comp.num_long_cls());
    BUDDY_DO(num_bin_cls = 0);
  }

  const Comp &super_comp() {
    return *super_comp_ptr;
  }

  StackLevel & stack_level() {
    return *stack_lvl_ptr;
  }

  void set_var_in_sup_comp_unvisited(const uint32_t v) {
    data[v] = CA_VAR_IN_SUP_COMP_UNVISITED | (data[v] & CA_CL_MASK);
  }

  void set_var_in_sup_comp_unvisited_raw(const uint32_t v) {
    data[v] = CA_VAR_IN_SUP_COMP_UNVISITED;
  }

  void set_clause_in_sup_comp_unvisited(const ClauseIndex cl) {
    data[cl] = CA_CL_IN_SUP_COMP_UNVISITED | (data[cl] & CA_VAR_MASK);
  }

  void clear_var(const uint32_t v) {
    data[v] &= CA_CL_MASK;
  }

  void clear_cl(const ClauseIndex cl) {
    data[cl] &= CA_VAR_MASK;
  }

  // REMOVES from unseen of super, sets visited
  void set_var_visited(const uint32_t v) {
    data[v] = CA_VAR_VISITED | (data[v] & CA_CL_MASK);
  }

  void set_clause_visited(const ClauseIndex cl) {
    clear_cl(cl);
    data[cl] = CA_CL_VISITED | (data[cl] & CA_VAR_MASK);
  }

  void set_clause_visited(const ClauseIndex cl, const bool all_lits_unkn) {
    clear_cl(cl);
    data[cl] = CA_CL_VISITED | (all_lits_unkn?CA_CL_ALL_LITS_UNK:0) | (data[cl] & CA_VAR_MASK);
  }

  void set_var_in_peer_comp(const uint32_t v) {
    data[v] = CA_VAR_IN_PEER_COMP | (data[v] & CA_CL_MASK);
  }

  void set_clause_in_peer_comp(const ClauseIndex cl) {
    data[cl] = CA_CL_IN_PEER_COMP | (data[cl] & CA_VAR_MASK);
  }

  bool var_visited(const uint32_t v) const {
    return data[v] & CA_VAR_VISITED;
  }

  bool clause_visited(const ClauseIndex cl) const {
    return data[cl] & CA_CL_VISITED;
  }

  bool clause_all_lits_unkn(const ClauseIndex cl) const {
    return data[cl] & CA_CL_ALL_LITS_UNK;
  }

  bool var_nil(const uint32_t v) const {
    return (data[v] & CA_VAR_MASK) == 0;
  }

  bool clause_nil(const ClauseIndex cl) const {
    return (data[cl] & CA_CL_MASK) == 0;
  }

  bool var_unvisited_in_sup_comp(uint32_t v) const {
    return data[v] & CA_VAR_IN_SUP_COMP_UNVISITED;
  }

  bool clause_unvisited_in_sup_comp(ClauseIndex cl) const {
    return data[cl] & CA_CL_IN_SUP_COMP_UNVISITED;
  }

  // called exactly once during lifetime of counter
  void init_data(uint32_t max_var_id, uint32_t max_cl_id) {
    data_sz = std::max(max_var_id,max_cl_id)  + 1;
    debug_print("Creating new data[] of size: " << data_size << " and zeroing it.");
    data = new uint8_t[data_sz];
    clear_data();
  }

  void clear_data() { memset(data, 0, data_sz); }

  // At this point exploreRemainingCompOf has been called already which
  // set up search_stack_, data[] etc. so this is now quite easy.
  Comp* make_comp(const uint32_t comp_vars_size) {
    debug_print(COLREDBG << __PRETTY_FUNCTION__ << " start.");
    Comp *p_new_comp = new Comp();
    p_new_comp->reserve_space(comp_vars_size, super_comp().num_long_cls());
    curr_comp.clear();

    // Fill variables in new comp
    for (auto v_it = super_comp().vars_begin(); *v_it != sentinel;  v_it++)
      if (var_visited(*v_it)) { //we have to put a var into our comp
        p_new_comp->add_var(*v_it);
        curr_comp.add_var(*v_it);
        set_var_in_peer_comp(*v_it);
      }
    p_new_comp->close_vars_data();
    curr_comp.close_vars_data();

    // Fill (long) clause IDs in new comp
    for (auto it_cl = super_comp().cls_begin(); *it_cl != sentinel; it_cl++)
      if (clause_visited(*it_cl)) {
        p_new_comp->add_cl(*it_cl);
        if (!clause_all_lits_unkn(*it_cl)) curr_comp.add_cl(*it_cl);
        set_clause_in_peer_comp(*it_cl);
      }
    p_new_comp->close_cls_data();
    curr_comp.close_cls_data();
    BUDDY_DO(p_new_comp->setNumBinCls(num_bin_cls/2));

    debug_print(COLREDBG << __PRETTY_FUNCTION__ << " finish." <<
        " New comp vars: " << p_new_comp->nVars() <<
        " long cls:" << p_new_comp->num_long_cls());
    return p_new_comp;
  }

  Comp curr_comp;
#ifdef BUDDY_ENABLED
  uint32_t num_bin_cls = 0;
#endif

private:
  Comp const* super_comp_ptr;
  StackLevel *stack_lvl_ptr;
  uint8_t* data = nullptr; // all variables and all clause IDXs can be indexed here
  uint32_t data_sz = 0;
};
