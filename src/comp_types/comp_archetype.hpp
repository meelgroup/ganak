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

#include "../common.hpp"
#include "comp.hpp"
#include "cacheable_comp.hpp"

using std::cout;
using std::endl;

namespace GanakInt {

// State values for variables found during comp analysis (CA)
#define   CA_VAR_IN_SUP_COMP_UNVISITED  1
#define   CA_VAR_VISITED 2
#define   CA_VAR_IN_PEER_COMP  4
// 1 + 2 + 4 == 7
#define   CA_VAR_MASK  7

#define   CA_CL_IN_SUP_COMP_UNVISITED  8
#define   CA_CL_VISITED 16
#define   CA_CL_IN_PEER_COMP  32
/* #define   WHATEVER  64 */
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
  }

  const Comp &super_comp() {
    return *super_comp_ptr;
  }

  StackLevel& stack_level() {
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

  void set_var_in_peer_comp(const uint32_t v) {
    data[v] = CA_VAR_IN_PEER_COMP | (data[v] & CA_CL_MASK);
  }

  bool var_in_peer_comp(const uint32_t v) const {
    return CA_VAR_IN_PEER_COMP & (data[v] & CA_CL_MASK);
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
    debug_print("Creating new data[] of size: " << data_sz << " and zeroing it.");
    data = new uint8_t[data_sz];
    clear_data();
  }

  void clear_data() {
    num_long_cls = 0;
    num_bin_cls = 0;
    memset(data, 0, data_sz);
  }

  // At this point explore_comp has been called already which
  // set up search_stack_, data[] etc. so this is now quite easy.
  Comp* make_comp(const uint32_t comp_vars_size);

  uint32_t num_long_cls = 0;
  uint32_t num_bin_cls = 0;
private:
  Comp const* super_comp_ptr;
  StackLevel *stack_lvl_ptr;
  uint8_t* data = nullptr; // all variables and all clause IDXs can be indexed here
  uint32_t data_sz = 0;
};

}
