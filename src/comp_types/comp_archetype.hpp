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

#include <cstring>
#include <iostream>
#include <limits>

#include "../common.hpp"
#include "comp.hpp"
#include "cacheable_comp.hpp"

using std::numeric_limits;
using std::cout;
using std::endl;

namespace GanakInt {

class StackLevel;

// There is exactly ONE of this. Inside CompAnalyzer, which is inside CompManager, which is inside Solver
class CompArchetype {
public:
  CompArchetype() = default;
  ~CompArchetype() { delete[] raw_data; }

  // called every time we want to deal with a new component
  void re_initialize(StackLevel &stack_level, const Comp &super_comp) {
    debug_print("Reinitializing data to all-zero in CompArchetype");
    super_comp_ptr = &super_comp;
    stack_lvl_ptr = &stack_level;
    clear_data();
  }

  const Comp &super_comp() { return *super_comp_ptr; }
  StackLevel& stack_level() { return *stack_lvl_ptr; }

  void set_var_visited(const uint32_t v) { v_data[v] = tstamp; }
  void set_clause_visited(const ClauseIndex cl) { cl_data[cl] = tstamp; }

  void set_var_clear(const uint32_t v) { v_data[v] = 0; }
  void set_cl_clear(const uint32_t cl) { cl_data[cl] = 0; }

  void set_var_in_sup_comp_unvisited(const uint32_t v) { v_data[v] = tstamp-1; }
  void set_clause_in_sup_comp_unvisited(const ClauseIndex cl) { cl_data[cl] = tstamp-1; }
  bool var_unvisited_in_sup_comp(uint32_t v) const { return v_data[v] == tstamp-1; }
  bool clause_unvisited_in_sup_comp(ClauseIndex cl) const { return cl_data[cl] == tstamp-1; }

  bool var_visited(const uint32_t v) const { return v_data[v] == tstamp; }
  bool clause_visited(const ClauseIndex cl) const { return cl_data[cl] == tstamp; }

  // called exactly once during lifetime of counter
  void init_data(uint32_t max_var_id, uint32_t max_cl_id) {
    assert(tstamp == 0);
    data_sz = max_var_id +1 + max_cl_id + 1;
    debug_print("Creating new data[] of size: " << data_sz << " and zeroing it.");
    raw_data = new uint64_t[data_sz];
    memset(raw_data, 0, data_sz * sizeof(uint64_t));
    v_data = raw_data;
    cl_data = raw_data + max_var_id + 1;
    clear_data();
  }

  void clear_data() {
    num_long_cls = 0;
    num_bin_cls = 0;
    tstamp+=2;
  }

  // At this point explore_comp has been called already which
  // set up search_stack_, data[] etc. so this is now quite easy.
  Comp* make_comp(const uint32_t comp_vars_size);

  uint32_t num_long_cls = 0;
  uint32_t num_bin_cls = 0;
private:
  uint64_t tstamp = 0;
  Comp const* super_comp_ptr;
  StackLevel *stack_lvl_ptr;
  uint64_t* raw_data = nullptr;
  uint64_t* cl_data = nullptr;
  uint64_t* v_data = nullptr;
  uint32_t data_sz = 0;
};

}
