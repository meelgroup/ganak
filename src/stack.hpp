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

#include <gmpxx.h>
#include <cassert>
#include <vector>
#include <iostream>
#include "common.hpp"
using std::vector;
using std::cout;
using std::endl;

namespace GanakInt {

class StackLevel {
public:
  StackLevel(uint32_t super_comp, uint32_t comp_stack_ofs, bool _is_indep, uint64_t _tstamp,
      const FG& _fg) :
      fg(_fg),
      tstamp(_tstamp),
      is_indep(_is_indep),
      super_comp_(super_comp),
      remaining_comps_ofs_(comp_stack_ofs),
      unprocessed_comps_end_(comp_stack_ofs) {
    branch_mc[0] = nullptr;
    branch_mc[1] = nullptr;
    assert(super_comp < comp_stack_ofs);
  }
  const FG& fg;
  uint64_t tstamp;
  bool is_indep;

  inline const FF val_or_zero(const bool b) const {
    if (branch_mc[b] == nullptr) return fg->zero();
    return branch_mc[b]->dup();

  }
  inline bool is_zero(const bool b) const {
    return (branch_mc[b] == nullptr || branch_mc[b]->is_zero());
  }
  inline bool is_one(const bool b) const {
    return (branch_mc[b] && branch_mc[b]->is_zero());
  }

  uint32_t var = 0;
  void reset() {
    act_branch = 0;
    branch_unsat[0] = false;
    branch_unsat[1] = false;
    branch_mc[0] = nullptr;
    branch_mc[1] = nullptr;
  }
private:

  /// active Comp, once initialized, it does not change
  const uint32_t super_comp_ = 0;

  // branch (i.e. left = false/right = true)
  bool act_branch = false;

  //  Solution count
  FF branch_mc[2];
  bool branch_unsat[2] = {false,false};

  /// remaining Comps

  // the start offset in the comp stack for
  // the remaining comps in this decision level
  // all remaining comps can hence be found in
  // [remaining_comps_ofs_, "nextLevel".remaining_comps_begin_)
  // SET ONCE, NEVER TOUCHED
  const uint32_t remaining_comps_ofs_ = 0;

  // boundary of the stack marking which comps still need to be processed
  // all comps to be processed can be found in
  // [remaining_comps_ofs_, unprocessed_comps_end_)
  // also, all processed, can be found
  // in [unprocessed_comps_end_, comp_stack.size())
  // KEEPS BEING DECREMENTED, until it reaches remaining_comps_ofs_
  uint32_t unprocessed_comps_end_ = 0;

public:
  bool has_unproc_comps() const {
    assert(unprocessed_comps_end_ >= remaining_comps_ofs_);
    return unprocessed_comps_end_ > remaining_comps_ofs_;
  }
  uint32_t num_unproc_comps() const {
    assert(unprocessed_comps_end_ >= remaining_comps_ofs_);
    return unprocessed_comps_end_ - remaining_comps_ofs_;
  }
  void next_unproc_comp() {
    assert(unprocessed_comps_end_ > remaining_comps_ofs_);
    unprocessed_comps_end_--;
  }
  void reset_remain_comps() { unprocessed_comps_end_ = remaining_comps_ofs_; }
  auto get_unprocessed_comps_end() const { return unprocessed_comps_end_; }
  uint32_t super_comp() const { return super_comp_; }
  bool is_right_branch() const { return act_branch; }
  uint32_t get_unproc_comps_end() const { return unprocessed_comps_end_; }
  uint32_t remaining_comps_ofs() const { return remaining_comps_ofs_; }
  void set_unprocessed_comps_end(uint32_t end) {
    unprocessed_comps_end_ = end;
    assert(remaining_comps_ofs_ <= unprocessed_comps_end_);
  }

  uint32_t curr_remain_comp() const {
    assert(remaining_comps_ofs_ <= unprocessed_comps_end_ - 1);
    return unprocessed_comps_end_ - 1;
  }

  void change_to_right_branch() {
    assert(act_branch == false);
    act_branch = true;
    SLOW_DEBUG_DO(assert(is_zero(act_branch)));
  }

  bool another_comp_possible() const {
    return (!branch_found_unsat()) && has_unproc_comps();
  }

  inline void common_print(const FF& before) {
    cout << "now "
        << ((act_branch) ? "right" : "left")
        << " count is: " << val_or_zero(act_branch)
        << " before it was: " << before
        << " var: " << var
        << " while " << ((!act_branch) ? "right" : "left")
        << " count is: " << branch_mc[!act_branch]
        << endl;
  }

  void include_one_sol() {
    VERBOSE_DEBUG_DO(cout << COLRED << "incl sol: ONE" << COLDEF << " ");
#ifdef VERBOSE_DEBUG
    auto before = val_or_zero(act_branch);
#endif
    if (branch_unsat[act_branch]) {
      VERBOSE_DEBUG_DO(cout << "-> incl sol unsat branch, doing  nothing." << endl);
      assert(is_zero(act_branch));
      return;
    }
    if (!is_indep || is_zero(act_branch))
      branch_mc[act_branch] = fg->one();
    VERBOSE_DEBUG_DO(common_print(before));
  }

  void include_solution(const FF& solutions) {
    VERBOSE_DEBUG_DO(cout << COLRED << "incl sol: " << *solutions << COLDEF << " ");
#ifdef VERBOSE_DEBUG
    auto before = val_or_zero(act_branch);
#endif
    if (branch_unsat[act_branch]) {
      VERBOSE_DEBUG_DO(cout << "-> incl sol unsat branch, doing  nothing." << endl);
      assert(is_zero(act_branch));
      return;
    }

    if (solutions->is_zero()) mark_branch_unsat();
    else {
      if (!is_indep) branch_mc[act_branch] = fg->one();
      else {
        if (!is_zero(act_branch)) {
          if (is_one(act_branch)) branch_mc[act_branch] = solutions->dup();
          else *branch_mc[act_branch] *= *solutions;
        } else branch_mc[act_branch] = solutions->dup();
      }
    }
    VERBOSE_DEBUG_DO(common_print(before));
  }

  void div_solution_left_side(const FF& div_by) {
    VERBOSE_DEBUG_DO(cout << COLRED << "left side div sol: " << *div_by << COLDEF << " " << endl;);
    if (act_branch == 0) return;
#ifdef VERBOSE_DEBUG
    auto before = val_or_zero(0);
#endif
    if (branch_unsat[0]) {
      VERBOSE_DEBUG_DO(cout << "-> left side incl sol unsat branch, doing  nothing." << endl);
      assert(is_zero(0));
      return;
    }

    if (is_indep) {
      assert(!is_zero(0));
      *branch_mc[0] /= *div_by;
    }
    VERBOSE_DEBUG_DO(cout << "now "
        << ((0) ? "right" : "left")
        << " count is: " << *val_or_zero(0)
        << " before it was: " << *before
        << " var: " << var
        << endl);
  }

  bool branch_found_unsat() const { return branch_unsat[act_branch]; }
  void mark_branch_unsat() {
    branch_unsat[act_branch] = true;
    branch_mc[act_branch] = nullptr;
  }
  const FF& get_branch_sols() const { return branch_mc[act_branch]; }
  const FF& get_model_side(int side) const { return branch_mc[side]; }
  void zero_out_branch_sol() { branch_mc[act_branch] = nullptr; }
  FF total_model_count() const {
    if (is_indep) {
      if (is_zero(0)) return val_or_zero(1);
      else if (is_zero(1)) return val_or_zero(0);
      return branch_mc[0]->add(*branch_mc[1]);
    }
    else {
      if (is_zero(0)) return val_or_zero(1);
      else return val_or_zero(0);
    }
  }

  // for cube creation
  bool branch_found_unsat(int side) const { return branch_unsat[side]; }
  FF left_model_count() const { return val_or_zero(0); }
  FF right_model_count() const { return val_or_zero(1); }

  void zero_out_all_sol() {
    branch_mc[0] = nullptr;
    branch_mc[1] = nullptr;
  }
};

class DecisionStack: public vector<StackLevel> {
public:

  const StackLevel& top() const{
    assert(vector<StackLevel>::size() > 0);
    return vector<StackLevel>::back();
  }

  StackLevel& top(){
    assert(vector<StackLevel>::size() > 0);
    return vector<StackLevel>::back();
  }

  /// 0 means pre-1st-decision
  int32_t get_decision_level() const {
    assert(vector<StackLevel>::size() > 0);
    return (int)vector<StackLevel>::size() - 1;
  }
};

}
