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
#include "structures.hpp"
using std::vector;
using std::cout;
using std::endl;

template<typename T>
class StackLevel {
public:
  StackLevel(uint32_t super_comp, uint32_t comp_stack_ofs) :
      super_comp_(super_comp),
      remaining_comps_ofs_(comp_stack_ofs),
      unprocessed_comps_end_(comp_stack_ofs) {
    assert(super_comp < comp_stack_ofs);
  }
  uint32_t var = 0;
private:

  /// active Comp, once initialized, it should not change
  const uint32_t super_comp_ = 0;
  // branch (i.e. left = false/right = true)
  bool active_branch_ = false;

  //  Solution count
  T branch_model_count_[2] = {0,0};
  bool branch_found_unsat_[2] = {false,false};

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
  bool hasUnprocessedComps() const {
    assert(unprocessed_comps_end_ >= remaining_comps_ofs_);
    return unprocessed_comps_end_ > remaining_comps_ofs_;
  }
  uint32_t numUnprocessedComps() const {
    assert(unprocessed_comps_end_ >= remaining_comps_ofs_);
    return unprocessed_comps_end_ - remaining_comps_ofs_;
  }
  void nextUnprocessedComp() {
    assert(unprocessed_comps_end_ > remaining_comps_ofs_);
    unprocessed_comps_end_--;
  }

  void resetRemainingComps() {
    unprocessed_comps_end_ = remaining_comps_ofs_;
  }

  uint32_t super_comp() const {
    return super_comp_;
  }
  uint32_t getUnprocessedCompsEnd() const {
    return unprocessed_comps_end_;
  }
  uint32_t remaining_comps_ofs() const {
    return remaining_comps_ofs_;
  }
  void set_unprocessed_comps_end(uint32_t end) {
    unprocessed_comps_end_ = end;
    assert(remaining_comps_ofs_ <= unprocessed_comps_end_);
  }

  uint32_t currentRemainingComp() const {
    assert(remaining_comps_ofs_ <= unprocessed_comps_end_ - 1);
    return unprocessed_comps_end_ - 1;
  }
  bool is_right_branch() const {
    return active_branch_;
  }

  void change_to_right_branch() {
    assert(active_branch_ == false);
    active_branch_ = true;
    SLOW_DEBUG_DO(assert(branch_model_count_[active_branch_] == 0));
  }

  bool another_comp_possible() const {
    return (!branch_found_unsat()) && hasUnprocessedComps();
  }

  template<class T2>
  void include_solution(const T2& solutions) {
    VERBOSE_DEBUG_DO(cout << "incl sol: " << solutions << endl);
#ifdef VERBOSE_DEBUG
    auto before = branch_model_count_[active_branch_];
#endif
    if (branch_found_unsat_[active_branch_]) {
      assert(branch_model_count_[active_branch_] == 0);
      return;
    }
    if (solutions == 0) branch_found_unsat_[active_branch_] = true;
    if (branch_model_count_[active_branch_] == 0) {
      branch_model_count_[active_branch_] = solutions;
    } else {
      branch_model_count_[active_branch_] *= solutions;
    }
    VERBOSE_DEBUG_DO(cout << "now "
        << ((active_branch_) ? "right" : "left")
        << " count is: " << branch_model_count_[active_branch_]
        << " before it was: " << before
        << " var: " << var
        << endl);
  }

  bool branch_found_unsat() const {
    return branch_found_unsat_[active_branch_];
  }
  void mark_branch_unsat() {
    branch_found_unsat_[active_branch_] = true;
  }

  const T& getBranchSols() const {
    return branch_model_count_[active_branch_];
  }

  const T& get_model_side(int side) const {
    return branch_model_count_[side];
  }

  void zero_out_branch_sol() {
    branch_model_count_[active_branch_] = 0;
  }

  void zero_out_all_sol() {
    branch_model_count_[0] = 0;
    branch_model_count_[1] = 0;
  }

  const T getTotalModelCount() const {
    return branch_model_count_[0] + branch_model_count_[1];
  }

  const T& get_left_model_count() const {
    return branch_model_count_[0];
  }
  const T& get_right_model_count() const {
    return branch_model_count_[1];
  }
};

template<typename T>
class DecisionStack: public vector<StackLevel<T>> {
public:

  const StackLevel<T>& top() const{
    assert(vector<StackLevel<T>>::size() > 0);
    return vector<StackLevel<T>>::back();
  }

  StackLevel<T>& top(){
    assert(vector<StackLevel<T>>::size() > 0);
    return vector<StackLevel<T>>::back();
  }

  /// 0 means pre-1st-decision
  int32_t get_decision_level() const {
    assert(vector<StackLevel<T>>::size() > 0);
    return (int)vector<StackLevel<T>>::size() - 1;
  }
};
