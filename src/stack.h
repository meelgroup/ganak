/*
 * stack.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */
#pragma once

#include <gmpxx.h>
#include <cassert>
#include <vector>
#include "common.h"
using std::vector;

class StackLevel {
public:
  StackLevel(uint32_t super_comp, uint32_t trail_ofs, uint32_t comp_stack_ofs) :
      super_comp_(super_comp),
      trail_ofs_(trail_ofs),
      remaining_comps_ofs_(comp_stack_ofs),
      unprocessed_comps_end_(comp_stack_ofs) {
    assert(super_comp < comp_stack_ofs);
  }
private:

  /// active Component, once initialized, it should not change
  const uint32_t super_comp_ = 0;
  // branch (i.e. left/right)
  bool active_branch_ = false;

  // offset in the literal stack where to store set lits
  const uint32_t trail_ofs_ = 0;

  //  Solution count
  mpz_class branch_model_count_[2] = {0,0};
  bool branch_found_unsat_[2] = {false,false};

  /// remaining Components

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

  uint32_t branch_variable_ = 0;
public:
  bool on_path_to_target_ = false;

  bool hasUnprocessedComponents() const {
    assert(unprocessed_comps_end_ >= remaining_comps_ofs_);
    return unprocessed_comps_end_ > remaining_comps_ofs_;
  }
  uint32_t numUnprocessedComponents() const {
    assert(unprocessed_comps_end_ >= remaining_comps_ofs_);
    return unprocessed_comps_end_ - remaining_comps_ofs_;
  }
  void nextUnprocessedComponent() {
    assert(unprocessed_comps_end_ > remaining_comps_ofs_);
    unprocessed_comps_end_--;
  }

  void resetRemainingComps() {
    unprocessed_comps_end_ = remaining_comps_ofs_;
  }

  uint32_t super_comp() const {
    return super_comp_;
  }
  uint32_t getUnprocessedComponentsEnd() const {
    return unprocessed_comps_end_;
  }
  uint32_t remaining_comps_ofs() const {
    return remaining_comps_ofs_;
  }
  void set_unprocessed_comps_end(uint32_t end) {
    unprocessed_comps_end_ = end;
    assert(remaining_comps_ofs_ <= unprocessed_comps_end_);
  }

  uint32_t currentRemainingComponent() const {
    assert(remaining_comps_ofs_ <= unprocessed_comps_end_ - 1);
    return unprocessed_comps_end_ - 1;
  }
  bool is_right_branch() const {
    return active_branch_;
  }

  void change_to_right_branch() {
    assert(active_branch_ == false);
    active_branch_ = true;
  }

  bool anotherCompProcessible() const {
    return (!branch_found_unsat()) && hasUnprocessedComponents();
  }

  uint32_t trail_ofs() const {
    return trail_ofs_;
  }
  void includeSolution(const mpz_class &solutions) {
    if (branch_found_unsat_[active_branch_]) {
      assert(branch_model_count_[active_branch_] == 0);
      return;
    }
    if (solutions == 0)
      branch_found_unsat_[active_branch_] = true;
    if (branch_model_count_[active_branch_] == 0)
      branch_model_count_[active_branch_] = solutions;
    else
      branch_model_count_[active_branch_] *= solutions;

  }
  void includeSolution(uint32_t solutions) {
    if (branch_found_unsat_[active_branch_]) {
      assert(branch_model_count_[active_branch_] == 0);
      return;
    }
    if (solutions == 0)
      branch_found_unsat_[active_branch_] = true;
    if (branch_model_count_[active_branch_] == 0)
      branch_model_count_[active_branch_] = solutions;
    else
      branch_model_count_[active_branch_] *= solutions;
  }

  bool branch_found_unsat() const {
    return branch_found_unsat_[active_branch_];
  }
  void mark_branch_unsat() {
    branch_found_unsat_[active_branch_] = true;
  }

  const mpz_class getBranchSols() const {
    return branch_model_count_[active_branch_];
  }

  uint32_t getbranchvar() const {
    return branch_variable_;
  }

  void setbranchvariable(const uint32_t v){
    branch_variable_ = v;
  }

  void setonpath(const bool on_path) {
    on_path_to_target_ = on_path;
  }

  const mpz_class getTotalModelCount() const {
    return branch_model_count_[0] + branch_model_count_[1];
  }
};

class DecisionStack: public vector<StackLevel> {
  uint32_t failed_lit_test_active = 0;
public:

  void startFailedLitTest() {
    failed_lit_test_active = true;
  }
  void stopFailedLitTest() {
    failed_lit_test_active = false;
  }

  const StackLevel &top() const{
    assert(size() > 0);
    return back();
  }

  StackLevel &top(){
    assert(size() > 0);
    return back();
  }

  /// 0 means pre-1st-decision
  uint32_t get_decision_level() const {
    assert(size() > 0);
    return size() - 1 + failed_lit_test_active;
  }
};
