/*
 * stack.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#ifndef STACK_H_
#define STACK_H_

#include <gmpxx.h>

class StackLevel {
  /// active Component, once initialized, it should not change
  const unsigned super_component_ = 0;
  // branch
  bool active_branch_ = false;

  // offset in the literal stack where to store set lits
  const unsigned literal_stack_ofs_ = 0;

  //  Solutioncount
  mpz_class branch_model_count_[2] = {0,0};
  bool branch_found_unsat_[2] = {false,false};

  /// remaining Components

  // the start offset in the component stack for
  // the remaining components in this decision level
  // all remaining components can hence be found in
  // [remaining_components_ofs_, "nextLevel".remaining_components_begin_)
  const unsigned remaining_components_ofs_ = 0;

  // boundary of the stack marking which components still need to be processed
  // all components to be processed can be found in
  // [remaining_components_ofs_, unprocessed_components_end_)
  // also, all processed, can be found
  // in [unprocessed_components_end_, component_stack.size())
  unsigned unprocessed_components_end_ = 0;

  unsigned branch_variable_ = 0;
public:

  bool hasUnprocessedComponents() {
    assert(unprocessed_components_end_ >= remaining_components_ofs_);
    return unprocessed_components_end_ > remaining_components_ofs_;
  }
  void nextUnprocessedComponent() {
    assert(unprocessed_components_end_ > remaining_components_ofs_);
    unprocessed_components_end_--;
  }

  void resetRemainingComps() {
    unprocessed_components_end_ = remaining_components_ofs_;
  }
  unsigned super_component() {
    return super_component_;
  }
  unsigned remaining_components_ofs() {
    return remaining_components_ofs_;
  }
  void set_unprocessed_components_end(unsigned end) {
    unprocessed_components_end_ = end;
    assert(remaining_components_ofs_ <= unprocessed_components_end_);
  }

  StackLevel(unsigned super_comp, unsigned lit_stack_ofs,
      unsigned comp_stack_ofs) :
      super_component_(super_comp),
      literal_stack_ofs_(lit_stack_ofs),
      remaining_components_ofs_(comp_stack_ofs),
      unprocessed_components_end_(comp_stack_ofs) {
    assert(super_comp < comp_stack_ofs);
  }

  unsigned currentRemainingComponent() {
    assert(remaining_components_ofs_ <= unprocessed_components_end_ - 1);
    return unprocessed_components_end_ - 1;
  }
  bool isSecondBranch() {
    return active_branch_;
  }

  void changeBranch() {
    active_branch_ = true;
  }

  bool anotherCompProcessible() {
    return (!branch_found_unsat()) && hasUnprocessedComponents();
  }

  unsigned literal_stack_ofs() {
    return literal_stack_ofs_;
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
  void includeSolution(unsigned solutions) {
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

  bool branch_found_unsat() {
    return branch_found_unsat_[active_branch_];
  }
  void mark_branch_unsat() {
    branch_found_unsat_[active_branch_] = true;
  }

  const mpz_class getBranchSols() const {
    return branch_model_count_[active_branch_];
  }

  unsigned getbranchvar(){
    return branch_variable_;
  }

  void setbranchvariable(unsigned max_score_var){
    branch_variable_ = max_score_var;
  }
//  void set_both_branches_unsat(){
//	  branch_found_unsat_[0] =
//			  branch_found_unsat_[1] = true;
//	  branch_model_count_[0] = branch_model_count_[1] = 0;
//	  active_branch_ = 1;
//  }
  const mpz_class getTotalModelCount() const {
    return branch_model_count_[0] + branch_model_count_[1];
  }
};

class DecisionStack: public vector<StackLevel> {
  unsigned int failed_literal_test_active = 0;
public:

  //begin for implicit BCP
  void startFailedLitTest() {
    failed_literal_test_active = true;
  }
  void stopFailedLitTest() {
    failed_literal_test_active = false;
  }
  //end for implicit BCP

  StackLevel &top() {
    assert(size() > 0);
    return back();
  }
  int get_decision_level() const {
    assert(size() > 0);
    return size() - 1 + failed_literal_test_active;
  } // 0 means pre-1st-decision

};



#endif /* STACK_H_ */
