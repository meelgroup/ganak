/*
 * comp_management.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#include "comp_management.h"

void ComponentManager::initialize(LiteralIndexedVector<LitWatchList> & literals,
    vector<Lit> &lit_pool, uint32_t nVars){
  assert(comp_stack_.empty());

  ana_.initialize(literals, lit_pool);
  CacheableComponent::adjustPackSize(ana_.max_variable_id(), ana_.max_clause_id());

  //Add dummy comp
  comp_stack_.push_back(new Component());

  //Add full comp
  comp_stack_.push_back(new Component());
  assert(comp_stack_.size() == 2);
  comp_stack_.back()->createStartingComponent(
      ana_.max_variable_id()
      , ana_.max_clause_id());


  cache_.init(*comp_stack_.back(), seedforCLHASH);
  for (uint32_t i = 0 ; i < (nVars + 5); i++){
    cachescore_.push_back(0);
  }
}

void ComponentManager::removeAllCachePollutionsOf(StackLevel &top) {
  // all processed comps are found in
  // [top.currentRemainingComponent(), comp_stack_.size())
  // first, remove the list of descendants from the father
  assert(top.remaining_comps_ofs() <= comp_stack_.size());
  assert(top.super_comp() != 0);
  assert(cache_.hasEntry(getSuperComponentOf(top).id()));

  if (top.remaining_comps_ofs() == comp_stack_.size()) return;

  for (uint32_t u = top.remaining_comps_ofs(); u < comp_stack_.size(); u++) {
    assert(cache_.hasEntry(comp_stack_[u]->id()));
    cache_.cleanPollutionsInvolving(comp_stack_[u]->id());
  }

  SLOW_DEBUG_DO(cache_.test_descendantstree_consistency());
}
