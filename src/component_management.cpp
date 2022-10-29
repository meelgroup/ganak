/*
 * component_management.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#include "component_management.h"

void ComponentManager::initialize(LiteralIndexedVector<LitWatchList> & literals,
    vector<Lit> &lit_pool, unsigned num_variables){
  assert(component_stack_.empty());

  ana_.initialize(literals, lit_pool);
  CacheableComponent::adjustPackSize(ana_.max_variable_id(), ana_.max_clause_id());

  //Add dummy component
  component_stack_.push_back(new Component());

  //Add full component
  component_stack_.push_back(new Component());
  assert(component_stack_.size() == 2);
  component_stack_.back()->createAsDummyComponent(
      ana_.max_variable_id()
      , ana_.max_clause_id());


  cache_.init(*component_stack_.back(), seedforCLHASH);
  for (unsigned i = 0 ; i < (num_variables + 5); i++){
    cachescore_.push_back(0);
  }
}

void ComponentManager::removeAllCachePollutionsOf(StackLevel &top) {
  // all processed components are found in
  // [top.currentRemainingComponent(), component_stack_.size())
  // first, remove the list of descendants from the father
  assert(top.remaining_components_ofs() <= component_stack_.size());
  assert(top.super_component() != 0);
  assert(cache_.hasEntry(superComponentOf(top).id()));

  if (top.remaining_components_ofs() == component_stack_.size())
    return;

  for (unsigned u = top.remaining_components_ofs(); u < component_stack_.size();
      u++) {
    assert(cache_.hasEntry(component_stack_[u]->id()));
    cache_.cleanPollutionsInvolving(component_stack_[u]->id());
  }

#ifdef DEBUG
  cache_.test_descendantstree_consistency();
#endif
}
