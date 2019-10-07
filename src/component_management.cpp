/*
 * component_management.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#include "component_management.h"

void ComponentManager::initialize(LiteralIndexedVector<Literal> & literals,
    vector<LiteralID> &lit_pool, unsigned int num_variables) {
      cout << "Doing Initialization "<< endl;

  ana_.initialize(literals, lit_pool);
  // BEGIN CACHE INIT
  CacheableComponent::adjustPackSize(ana_.max_variable_id(), ana_.max_clause_id());

  component_stack_.clear();
  component_stack_.reserve(ana_.max_variable_id() + 2);
  component_stack_.push_back(new Component());
  component_stack_.push_back(new Component());
  assert(component_stack_.size() == 2);
  component_stack_.back()->createAsDummyComponent(ana_.max_variable_id(),
      ana_.max_clause_id());
  if (config_.usecachetencoding || config_.useIsomorphicComponentCaching){
    component_stack_.back()->setisomorphism(false);
    component_stack_.back()->createAsDummyComponentEncoding(ana_.max_variable_id(),
    ana_.max_clause_id(),ana_.clsidtoofs(),lit_pool);
  }
  cache_.init(keyforHash,*component_stack_.back());
  cachescore_.reserve(num_variables + 5);
  for (unsigned i = 0 ; i < (num_variables + 5); i++){
    cachescore_.push_back(0);
  }
  // cout << "cachescore size "<< cachescore_.size()<<endl;
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
