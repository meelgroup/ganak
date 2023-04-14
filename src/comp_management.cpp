/*
 * comp_management.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#include "comp_management.h"
#include "solver.h"
#include <iomanip>

// Initialized exactly once when Solver is created.
//   it also inits the included analyzer called "ana_"
void ComponentManager::initialize(LiteralIndexedVector<LitWatchList> & literals,
    vector<Lit> &lit_pool, uint32_t nVars){
  assert(comp_stack_.empty());

  ana_.initialize(literals, lit_pool);
  sz = BasePackedComponent::calcPackSize(ana_.max_variable_id(), ana_.max_clause_id());

  //Add dummy comp
  comp_stack_.push_back(new Component());

  //Add full comp
  comp_stack_.push_back(new Component());
  assert(comp_stack_.size() == 2);
  comp_stack_.back()->createStartingComponent(
      ana_.max_variable_id()
      , ana_.max_clause_id());

  cache_.init(*comp_stack_.back(), randomseedforCLHASH);
  for (uint32_t i = 0 ; i < nVars + 1; i++) cachescore_.push_back(0);
}

void ComponentManager::removeAllCachePollutionsOf(const StackLevel &top) {
  // all processed comps are found in
  // [top.currentRemainingComponent(), comp_stack_.size())
  // first, remove the list of descendants from the father
  assert(top.remaining_comps_ofs() <= comp_stack_.size());
  assert(top.super_comp() != 0);
  assert(cache_.hasEntry(getSuperComponentOf(top).id()));

  if (top.remaining_comps_ofs() == comp_stack_.size()) return;

  for (uint32_t u = top.remaining_comps_ofs(); u < comp_stack_.size(); u++) {
    assert(cache_.hasEntry(comp_stack_[u]->id()));
    stats.cache_pollutions_removed += cache_.cleanPollutionsInvolving(comp_stack_[u]->id());
  }
  stats.cache_pollutions_called++;

  SLOW_DEBUG_DO(cache_.test_descendantstree_consistency());
}

double ComponentManager::get_comp_score(StackLevel &top)
{
  const Component& super_comp = getSuperComponentOf(top);

  // This reinitializes archetype, sets up seen[] or all cls&vars unseen (if unset), etc.
  ana_.setupAnalysisContext(top, super_comp);

  double score = 1;
  uint32_t num_comps = 1;
  double sizes_mult = 1;
  uint32_t max_comp_size = 1;
  for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++) {
    if (ana_.isUnseenAndActive(*vt)) {
      bool multi_var_comp = ana_.exploreRemainingCompOf(*vt);
      if (multi_var_comp) {
        Component *p_new_comp = ana_.makeComponentFromArcheType();
        max_comp_size = std::max(p_new_comp->nVars(), max_comp_size);
        delete p_new_comp;
        sizes_mult *= p_new_comp->nVars();
      }
      num_comps++;
    }
  }
  score=num_comps*sizes_mult/(double)max_comp_size;
  /* cout << "v: " << std::setw(4) << top.getbranchvar() */
  /*   << " dec_lev: " << solver_->dec_level() << " Num comps: " */
  /*   << num_comps << " max_comp_size: " << max_comp_size */
  /*   << " score: " << score << endl; */
  return score;
}

// This creates comps
void ComponentManager::recordRemainingCompsFor(StackLevel &top)
{
  const Component& super_comp = getSuperComponentOf(top);
  const uint32_t new_comps_start_ofs = comp_stack_.size();

  // This reinitializes archetype, sets up seen[] or all cls&vars unseen (if unset), etc.
  // Also zeroes out frequency_scores(!)
  ana_.setupAnalysisContext(top, super_comp);

  for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++) {
    print_debug("Going to NEXT var that's unseen & active in this component... if it exists. Var: " << *vt);
    if (ana_.isUnseenAndActive(*vt) && ana_.exploreRemainingCompOf(*vt)) {
      Component *p_new_comp = ana_.makeComponentFromArcheType();
      CacheableComponent *packed_comp = NULL;
      if (config_.do_pcc) {
#ifdef DOPCC
        packed_comp = new CacheableComponent(randomseedforCLHASH, ana_.getArchetype().current_comp_for_caching_, sz);
        packed_comp->finish_hashing(packed_comp->SizeInBytes(sz), packed_comp->nVars(sz));
#else
        exit(-1);
#endif
      } else {
        packed_comp = new CacheableComponent(ana_.getArchetype().current_comp_for_caching_, sz);
      }
      solver_->comp_size_queue.push(packed_comp->nVars(sz));

      // Check if new comp is already in cache
      if (solver_->dec_level() > 0 && p_new_comp->nVars() > 0)
        stats.comp_size_per_depth.push((double)p_new_comp->nVars()/(double)solver_->dec_level());
      if (!cache_.manageNewComponent(top, *packed_comp)) {
        stats.cache_hits_misses.push(p_new_comp->nVars());
        comp_stack_.push_back(p_new_comp);
        p_new_comp->set_id(cache_.storeAsEntry(*packed_comp, super_comp.id()));
#ifdef VERBOSE_DEBUG
        cout << COLYEL2 "New comp. ID: " << p_new_comp->id()
            << " num vars: " << p_new_comp->nVars() << " vars: ";
        for(auto v = p_new_comp->varsBegin(); *v != varsSENTINEL; v++) cout << *v << " ";
        cout << endl;
#endif
      } else {
        stats.cache_hits_misses.push(0);
        if (config_.do_cache_score) {
          stats.numcachedec_++;
          if (stats.numcachedec_ % 128 == 0) rescale_cache_scores();
          for (vector<VariableIndex>::const_iterator it = p_new_comp->varsBegin();
              *it != varsSENTINEL; it++) {
            cachescore_[*it] -= 1;
          }
        }

#ifdef VERBOSE_DEBUG
        cout << COLYEL2 "Component already in cache."
            << " num vars: " << p_new_comp->nVars() << " vars: ";
        for(auto v = p_new_comp->varsBegin(); *v != varsSENTINEL; v++) cout << *v << " ";
        cout << endl;
#endif
        delete packed_comp;
        delete p_new_comp;
      }
    }
  }

  print_debug("We now set the unprocessed_comps_end_ in 'top' to comp_stack_.size(): " << comp_stack_.size() << ", while top.remaining_comps_ofs(): " << top.remaining_comps_ofs());
  top.set_unprocessed_comps_end(comp_stack_.size());
  sortComponentStackRange(new_comps_start_ofs, comp_stack_.size());
}
