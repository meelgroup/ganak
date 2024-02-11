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

#include "comp_management.hpp"
#include "counter.hpp"
#include <iomanip>

// Initialized exactly once when Counter is created.
//   it also inits the included analyzer called "ana_"
void CompManager::initialize(LiteralIndexedVector<LitWatchList> & watches,
    const ClauseAllocator* _alloc, const vector<ClauseOfs>& long_irred_cls, uint32_t nVars){
  assert(comp_stack_.empty());

  ana_.initialize(watches, _alloc, long_irred_cls);

  //Add dummy comp
  comp_stack_.push_back(new Comp());

  //Add full comp
  comp_stack_.push_back(new Comp());
  assert(comp_stack_.size() == 2);
  comp_stack_.back()->create_init_comp(ana_.get_max_var() , ana_.get_max_clid());
  cache_.init(*comp_stack_.back(), hash_seed);
  for (uint32_t i = 0 ; i < nVars + 1; i++) cache_hit_score.push_back(0);
}

void CompManager::removeAllCachePollutionsOfIfExists(const StackLevel &top) {
  assert(top.remaining_comps_ofs() <= comp_stack_.size());
  assert(top.super_comp() != 0);
  if (cache_.hasEntry(getSuperCompOf(top).id())) removeAllCachePollutionsOf(top);
}

void CompManager::removeAllCachePollutionsOf(const StackLevel &top) {
  // all processed comps are found in
  // [top.currentRemainingComp(), comp_stack_.size())
  // first, remove the list of descendants from the father
  assert(top.remaining_comps_ofs() <= comp_stack_.size());
  assert(top.super_comp() != 0);
  assert(cache_.hasEntry(getSuperCompOf(top).id()));

  if (top.remaining_comps_ofs() == comp_stack_.size()) return;

  for (uint32_t u = top.remaining_comps_ofs(); u < comp_stack_.size(); u++) {
    assert(cache_.hasEntry(comp_stack_[u]->id()));
    stats.cache_pollutions_removed += cache_.cleanPollutionsInvolving(comp_stack_[u]->id());
  }
  stats.cache_pollutions_called++;

  /* SLOW_DEBUG_DO(cache_.test_descendantstree_consistency()); */
}

// This creates potential component, checks if it's already in the
// cache, and if so, uses that, otherwise, it creates it
// and adds it to the component stack
void CompManager::recordRemainingCompsFor(StackLevel &top)
{
  const Comp& super_comp = getSuperCompOf(top);
  const uint32_t new_comps_start_ofs = comp_stack_.size();

  // This reinitializes archetype, sets up seen[] or all cls&vars unseen (if unset), etc.
  // Also zeroes out frequency_scores(!)
  ana_.setupAnalysisContext(top, super_comp);

  for (auto vt = super_comp.vars_begin(); *vt != sentinel; vt++) {
    debug_print("Going to NEXT var that's unseen & set in this component... if it exists. Var: " << *vt);
    if (ana_.isUnseenAndSet(*vt) &&
        ana_.exploreRemainingCompOf(*vt)) {

      // Actually makes both a component returned, AND an current_comp_for_caching_ in
      //        Archetype -- BUT, this current_comp_for_caching_ only contains a clause
      //        in case  at least one lit in it is unknown
      Comp *p_new_comp = ana_.makeCompFromArcheType();
      CacheableComp packed_comp(hash_seed, ana_.get_archetype().current_comp_for_caching_);

      // Update stats
      solver_->comp_size_q.push(p_new_comp->nVars());
      stats.comp_size_times_depth_q.push(p_new_comp->nVars()*(solver_->dec_level()/20U+1));

      // Check if new comp is already in cache
      if (!cache_.manageNewComp(top, p_new_comp->nVars(), packed_comp)) {
        stats.cache_hits_misses_q.push(0);
        comp_stack_.push_back(p_new_comp);
        p_new_comp->set_id(cache_.storeAsEntry(packed_comp, super_comp.id()));
#ifdef VERBOSE_DEBUG
        cout << COLYEL2 "New comp. ID: " << p_new_comp->id()
            << " num vars: " << p_new_comp->nVars() << " vars: ";
        all_vars_in_comp(p_new_comp, v) cout << *v << " ";
        cout << endl;
#endif
      } else {
        stats.cache_hits_misses_q.push(p_new_comp->nVars());
        if (conf.do_cache_hit_scores) {
          stats.numcachedec_++;
          bump_cache_hit_score(*p_new_comp);
        }

#ifdef VERBOSE_DEBUG
        cout << COLYEL2 "Comp already in cache."
            << " num vars: " << p_new_comp->nVars() << " vars: ";
        all_vars_in_comp(p_new_comp, v) cout << *v << " ";
        cout << endl;
#endif
        delete p_new_comp;
      }
    }
  }

  debug_print("We now set the unprocessed_comps_end_ in 'top' to comp_stack_.size(): " << comp_stack_.size() << ", while top.remaining_comps_ofs(): " << top.remaining_comps_ofs());
  top.set_unprocessed_comps_end(comp_stack_.size());
  sortCompStackRange(new_comps_start_ofs, comp_stack_.size());
}
