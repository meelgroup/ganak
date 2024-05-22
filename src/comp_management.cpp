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
#include "common.hpp"
#include "counter.hpp"
#include "mpreal.h"

template class CompManager<mpz_class>;
template class CompManager<mpfr::mpreal>;

template<typename T>
void CompManager<T>::removeAllCachePollutionsOfIfExists(const StackLevel<T> &top) {
  assert(top.remaining_comps_ofs() <= comp_stack.size());
  assert(top.super_comp() != 0);

  for (uint32_t u = top.remaining_comps_ofs(); u < comp_stack.size(); u++) {
    if (!cache.exists(comp_stack[u]->id())) continue;
    stats.cache_pollutions_removed += cache.clean_pollutions_involving(comp_stack[u]->id());
  }
  stats.cache_pollutions_called++;
}

template<typename T>
void CompManager<T>::removeAllCachePollutionsOf(const StackLevel<T> &top) {
  // all processed comps are found in
  // [top.currentRemainingComp(), comp_stack.size())
  // first, remove the list of descendants from the father
  assert(top.remaining_comps_ofs() <= comp_stack.size());
  assert(top.super_comp() != 0);
  assert(cache.exists(get_super_comp(top).id()));

  for (uint32_t u = top.remaining_comps_ofs(); u < comp_stack.size(); u++) {
    assert(cache.exists(comp_stack[u]->id()));
    stats.cache_pollutions_removed += cache.clean_pollutions_involving(comp_stack[u]->id());
  }
  stats.cache_pollutions_called++;

  /* SLOW_DEBUG_DO(cache.test_descendantstree_consistency()); */
}

// This creates potential component, checks if it's already in the
// cache, and if so, uses that, otherwise, it creates it
// and adds it to the component stack
template<typename T>
void CompManager<T>::recordRemainingCompsFor(StackLevel<T> &top)
{
  const Comp& super_comp = get_super_comp(top);
  const uint32_t new_comps_start_ofs = comp_stack.size();

  // This reinitializes archetype, sets up seen[] or all cls&vars unvisited (if unset), etc.
  // Also zeroes out frequency_scores(!)
  ana.setup_analysis_context(top, super_comp);

  all_vars_in_comp(super_comp, vt) {
    debug_print("Going to NEXT var that's unvisited & set in this component... if it exists. Var: " << *vt);
    if (ana.var_unvisited_sup_comp(*vt) && ana.explore_comp(*vt)) {
      // Actually makes both a component returned, AND an current_comp_for_caching_ in
      //        Archetype -- BUT, this current_comp_for_caching_ only contains a clause
      //        in case  at least one lit in it is unknown
      Comp *p_new_comp = ana.make_comp_from_archetype();
      CacheableComp<T> packed_comp(hash_seed, ana.get_archetype().cache_comp);

      // Update stats
      counter->depth_q.push(counter->decision_level());
      counter->comp_size_q.push(p_new_comp->nVars());

      // TODO Yash: count it 1-by-1 in case the number of variables & clauses is small
      //       essentially, brute-forcing the count
      if (!cache.find_comp_and_incorporate_cnt(top, p_new_comp->nVars(), packed_comp)) {
        // Cache miss
        comp_stack.push_back(p_new_comp);

        p_new_comp->set_id(cache.new_comp(packed_comp, super_comp.id()));
        stats.incorporate_cache_store(packed_comp, p_new_comp->nVars());
#ifdef VERBOSE_DEBUG
        cout << COLYEL2 "New comp. ID: " << p_new_comp->id()
            << " num vars: " << p_new_comp->nVars() << " vars: ";
        all_vars_in_comp(*p_new_comp, v) cout << *v << " ";
        cout << endl;
#endif
      } else {
        // Cache hit
#ifdef VERBOSE_DEBUG
        cout << COLYEL2 "Comp already in cache."
            << " num vars: " << p_new_comp->nVars() << " vars: ";
        all_vars_in_comp(*p_new_comp, v) cout << *v << " ";
        cout << endl;
#endif
        delete p_new_comp;
      }
    }
  }

  debug_print("We now set the unprocessed_comps_end_ in 'top' to comp_stack.size(): "
      << comp_stack.size() << ", while top.remaining_comps_ofs(): " << top.remaining_comps_ofs());
  top.set_unprocessed_comps_end(comp_stack.size());
  sortCompStackRange(new_comps_start_ofs, comp_stack.size());
}
