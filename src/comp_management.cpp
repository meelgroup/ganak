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
#include <gmpxx.h>

using namespace GanakInt;

CompManager::CompManager(const CounterConfiguration& config,
    DataAndStatistics& statistics,
    const LiteralIndexedVector<TriValue>& lit_values, Counter* _counter) :
    fg(_counter->get_fg()->dup()), conf(config), stats(statistics), cache(statistics, conf),
    ana(lit_values, _counter)
{
  get_random_seed_for_hash();
}

void CompManager::remove_cache_pollutions_of_if_exists(const StackLevel &top) {
  assert(top.remaining_comps_ofs() <= comp_stack.size());
  assert(top.super_comp() != 0);

  for (uint32_t u = top.remaining_comps_ofs(); u < comp_stack.size(); u++) {
    if (!cache.exists(comp_stack[u]->id())) continue;
    stats.cache_pollutions_removed += cache.clean_pollutions_involving(comp_stack[u]->id());
  }
  stats.cache_pollutions_called++;
}

void CompManager::remove_cache_pollutions_of(const StackLevel &top) {
  // all processed comps are found in
  // [top.curr_remain_comp(), comp_stack.size())
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
// Runs a LOT of the time, like 40%+
void CompManager::record_remaining_comps_for(StackLevel &top) {
  const Comp& super_comp = get_super_comp(top);
  const uint32_t new_comps_start_ofs = comp_stack.size();

  // This reinitializes archetype, sets up seen[] or all cls&vars unvisited (if unset), etc.
  // Sets all unknown vars in seen[] and sets all clauses in seen[] to unvisited
  // Also zeroes out frequency_scores. Sets num_long_cls and num_bin_cls to 0
  ana.setup_analysis_context(top, super_comp);

  all_vars_in_comp(super_comp, vt) {
    debug_print("Going to NEXT var that's unvisited & set in this component... if it exists. Var: " << *vt);
    if (ana.var_unvisited_in_sup_comp(*vt) &&
        ana.explore_comp(*vt, super_comp.num_long_cls(), super_comp.num_bin_cls())) {
      Comp *p_new_comp = ana.make_comp_from_archetype();
      CacheableComp packed_comp(hash_seed, *p_new_comp, ana.get_long_clauses_data_ptrs(),
          ana.get_vals());

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
        free(p_new_comp);
      }
    }
  }

  debug_print("We now set the unprocessed_comps_end_ in 'top' to comp_stack.size(): "
      << comp_stack.size() << ", while top.remaining_comps_ofs(): " << top.remaining_comps_ofs());
  top.set_unprocessed_comps_end(comp_stack.size());
  sort_comp_stack_range(new_comps_start_ofs, comp_stack.size());
}
