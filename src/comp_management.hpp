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

#include "comp_types/comp.hpp"
#include "comp_cache.hpp"
#include "comp_analyzer.hpp"

#include <random>
#include <gmpxx.h>
#include "containers.hpp"
#include "stack.hpp"
#include "clhash/clhash.h"
#include "counter_config.hpp"

template<typename T> class Counter;

// There is exactly ONE of this, inside Counter
template<typename T>
class CompManager {
public:

  CompManager(const CounterConfiguration &config, DataAndStatistics<T> &statistics,
                   const LiteralIndexedVector<TriValue> &lit_values,
                   const uint32_t& indep_support_end, Counter<T>* _solver);
  ~CompManager() {
    free(hash_seed);
    for(auto& comp: comp_stack) delete comp;
    comp_stack.clear();
  }

  auto freq_score_of(uint32_t v) const { return ana.freq_score_of(v); }
  void initialize(const LiteralIndexedVector<LitWatchList> &watches,
    const ClauseAllocator<T>* _alloc, const vector<ClauseOfs>& long_irred_cls);
  const CompCache<T>& get_cache() const { return cache; }
  const CompAnalyzer<T>& get_ana() const { return ana; }

  uint64_t get_num_cache_entries_used() const { return cache.get_num_entries_used(); }
  void save_count(const uint32_t stack_comp_id, const T& value) {
    debug_print(COLYEL2 << "Store. comp ID: " << stack_comp_id
        << " cache ID: " << comp_stack[stack_comp_id]->id() << " cnt: " << value);
    cache.store_value(comp_stack[stack_comp_id]->id(), value);
  }

  const auto& get_comp_stack() const { return comp_stack; }
  Comp& get_super_comp(const StackLevel<T>& lev) {
    assert(comp_stack.size() > lev.super_comp());
    return *comp_stack[lev.super_comp()];
  }

  uint32_t comp_stack_size() { return comp_stack.size(); }
  const Comp* at(const size_t at) const { return comp_stack.at(at); }
  void clean_remain_comps_of(const StackLevel<T>& top) {
    debug_print(COLYEL2 "cleaning (all remaining) comps of var: " << top.var);
    while (comp_stack.size() > top.remaining_comps_ofs()) {
      if (cache.exists(comp_stack.back()->id())) cache.entry(comp_stack.back()->id()).set_deletable();

      debug_print(COLYEL2 "-> deleting comp ID: " << comp_stack.back()->id());
      delete comp_stack.back();
      comp_stack.pop_back();
    }
    assert(top.remaining_comps_ofs() <= comp_stack.size());
  }

  // checks for the next yet to explore remaining comp of top
  // returns true if a non-trivial non-cached comp
  // has been found and is now stack_.TOS_NextComp()
  // returns false if all comps have been processed
  inline bool find_next_remain_comp_of(StackLevel<T> &top);
  void record_remaining_comps_for(StackLevel<T> &top);
  inline void sort_comp_stack_range(uint32_t start, uint32_t end);
  inline double get_alternate_score_comps(uint32_t start, uint32_t end) const;

  void removeAllCachePollutionsOfIfExists(const StackLevel<T> &top);
  void removeAllCachePollutionsOf(const StackLevel<T> &top);
  void* hash_seed; //stores a bunch of __m128 aligned data pieces, each
                   //133*8 long, see: RANDOM_BYTES_NEEDED_FOR_CLHASH
  void getrandomseedforclhash() {
    std::mt19937_64 eng(conf.seed);
    std::uniform_int_distribution<uint64_t> distr;
    hash_seed = get_random_key_for_clhash(distr(eng), distr(eng));
  }

private:
  static constexpr bool weighted = std::is_same<T, mpfr::mpreal>::value || std::is_same<T, mpq_class>::value;
  const CounterConfiguration &conf;
  DataAndStatistics<T> &stats;

  // components thus far found. There is one at pos 0 that's DUMMY (empty!)
  vector<Comp*> comp_stack;
  CompCache<T> cache;
  CompAnalyzer<T> ana;
  Counter<T>* counter;
};

template<typename T>
void CompManager<T>::sort_comp_stack_range(uint32_t start, uint32_t end) {
  debug_print(COLYEL2 "sorting comp stack range");
  assert(start <= end);
  if (start == end) return;
  // sort the remaining comps for processing
  stats.comp_sorts++;
  stats.comp_sizes+= end - start;
  for (uint32_t i = start; i < end; i++)
    for (uint32_t j = i + 1; j < end; j++) {
      if (comp_stack[i]->nVars()
                  < comp_stack[j]->nVars())
        std::swap(comp_stack[i], comp_stack[j]);
    }
}

template<typename T>
double CompManager<T>::get_alternate_score_comps(uint32_t start, uint32_t end) const
{
  double score = 1;
  assert(start <= end);
  // sort the remaining comps for processing
  for (uint32_t i = start; i < end; i++) score *= comp_stack[i]->nVars();
  return score;
}

template<typename T>
bool CompManager<T>::find_next_remain_comp_of(StackLevel<T>& top)
{
  debug_print(COLREDBG"-*-> Running find_next_remain_comp_of");
  debug_print("top.remaining_comps_ofs():" << top.remaining_comps_ofs()
      << " comp_stack.size(): " << comp_stack.size());
  if (comp_stack.size() <= top.remaining_comps_ofs()) record_remaining_comps_for(top);
  else {
    debug_print("Not running record_remaining_comps_for, comp_stack.size() > top.remaining_comps_ofs()."
        " comp_stack.size(): " << comp_stack.size()
        << " top.reimaining_comps_ofs(): " << top.remaining_comps_ofs());
  }

  if (top.branch_found_unsat()) {
    // can ONLY happen with negative weights
    if constexpr (!weighted) assert(false);
    return false;
  }
  if (top.has_unproc_comps()) {
    debug_print(COLREDBG"-*-> Finished find_next_remain_comp_of, has_unproc_comps.");
    return true;
  }

  // if no component remains
  // make sure, at least that the current branch is considered SAT
  top.include_solution(1);
  debug_print(COLREDBG "-*-> Finished find_next_remain_comp_of, no more remaining comps. "
      "top.branchvar() was: "
      << top.var  <<" include_solution(1) fired. "
      " New Model cnt: " << top.total_model_count()
      << " left: " << top.left_model_count() << " right: "
      << top.right_model_count() << " , returning.");
  return false;
}

template<typename T>
CompManager<T>::CompManager(const CounterConfiguration& config, DataAndStatistics<T>& statistics,
                 const LiteralIndexedVector<TriValue>& lit_values,
                 const uint32_t& indep_support_end, Counter<T>* _counter) :
    conf(config), stats(statistics), cache(_counter->nVars(), statistics, conf),
    ana(lit_values, indep_support_end, _counter), counter(_counter)
{ }


// Initialized exactly once when Counter is created.
//   it also inits the included analyzer called "ana"
template<typename T>
void CompManager<T>::initialize(const LiteralIndexedVector<LitWatchList> & watches,
    const ClauseAllocator<T>* _alloc, const vector<ClauseOfs>& long_irred_cls){
  assert(comp_stack.empty());

  ana.initialize(watches, _alloc, long_irred_cls);

  //Add dummy comp
  comp_stack.push_back(new Comp());

  //Add full comp
  comp_stack.push_back(new Comp());
  assert(comp_stack.size() == 2);
  comp_stack.back()->create_init_comp(ana.get_max_var() , ana.get_max_clid());
  cache.init(*comp_stack.back(), hash_seed);
}

