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

#include "comp_types/base_packed_comp.hpp"
#include "comp_types/comp.hpp"
#include "comp_cache.hpp"
#include "comp_analyzer.hpp"

#include <unordered_map>
#include <random>
#include <gmpxx.h>
#include "containers.hpp"
#include "stack.hpp"
#include "clhash/clhash.h"
/* #include "clhash/minim.hpp" */
#include "counter_config.hpp"

class Counter;

// There is exactly ONE of this, inside Counter
class CompManager
{
public:
  CompManager(const CounterConfiguration &config, DataAndStatistics &statistics,
                   const LiteralIndexedVector<TriValue> &lit_values,
                   const uint32_t& indep_support_end, Counter* _solver) :
      conf(config), stats(statistics), cache(statistics, conf),
      ana(lit_values, indep_support_end, _solver), solver_(_solver)
  { }

  ~CompManager() {
    free(hash_seed);
    for(auto& comp: comp_stack_) delete comp;
    comp_stack_.clear();
  }

  double freq_score_of(uint32_t v) const { return ana.freq_score_of(v); }

  void initialize(const LiteralIndexedVector<LitWatchList> &watches,
    const ClauseAllocator* _alloc, const vector<ClauseOfs>& long_irred_cls, uint32_t nVars);
  void delete_comps_with_vars(const set<uint32_t>& vars) {
    cache.delete_comps_with_vars(vars);
  }
  const CompCache& get_cache() const { return cache; }
  const CompAnalyzer& get_ana() const { return ana; }

  uint64_t get_num_cache_entries_used() const { return cache.get_num_entries_used(); }
  void save_count(uint32_t stack_comp_id, const mpz_class &value) {
    cache.storeValueOf(comp_stack_[stack_comp_id]->id(), value);
  }

  Comp& getSuperCompOf(const StackLevel &lev) {
    assert(comp_stack_.size() > lev.super_comp());
    return *comp_stack_[lev.super_comp()];
  }

  uint32_t comp_stack_size() { return comp_stack_.size(); }
  const Comp* at(const size_t at) const { return comp_stack_.at(at); }
  void cleanRemainingCompsOf(const StackLevel &top)
  {
    debug_print(COLYEL2 "cleaning (all remaining) comps of var: " << top.var);
    while (comp_stack_.size() > top.remaining_comps_ofs())
    {
      if (cache.hasEntry(comp_stack_.back()->id()))
        cache.entry(comp_stack_.back()->id()).set_deletable();

      debug_print(COLYEL2 "-> deleting comp ID: " << comp_stack_.back()->id());
      delete comp_stack_.back();
      comp_stack_.pop_back();
    }
    assert(top.remaining_comps_ofs() <= comp_stack_.size());
  }

  // checks for the next yet to explore remaining comp of top
  // returns true if a non-trivial non-cached comp
  // has been found and is now stack_.TOS_NextComp()
  // returns false if all comps have been processed
  inline bool findNextRemainingCompOf(StackLevel &top);
  void recordRemainingCompsFor(StackLevel &top);
  inline void sortCompStackRange(uint32_t start, uint32_t end);
  inline double get_alternate_score_comps(uint32_t start, uint32_t end) const;

  void removeAllCachePollutionsOfIfExists(const StackLevel &top);
  void removeAllCachePollutionsOf(const StackLevel &top);
  void* hash_seed; //stores a bunch of __m128 aligned data pieces, each
                                //133*8 long, see: RANDOM_BYTES_NEEDED_FOR_CLHASH
  void getrandomseedforclhash()
  {
    std::mt19937_64 eng(conf.seed); //Use the 64-bit Mersenne Twister 19937 generator
                               //and seed it with entropy.
    std::uniform_int_distribution<uint64_t> distr;
    hash_seed = get_random_key_for_clhash(distr(eng), distr(eng));
  }

  double get_cache_hit_score(const uint32_t v) const { return cache_hit_score[v]; }
  void bump_cache_hit_score(Comp &comp) {
    for (vector<uint32_t>::const_iterator it = comp.vars_begin(); *it != sentinel; it++) {
      cache_hit_score[*it] += act_inc;
      if (cache_hit_score[*it] > 1e100) {
        for (auto& s: cache_hit_score) s *= 1e-90;
        act_inc *= 1e-90;
      }
    }
  }

private:
  const CounterConfiguration &conf;
  DataAndStatistics &stats;

  // components thus far found. There is one at pos 0 that's DUMMY (empty!)
  vector<Comp *> comp_stack_;
  CompCache cache;
  CompAnalyzer ana;

  // indexed by variable, decremented when a variable is in a component,
  // and halved once in a while. The LARGER it is, the more likely the
  // variable gets picked for branching. So basically, the fewer times a
  // variable is in a component, the more likely the branch
  vector<double> cache_hit_score;
  double act_inc = 1.0;
  Counter* solver_;
};

void CompManager::sortCompStackRange(uint32_t start, uint32_t end) {
  debug_print(COLYEL2 "sorting comp stack range");
  assert(start <= end);
  // sort the remaining comps for processing
  for (uint32_t i = start; i < end; i++)
    for (uint32_t j = i + 1; j < end; j++) {
      if (conf.do_comp_reverse_sort && comp_stack_[i]->nVars() > comp_stack_[j]->nVars())
        std::swap(comp_stack_[i], comp_stack_[j]);
      if (!conf.do_comp_reverse_sort && comp_stack_[i]->nVars() < comp_stack_[j]->nVars())
        std::swap(comp_stack_[i], comp_stack_[j]);
    }
}

double CompManager::get_alternate_score_comps(uint32_t start, uint32_t end) const
{
  double score = 1;
  assert(start <= end);
  // sort the remaining comps for processing
  for (uint32_t i = start; i < end; i++) score *= comp_stack_[i]->nVars();
  return score;
}


bool CompManager::findNextRemainingCompOf(StackLevel &top)
{
  debug_print(COLREDBG"-*-> Running findNextRemainingCompOf");
  debug_print("top.remaining_comps_ofs():" << top.remaining_comps_ofs() << " comp_stack_.size(): " << comp_stack_.size());
  if (comp_stack_.size() <= top.remaining_comps_ofs()) {
    recordRemainingCompsFor(top);
  } else {
    debug_print("Not running recordRemainingCompsFor, comp_stack_.size() > top.remaining_comps_ofs(). comp_stack_.size(): " << comp_stack_.size() << " top.reimaining_comps_ofs(): " << top.remaining_comps_ofs());
  }

  assert(!top.branch_found_unsat());
  if (top.hasUnprocessedComps()) {
    debug_print(COLREDBG"-*-> Finished findNextRemainingCompOf, hasUnprocessedComps.");
    return true;
  }

  // if no comp remains then there is exactly 1 solution left
  top.includeSolution(1);
  debug_print(COLREDBG "-*-> Finished findNextRemainingCompOf, no more remaining comps. top.branchvar() was: " << top.var  <<" includeSolution(1) fired, returning.");
  return false;
}
