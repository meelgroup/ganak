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

#include "comp_types/base_packed_comp.h"
#include "comp_types/comp.h"
#include "comp_cache.h"
#include "alt_comp_analyzer.h"

#include <unordered_map>
#include <random>
#include <gmpxx.h>
#include "containers.h"
#include "stack.h"
#include "clhash/clhash.h"
/* #include "clhash/minim.h" */
#include "counter_config.h"

class Counter;

// There is exactly ONE of this, inside Counter
class ComponentManager
{
public:
  ComponentManager(const CounterConfiguration &config, DataAndStatistics &statistics,
                   const LiteralIndexedVector<TriValue> &lit_values,
                   const uint32_t& indep_support_end, Counter* _solver) :
      config_(config), stats(statistics), cache_(statistics, config_, sz),
      ana_(lit_values, indep_support_end, _solver), solver_(_solver)
  { }

  ~ComponentManager() {
    free(randomseedforCLHASH);
    for(auto& comp: comp_stack_) delete comp;
    comp_stack_.clear();
  }

  unsigned& scoreOf(VariableIndex v) { return ana_.scoreOf(v); }
  unsigned scoreOf(VariableIndex v) const { return ana_.scoreOf(v); }

  void initialize(LiteralIndexedVector<LitWatchList> &watches_,
    const ClauseAllocator* _alloc, const vector<ClauseOfs>& longIrredCls, uint32_t nVars);
  void delete_comps_with_vars(const set<uint32_t>& vars) {
    cache_.delete_comps_with_vars(vars);
  }
  const ComponentCache& get_cache() const { return cache_; }


  uint64_t get_num_cache_entries_used() const { return cache_.get_num_entries_used(); }
  void cacheModelCountOf(uint32_t stack_comp_id, const mpz_class &value) {
    cache_.storeValueOf(comp_stack_[stack_comp_id]->id(), value);
  }

  Component& getSuperComponentOf(const StackLevel &lev) {
    assert(comp_stack_.size() > lev.super_comp());
    return *comp_stack_[lev.super_comp()];
  }

  uint32_t comp_stack_size() { return comp_stack_.size(); }
  const Component* at(const size_t at) const { return comp_stack_.at(at); }
  double cacheScoreOf(const VariableIndex v) const { return cachescore_[v]; }
  void cleanRemainingComponentsOf(const StackLevel &top)
  {
    print_debug(COLYEL2 "cleaning (all remaining) comps of var: " << top.var);
    while (comp_stack_.size() > top.remaining_comps_ofs())
    {
      if (cache_.hasEntry(comp_stack_.back()->id()))
        cache_.entry(comp_stack_.back()->id()).set_deletable();

      print_debug(COLYEL2 "-> deleting comp ID: " << comp_stack_.back()->id());
      delete comp_stack_.back();
      comp_stack_.pop_back();
    }
    assert(top.remaining_comps_ofs() <= comp_stack_.size());
  }

  // checks for the next yet to explore remaining comp of top
  // returns true if a non-trivial non-cached comp
  // has been found and is now stack_.TOS_NextComp()
  // returns false if all comps have been processed
  inline bool findNextRemainingComponentOf(StackLevel &top);
  void recordRemainingCompsFor(StackLevel &top);
  inline void sortComponentStackRange(uint32_t start, uint32_t end);
  inline double get_alternate_score_comps(uint32_t start, uint32_t end) const;

  void removeAllCachePollutionsOf(const StackLevel &top);
  void* randomseedforCLHASH; //stores a bunch of __m128 aligned data pieces, each
                                //133*8 long, see: RANDOM_BYTES_NEEDED_FOR_CLHASH
  void getrandomseedforclhash()
  {
    std::mt19937_64 eng(config_.seed); //Use the 64-bit Mersenne Twister 19937 generator
                               //and seed it with entropy.
    std::uniform_int_distribution<uint64_t> distr;
    randomseedforCLHASH = get_random_key_for_clhash(distr(eng), distr(eng));
  }

  void rescale_cache_scores() { for (auto& c: cachescore_) c *= 0.5; }
  void decreasecachescore(Component &comp) {
    for (vector<VariableIndex>::const_iterator it = comp.varsBegin();
         *it != varsSENTINEL; it++) {
      cachescore_[*it] -= 1;
    }
  }

private:
  const CounterConfiguration &config_;
  DataAndStatistics &stats;

  // components thus far found. There is one at pos 0 that's DUMMY (empty!)
  vector<Component *> comp_stack_;
  ComponentCache cache_;
  ComponentAnalyzer ana_;
  vector<double> cachescore_;
  Counter* solver_;
  BPCSizes sz;
  vector<uint32_t> tmp_data_for_pcc;
};

void ComponentManager::sortComponentStackRange(uint32_t start, uint32_t end)
{
  print_debug(COLYEL2 "sorting comp stack range");
  assert(start <= end);
  // sort the remaining comps for processing
  for (uint32_t i = start; i < end; i++)
    for (uint32_t j = i + 1; j < end; j++)
    {
      if (comp_stack_[i]->nVars() < comp_stack_[j]->nVars())
        std::swap(comp_stack_[i], comp_stack_[j]);
    }
}

double ComponentManager::get_alternate_score_comps(uint32_t start, uint32_t end) const
{
  double score = 1;
  assert(start <= end);
  // sort the remaining comps for processing
  for (uint32_t i = start; i < end; i++) score *= comp_stack_[i]->nVars();
  return score;
}


bool ComponentManager::findNextRemainingComponentOf(StackLevel &top)
{
  print_debug(COLREDBG"-*-> Running findNextRemainingComponentOf");
  print_debug("top.remaining_comps_ofs():" << top.remaining_comps_ofs() << " comp_stack_.size(): " << comp_stack_.size());
  if (comp_stack_.size() <= top.remaining_comps_ofs()) {
    recordRemainingCompsFor(top);
  } else {
    print_debug("Not running recordRemainingCompsFor, comp_stack_.size() > top.remaining_comps_ofs(). comp_stack_.size(): " << comp_stack_.size() << " top.reimaining_comps_ofs(): " << top.remaining_comps_ofs());
  }

  assert(!top.branch_found_unsat());
  if (top.hasUnprocessedComponents()) {
    print_debug(COLREDBG"-*-> Finished findNextRemainingComponentOf, hasUnprocessedComponents.");
    return true;
  }

  // if no comp remains then there is exactly 1 solution left
  top.includeSolution(1);
  print_debug(COLREDBG "-*-> Finished findNextRemainingComponentOf, no more remaining comps. top.branchvar() was: " << top.var  <<" includeSolution(1) fired, returning.");
  return false;
}
