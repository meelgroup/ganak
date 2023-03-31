/*
 * comp_management.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#ifndef COMPONENT_MANAGEMENT_H_
#define COMPONENT_MANAGEMENT_H_

#include "comp_types/comp.h"
#include "comp_cache.h"
#include "alt_comp_analyzer.h"

#include <unordered_map>
#include <random>
#include <gmpxx.h>
#include "containers.h"
#include "stack.h"
#ifdef DOPCC
#include "clhash/clhash.h"
#else
#include "clhash/minim.h"
#endif
#include "solver_config.h"

// There is exactly ONE of this, inside Solver
class ComponentManager
{
public:
  ComponentManager(const SolverConfiguration &config, DataAndStatistics &statistics,
                   const LiteralIndexedVector<TriValue> &lit_values,
                   const set<uint32_t> &indep_support_) :
      config_(config), stats(statistics), cache_(statistics, config_),
      ana_(lit_values, indep_support_)
  {
  }

  ~ComponentManager() {
    for (const auto& c: seedforCLHASH) {
      assert(c != NULL);
      free(c);
    }
  }

  void initialize(LiteralIndexedVector<LitWatchList> &literals, vector<Lit> &lit_pool);
  const ComponentCache& get_cache() const { return cache_; }


  uint64_t get_num_cache_entries_used() const
  {
    return cache_.get_num_entries_used();
  }

  void cacheModelCountOf(uint32_t stack_comp_id, const mpz_class &value)
  {
    if (config_.do_comp_caching)
      cache_.storeValueOf(comp_stack_[stack_comp_id]->id(), value);
  }

  Component& getSuperComponentOf(const StackLevel &lev)
  {
    assert(comp_stack_.size() > lev.super_comp());
    return *comp_stack_[lev.super_comp()];
  }

  uint32_t comp_stack_size()
  {
    return comp_stack_.size();
  }

  const Component* at(const size_t at) const {
    return comp_stack_.at(at);
  }

  void cleanRemainingComponentsOf(const StackLevel &top)
  {
    print_debug(COLYEL2 "cleaning (all remaining) comps of var: " << top.getbranchvar());
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
  inline void recordRemainingCompsFor(StackLevel &top);
  inline void sortComponentStackRange(uint32_t start, uint32_t end);

  void gatherStatistics()
  {
    cache_.compute_size_used();
  }

  void removeAllCachePollutionsOf(const StackLevel &top);
  vector<void *> seedforCLHASH; //stores a bunch of __m128 aligned data pieces, each
                                //133*8 long, see: RANDOM_BYTES_NEEDED_FOR_CLHASH
  void getrandomseedforclhash()
  {
    std::mt19937_64 eng(config_.randomseed); //Use the 64-bit Mersenne Twister 19937 generator
                               //and seed it with entropy.
    std::uniform_int_distribution<uint64_t> distr;
    assert(seedforCLHASH.empty());
    seedforCLHASH.resize(config_.hashrange);
    for (uint32_t i = 0; i < config_.hashrange; i++) {
      seedforCLHASH[i] = get_random_key_for_clhash(distr(eng), distr(eng));
    }
  }

private:
  const SolverConfiguration &config_;
  DataAndStatistics &stats;

  // components thus far found. There is one at pos 0 that's DUMMY (empty!)
  vector<Component *> comp_stack_;
  ComponentCache cache_;
  ComponentAnalyzer ana_;
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
  print_debug(COLREDBG "-*-> Finished findNextRemainingComponentOf, no more remaining comps. top.branchvar() was: " << top.getbranchvar()  <<" includeSolution(1) fired, returning.");
  return false;
}

// This creates comps
void ComponentManager::recordRemainingCompsFor(StackLevel &top)
{
  const Component& super_comp = getSuperComponentOf(top);
  const uint32_t new_comps_start_ofs = comp_stack_.size();

  // This reinitializes archetype, sets up seen[] or all cls&vars unseen (if unset), etc.
  ana_.setupAnalysisContext(top, super_comp);

  for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++) {
    print_debug("Going to NEXT var that's unseen & active in this component... if it exists. Var: " << *vt);
    if (ana_.isUnseenAndActive(*vt) && ana_.exploreRemainingCompOf(*vt)) {
      Component *p_new_comp = ana_.makeComponentFromArcheType();
      CacheableComponent *packed_comp = NULL;
      if (config_.do_pcc) {
#ifdef DOPCC
        packed_comp = new CacheableComponent(seedforCLHASH, ana_.getArchetype().current_comp_for_caching_);
        packed_comp->finish_hashing(packed_comp->SizeInBytes(), packed_comp->nVars());
#else
        exit(-1);
#endif
      } else {
        packed_comp = new CacheableComponent(ana_.getArchetype().current_comp_for_caching_);
      }

      // Check if new comp is already in cache
      if (!cache_.manageNewComponent(top, *packed_comp)) {
        comp_stack_.push_back(p_new_comp);
        p_new_comp->set_id(cache_.storeAsEntry(*packed_comp, super_comp.id()));
#ifdef VERBOSE_DEBUG
        cout << COLYEL2 "New comp. ID: " << p_new_comp->id()
            << " num vars: " << p_new_comp->nVars() << " vars: ";
        for(auto v = p_new_comp->varsBegin(); *v != varsSENTINEL; v++) cout << *v << " ";
        cout << endl;
#endif
      } else {
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

#endif /* COMPONENT_MANAGEMENT_H_ */
