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

class ComponentManager
{
public:
  ComponentManager(const SolverConfiguration &config, DataAndStatistics &statistics,
                   const LiteralIndexedVector<TriValue> &lit_values,
                   const set<unsigned> &independent_support_) :
      config_(config), statistics_(statistics), cache_(statistics, config_),
      ana_(lit_values, independent_support_)
  {
  }

  ~ComponentManager() {
    for (const auto& c: seedforCLHASH) {
      assert(c != NULL);
      free(c);
    }
  }

  void initialize(LiteralIndexedVector<LitWatchList> &literals,
                  vector<Lit> &lit_pool, unsigned num_variables);

  unsigned scoreOf(VariableIndex v)
  {
    return ana_.scoreOf(v);
  }

  void cacheModelCountOf(unsigned stack_comp_id, const mpz_class &value)
  {
    if (config_.perform_comp_caching)
      cache_.storeValueOf(comp_stack_[stack_comp_id]->id(), value);
  }

  Component &getSuperComponentOf(const StackLevel &lev)
  {
    assert(comp_stack_.size() > lev.super_comp());
    return *comp_stack_[lev.super_comp()];
  }

  unsigned comp_stack_size()
  {
    return comp_stack_.size();
  }

  const Component* at(const size_t at) const {
    return comp_stack_.at(at);
  }

  void cleanRemainingComponentsOf(StackLevel &top)
  {
    print_debug(COLYEL2 "cleaning remaining comps of var: " << top.getbranchvar());
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

  Component &currentRemainingComponentOf(StackLevel &top)
  {
    assert(comp_stack_.size() > top.currentRemainingComponent());
    return *comp_stack_[top.currentRemainingComponent()];
  }

  // checks for the next yet to explore remaining comp of top
  // returns true if a non-trivial non-cached comp
  // has been found and is now stack_.TOS_NextComp()
  // returns false if all comps have been processed
  inline bool findNextRemainingComponentOf(StackLevel &top);
  inline void recordRemainingCompsFor(StackLevel &top);
  inline void sortComponentStackRange(unsigned start, unsigned end);
  inline float cacheScoreOf(VariableIndex v);
  inline void increasecachescores();
  inline void decreasecachescore(Component &comp);

  void gatherStatistics()
  {
    cache_.compute_size_used();
  }

  void removeAllCachePollutionsOf(StackLevel &top);
  vector<void *> seedforCLHASH;
  void getrandomseedforclhash()
  {
    std::mt19937_64 eng(config_.randomseed); //Use the 64-bit Mersenne Twister 19937 generator
                               //and seed it with entropy.
    std::uniform_int_distribution<unsigned long long> distr;
    assert(seedforCLHASH.empty());
    seedforCLHASH.resize(config_.hashrange);
    for (unsigned i = 0; i < config_.hashrange; i++) {
      seedforCLHASH[i] = get_random_key_for_clhash(distr(eng), distr(eng));
    }
  }

private:
  const SolverConfiguration &config_;
  DataAndStatistics &statistics_;

  vector<Component *> comp_stack_;
  ComponentCache cache_;
  ComponentAnalyzer ana_;
  vector<float> cachescore_;
};

float ComponentManager::cacheScoreOf(VariableIndex v)
{
  return cachescore_[v];
}

void ComponentManager::increasecachescores()
{
  for (unsigned i = 0; i < cachescore_.size(); i++)
  {
    cachescore_[i] *= 0.5;
  }
}
void ComponentManager::decreasecachescore(Component &comp)
{
  for (vector<VariableIndex>::const_iterator it = comp.varsBegin();
       *it != varsSENTINEL; it++)
  {
    cachescore_[*it] -= 1;
  }
}

void ComponentManager::sortComponentStackRange(unsigned start, unsigned end)
{
  print_debug(COLYEL2 "sorting comp stack range");
  assert(start <= end);
  // sort the remaining comps for processing
  for (unsigned i = start; i < end; i++)
    for (unsigned j = i + 1; j < end; j++)
    {
      if (comp_stack_[i]->num_variables() < comp_stack_[j]->num_variables())
        std::swap(comp_stack_[i], comp_stack_[j]);
    }
}

bool ComponentManager::findNextRemainingComponentOf(StackLevel &top)
{
  print_debug(COLREDBG"-*-> Running findNextRemainingComponentOf");
  print_debug("top.remaining_comps_ofs():" << top.remaining_comps_ofs() );
  if (comp_stack_.size() <= top.remaining_comps_ofs())
    recordRemainingCompsFor(top);

  assert(!top.branch_found_unsat());
  if (top.hasUnprocessedComponents()) {
    print_debug(COLREDBG"-*-> Finished findNextRemainingComponentOf, hasUnprocessedComponents.");
    return true;
  }

  // if no comp remains
  // make sure, at least that the current branch is considered SAT

  top.includeSolution(1);
  print_debug(COLREDBG "-*-> Finished findNextRemainingComponentOf, no more remaining comps. top.branchvar() was: " << top.getbranchvar()  <<" includeSolution(1) fired, returning.");
  return false;
}

// This creates comps
void ComponentManager::recordRemainingCompsFor(StackLevel &top)
{
  const Component& super_comp = getSuperComponentOf(top);
  const unsigned new_comps_start_ofs = comp_stack_.size();

  ana_.setupAnalysisContext(top, super_comp);

  for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++) {
    print_debug("checking var: " << *vt << " which comp it's in");
    if (ana_.isUnseenAndActive(*vt) && ana_.exploreRemainingCompOf(*vt)) {
      // Create new comp
      Component *p_new_comp = ana_.makeComponentFromArcheType();
      CacheableComponent *packed_comp = NULL;
      if (config_.perform_pcc) {
#ifdef DOPCC
        packed_comp = new CacheableComponent(seedforCLHASH, ana_.getArchetype().current_comp_for_caching_);
#else
        assert(false);
        exit(-1);
#endif
      } else {
        packed_comp = new CacheableComponent(ana_.getArchetype().current_comp_for_caching_);
      }

      // Check if new comp is already in cache
      if (!cache_.manageNewComponent(top, *packed_comp)) {
        comp_stack_.push_back(p_new_comp);
        p_new_comp->set_id(cache_.storeAsEntry(*packed_comp, super_comp.id()));
        cout << COLYEL2 "New comp. ID: " << p_new_comp->id()
            << " num vars: " << p_new_comp->num_variables() << " vars: ";
        for(auto v = p_new_comp->varsBegin(); *v != varsSENTINEL; v++) cout << *v << " ";
        cout << endl;
      } else {
        cout << COLYEL2 "Component already in cache. ID: " << p_new_comp->id()
            << " num vars: " << p_new_comp->num_variables() << " vars: ";
        for(auto v = p_new_comp->varsBegin(); *v != varsSENTINEL; v++) cout << *v << " ";
        cout << endl;

        //cache score should be decreased since we have a cache hit
        if (config_.use_csvsads) {
          statistics_.numcachedec_++;
          if (statistics_.numcachedec_ % 128 == 0) increasecachescores();
          for (auto it = p_new_comp->varsBegin(); *it != varsSENTINEL; it++)
            cachescore_[*it] -= 1;
        }
        delete packed_comp;
        delete p_new_comp;
      }
    }
  }
  top.set_unprocessed_comps_end(comp_stack_.size());
  sortComponentStackRange(new_comps_start_ofs, comp_stack_.size());
}

#endif /* COMPONENT_MANAGEMENT_H_ */
