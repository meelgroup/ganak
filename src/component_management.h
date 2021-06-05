/*
 * component_management.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#ifndef COMPONENT_MANAGEMENT_H_
#define COMPONENT_MANAGEMENT_H_

#include "component_types/component.h"
#include "component_cache.h"
#include "alt_component_analyzer.h"
//#include "component_analyzer.h"

// #include <vector>
#include <unordered_map>
#include <random>
#include <set>
#include <gmpxx.h>
#include "containers.h"
#include "stack.h"
#include "clhash/clhash.h"
#include "solver_config.h"
#include "component_types/difference_packed_component.h"
#include "component_types/ClHashComponent.h"

#include <iostream>
#include <fstream>

using namespace std;

typedef AltComponentAnalyzer ComponentAnalyzer;

class ComponentManager
{
public:
  ComponentManager(SolverConfiguration &config, DataAndStatistics &statistics,
                   LiteralIndexedVector<TriValue> &lit_values,
                   LiteralIndexedVector<Literal> &literals,
                   std::set<unsigned> &independent_support,
                   bool perform_projected_model_counting) : config_(config),
                                                            statistics_(statistics),
                                                            cache_(statistics, config_),
                                                            ana_(statistics, lit_values, independent_support, perform_projected_model_counting, config_, literals)
  {
  }
  ~ComponentManager() {
      for (auto* ptr : component_stack_) {
          delete ptr;
      }
  }

  void initialize(
      LiteralIndexedVector<Literal> &literals,
      shared_ptr<vector<LiteralID>> literal_pool,
      unsigned num_variables);

  void destroy();

  unsigned scoreOf(VariableIndex v)
  {
    return ana_.scoreOf(v);
  }

  void cacheModelCountOf(unsigned stack_comp_id, const mpz_class &value)
  {
    if (config_.perform_component_caching)
      cache_.storeValueOf(component_stack_[stack_comp_id]->id(), value);
  }

  Component &superComponentOf(StackLevel &lev)
  {
    assert(component_stack_.size() > lev.super_component());
    return *component_stack_[lev.super_component()];
  }

  unsigned component_stack_size()
  {
    return component_stack_.size();
  }

  void cleanRemainingComponentsOf(StackLevel &top)
  {
    while (component_stack_.size() > top.remaining_components_ofs())
    {
      if (cache_.hasEntry(component_stack_.back()->id()))
        cache_.entry(component_stack_.back()->id()).set_deletable();
      delete component_stack_.back();
      component_stack_.pop_back();
    }
    assert(top.remaining_components_ofs() <= component_stack_.size());
  }

  Component &currentRemainingComponentOf(StackLevel &top)
  {
    assert(component_stack_.size() > top.currentRemainingComponent());
    return *component_stack_[top.currentRemainingComponent()];
  }

  // checks for the next yet to explore remaining component of top
  // returns true if a non-trivial non-cached component
  // has been found and is now stack_.TOS_NextComp()
  // returns false if all components have been processed;
  inline bool findNextRemainingComponentOf(StackLevel &top);

  inline void recordRemainingCompsFor(StackLevel &top);

  inline void sortComponentStackRange(unsigned start, unsigned end);

  inline float cacheScoreOf(VariableIndex v);

  inline void increasecachescores();

  inline void decreasecachescore(Component &comp);

  void gatherStatistics()
  {
    cache_.compute_byte_size_infrasture();
  }

    /**
     * Write all components currently in the entry base, to a file.
     * Each component C is printed as: |C| hit_count |vars|
     * with |C| the size of the component,
     * hit_count the number of times there was a cache hit on this component,
     * |vars| the number of variables that were stored in the component (only printed if > 0).
     * @param file_name The name of the file to write to..
     */
  void writeComponentsToFile(const string &file_name)
  {
    ofstream out(file_name, ios_base::app);
    out << statistics_.input_file_ << endl;
    out << "Final components C:" << endl;
    for (auto c : cache_.get_entry_base())
    {
      if (c != nullptr && c->get_cache_hit_count() > 0)
      {
        out << "|C|: " << c->num_variables();
        out << ", \thitcs: " << c->get_cache_hit_count();
        if (!c->variables_.empty())
            out << ", \t|vars|: " << c->variables_.size();
        out << endl;
      }
    }
    out << "\n"
        << endl;
  }

  /**
     * Write all components currently in the entry base, to a file.
     * Only prints components if |C| < |Vars| and hit_count > 0.
     * Each component C is printed as: |C| hit_count |vars|
     * with |C| the size of the component,
     * hit_count the number of times there was a cache hit on this component,
     * |vars| the number of variables that were stored in the component (only printed if > 0).
     * @param file_name The name of the file to write to..
     */
  void writeDiffComponentsToFile(const string &file_name)
  {
    ofstream out(file_name, ios_base::app);
    out << statistics_.input_file_ << endl;
    out << "Final components C:" << endl;
    for (auto c : cache_.get_entry_base())
    {
      if (c != nullptr && c->get_cache_hit_count() > 0 && c->num_variables() < c->variables_.size())
      {
        out << "|C|: " << c->num_variables();
        out << ", \thitcs: " << c->get_cache_hit_count();
        out << ", \t|vars|: " << c->variables_.size() << endl;
      }
    }
    out << "\n"
        << endl;
  }

  /**
     * Write the current cache score of each variable, to a file.
     * @param file_name The name of the file to write to..
     */
  void writeCacheScoresToFile(const string &file_name)
  {
    ofstream out(file_name, ios_base::app);
    out << statistics_.input_file_ << endl;
    out << "Variable scores:" << endl;
    for (unsigned i = 0; i < cachescore_.size(); i++)
    {
      out << "CacheStore: (" << i << "): \t" << cacheScoreOf(i) << endl;
    }
  }

  void removeAllCachePollutionsOf(StackLevel &top);

  vector<void *> seedforCLHASH;

  void getrandomseedforclhash()
  {
    std::random_device rd;     //Get a random seed from the OS entropy device, or whatever
    std::mt19937_64 eng(rd()); //Use the 64-bit Mersenne Twister 19937 generator
                               //and seed it with entropy.
    std::uniform_int_distribution<unsigned long long> distr;
    seedforCLHASH.reserve(config_.hashrange);
    for (unsigned int i = 0; i < config_.hashrange; i++)
    {
      seedforCLHASH[i] =
          get_random_key_for_clhash(distr(eng), distr(eng));
    }
  }

private:
  SolverConfiguration &config_;
  DataAndStatistics &statistics_;

  vector<Component *> component_stack_;
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
  assert(start <= end);
  // sort the remaining components for processing
  for (unsigned i = start; i < end; i++)
    for (unsigned j = i + 1; j < end; j++)
    {
      if (component_stack_[i]->num_variables() < component_stack_[j]->num_variables())
        swap(component_stack_[i], component_stack_[j]);
    }
}

bool ComponentManager::findNextRemainingComponentOf(StackLevel &top)
{
  // record Remaining Components if there are none!

  if (component_stack_.size() <= top.remaining_components_ofs())
  {

    recordRemainingCompsFor(top);
  }

  assert(!top.branch_found_unsat());

  if (top.hasUnprocessedComponents())
  {

    return true;
  }
  // if no component remains
  // make sure, at least that the current branch is considered SAT

  top.includeSolution(1);

  return false;
}

void ComponentManager::recordRemainingCompsFor(StackLevel &top)
{
  Component &super_comp = superComponentOf(top);
  unsigned new_comps_start_ofs = component_stack_.size();

  ana_.setupAnalysisContext(top, super_comp);

  for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++)
  {
    if (ana_.isUnseenAndActive(*vt) && ana_.exploreRemainingCompOf(*vt))
    {
      Component *p_new_comp = ana_.makeComponentFromArcheType();
      CacheableComponent *packed_comp = NULL;
      if (config_.perform_pcc)
      {
        if (config_.use_isocc)
        {
          if (p_new_comp->num_variables() >= config_.isocc_lb &&
              (!config_.isocc_ub_set || p_new_comp->num_variables() <= config_.isocc_ub))
          {
            auto temp_comp = ana_.makeIsoCompFromArcheType(p_new_comp);
            packed_comp = new ClHashComponent(seedforCLHASH, *temp_comp);
            delete temp_comp;
          }
          else
          {
            auto temp_comp = new DifferencePackedComponent(ana_.getArchetype().current_comp_for_caching_,
                                                           config_.use_icsvsads);
            packed_comp = new ClHashComponent(seedforCLHASH, *temp_comp);
            delete temp_comp;
          }
        }
        else if (config_.use_std)
        {
          //Directly create CL, skipping overhead of first constructing arrays of the standard encoding.
          packed_comp = ana_.makeSTDClHashFromArcheType(seedforCLHASH, p_new_comp);
        }
        else
        {
          auto temp_comp = new DifferencePackedComponent(ana_.getArchetype().current_comp_for_caching_);
          packed_comp = new ClHashComponent(seedforCLHASH, *temp_comp);
          delete temp_comp;
        }
      }
      else
      {
        if (config_.use_isocc)
        {
          if (p_new_comp->num_variables() >= config_.isocc_lb &&
              (!config_.isocc_ub_set || p_new_comp->num_variables() <= config_.isocc_ub))
          {
            packed_comp = ana_.makeIsoCompFromArcheType(p_new_comp);
          }
          else
          {
            packed_comp = new DifferencePackedComponent(ana_.getArchetype().current_comp_for_caching_,
                                                        config_.use_icsvsads);
          }
        }
        else if (config_.use_std)
        {
          std::cout << "STD without PCC is not implemented. It would be worse performance than Hybrid.";
          exit(EXIT_FAILURE);
          //                  packed_comp = new DifferencePackedComponent(
          //                          ana_.getArchetype().current_comp_for_caching_.getSTD());
        }
        else
        {
          packed_comp = new DifferencePackedComponent(
              ana_.getArchetype().current_comp_for_caching_);
        }
      }
      ana_.getArchetype().setComponentVarsSeen(p_new_comp);

      // Check cache
      auto cached_hit_comp = cache_.manageNewComponent(top, *packed_comp);
      if (cached_hit_comp == nullptr)
      {
        component_stack_.push_back(p_new_comp);
        p_new_comp->set_id(cache_.storeAsEntry(*packed_comp, super_comp.id()));
#ifdef VERB
        cout << "We have a cache miss for component: " << p_new_comp->id() << endl;
        p_new_comp->printcomp();
#endif
      }
      else
      {
#ifdef VERB
        cout << "We have a cache hit for component." << endl;
        p_new_comp->printcomp();
#endif
        //cache score should be decreased since we have a cache hit
        if (config_.use_csvsads)
        {
          statistics_.numcachedec_++;
          if (statistics_.numcachedec_ % 128 == 0)
          {
            increasecachescores();
          }
          if (config_.use_icsvsads)
          {
            cached_hit_comp->add_variables_of(*packed_comp);
            for (VariableIndex v : cached_hit_comp->variables_)
            {
              cachescore_[v] -= 1;
            }
          }
          else
          {
            for (auto it = p_new_comp->varsBegin();
                 *it != varsSENTINEL; it++)
            {
              cachescore_[*it] -= 1;
            }
          }
        }
        delete packed_comp;
        delete p_new_comp;
      }
    }
  }
  top.set_unprocessed_components_end(component_stack_.size());
  sortComponentStackRange(new_comps_start_ofs, component_stack_.size());
}

#endif /* COMPONENT_MANAGEMENT_H_ */
