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

#include <vector>
#include <unordered_map>
#include <gmpxx.h>
#include "containers.h"
#include "stack.h"
#include "clhash/clhash.h"
#include "solver_config.h"
using namespace std;

typedef AltComponentAnalyzer ComponentAnalyzer;

class ComponentManager {
public:
  ComponentManager(SolverConfiguration &config, DataAndStatistics &statistics,
        LiteralIndexedVector<TriValue> & lit_values, vector<Variable> & variables) :
        config_(config), statistics_(statistics), cache_(statistics, config),
        ana_(statistics,lit_values,config_, variables), variables_(variables),partial_solution_(1),saved_partial_solution_(1) {
  }

  void initialize(LiteralIndexedVector<Literal> & literals,
        vector<LiteralID> &lit_pool, unsigned int num_variables);

  float cacheScoreOf(VariableIndex v){
    return cachescore_[v];
  }

  float scoreOf(VariableIndex v) {
        return ana_.scoreOf(v);
    // if(config_.allowactivitydecrease){
    //     cout << "[component management]" << 5.0*ana_.scoreOf(v) << " "<< 0.1*cachescore_[v]<< endl;  
    //     return 5.0*ana_.scoreOf(v) + 0.05*cachescore_[v];
    // }
    // else{
    //     cout << "[component management]" << ana_.scoreOf(v) << endl;
    //   }
  }
  void increasecachescores(){
    // cout << "Doing increase "<<endl;
    for (unsigned i = 0; i< cachescore_.size();i++){
      cachescore_[i] *= 0.5;
    }
  }
  void decreasecachescore(Component &comp){
    // cout << "Doing decrease "<< endl;
    for(vector<VariableIndex>::const_iterator it = comp.varsBegin();
     *it != varsSENTINEL; it++){
        cachescore_[*it] -= 1;
     }
  }
  void printcomp(unsigned stack_comp_id){
    component_stack_[stack_comp_id]->printHash();
    cout << "Num vars :"<< component_stack_[stack_comp_id]->num_variables()<< endl;
    cout << "bliss hash :";
     component_stack_[stack_comp_id]->printblissHash();
  }
  void cacheModelCountOf(unsigned stack_comp_id, const mpf_class &value) {
    if (config_.perform_component_caching){
      // if (component_stack_[stack_comp_id]->getuseisomorphism()){
      //   cout << "hash: "<< endl;
      //   component_stack_[stack_comp_id]->printHash();
      //   cout << "blisshash: ";
      //   component_stack_[stack_comp_id]->printblissHash();
      //   cout << "Count: "<< value<< endl;
      // }
      // else{
      //   cout << "hash: "<< endl;
      //   component_stack_[stack_comp_id]->printHash();
      //   cout << "Count: "<< value<< endl;
      // }
      cache_.storeValueOf(component_stack_[stack_comp_id]->id(), value);
    }
  }

  Component & superComponentOf(StackLevel &lev) {
    assert(component_stack_.size() > lev.super_component());
    return *component_stack_[lev.super_component()];
  }

  unsigned component_stack_size() {
    return component_stack_.size();
  }

  void cleanRemainingComponentsOf(StackLevel &top) {
    while (component_stack_.size() > top.remaining_components_ofs()) {
      if (cache_.hasEntry(component_stack_.back()->id()))
        cache_.entry(component_stack_.back()->id()).set_deletable();
      delete component_stack_.back();
      component_stack_.pop_back();
    }
    assert(top.remaining_components_ofs() <= component_stack_.size());
  }

  Component & currentRemainingComponentOf(StackLevel &top) {
    assert(component_stack_.size() > top.currentRemainingComponent());
    return *component_stack_[top.currentRemainingComponent()];
  }
  
  // checks for the next yet to explore remaining component of top
  // returns true if a non-trivial non-cached component
  // has been found and is now stack_.TOS_NextComp()
  // returns false if all components have been processed;
  inline bool findNextRemainingComponentOf(StackLevel &top);

  inline void include_solution(int lit, StackLevel &top);

  inline void multiply_partial_solution(StackLevel &top);

  inline void include_partial_solution(int lit);

  inline void save_partial_solution();

  inline void remove_partial_solution(int lit);

  inline void recordRemainingCompsFor(StackLevel &top);

  inline void sortComponentStackRange(unsigned start, unsigned end);

  inline mpf_class get_saved_partial_sol();

  void gatherStatistics(){
//     statistics_.cache_bytes_memory_usage_ =
//	     cache_.recompute_bytes_memory_usage();
    cache_.compute_byte_size_infrasture();
  }

  void removeAllCachePollutionsOf(StackLevel &top);

  void * keyforHash;

  void getrandomkeyforhash(){
    keyforHash =
    get_random_key_for_clhash(UINT64_C(0x23a23cf5033c3c81),UINT64_C(0xb3816f6a2c68e530));
}

private:
  SolverConfiguration &config_;
  DataAndStatistics &statistics_;

  vector<Component *> component_stack_;
  ComponentCache cache_;
  ComponentAnalyzer ana_;
  vector<Variable> & variables_;
  vector<float> cachescore_;

  mpf_class partial_solution_;
  mpf_class saved_partial_solution_;
};


void ComponentManager::sortComponentStackRange(unsigned start, unsigned end){
    assert(start <= end);
    // sort the remaining components for processing
    for (unsigned i = start; i < end; i++)
      for (unsigned j = i + 1; j < end; j++) {
        if (component_stack_[i]->num_variables()
            < component_stack_[j]->num_variables())
          swap(component_stack_[i], component_stack_[j]);
      }
  }
bool ComponentManager::findNextRemainingComponentOf(StackLevel &top) {
    
    if (component_stack_.size() <= top.remaining_components_ofs())
      recordRemainingCompsFor(top);

    assert(!top.branch_found_unsat());
    // cout << "Number of components "<< top.remaining_components_ofs() << endl;
    if (top.hasUnprocessedComponents()) {
#ifdef VERB
        cout << "findNextRemainingComponentOf: Has unprocessed component" << endl;
#endif
      return true;
    }
    // if no component remains
    // make sure, at least that the current branch is considered SAT
    // cout << "Variable component management: " << top.getbranchvarsigned()<< " "<< partial_solution_ <<endl;
    // include_partial_solution(top.getbranchvarsigned());
    top.includeSolution(partial_solution_);
    partial_solution_ = 1;
    #ifdef VERB
    cout << "findNextRemainingComponentOf "<<partial_solution_<< endl;
    cout << "findNextRemainingComponentOf: no component remains" << endl;
    #endif
    return false;
  }

void ComponentManager::save_partial_solution(){
  saved_partial_solution_ = partial_solution_;
  cout << "Save Partial Solution: "<< partial_solution_<< endl;
  partial_solution_ = 1;
  // cout << "Partial Solution: "<< partial_solution_<< endl;
}

mpf_class ComponentManager::get_saved_partial_sol(){
  return saved_partial_solution_;
}

void ComponentManager::multiply_partial_solution(StackLevel &top){
  // cout << "Partial Solution "<< partial_solution_ << endl; 
  top.includeSolution(partial_solution_);
  partial_solution_ = 1;
}
void ComponentManager::remove_partial_solution(int lit){
  mpf_class w = 0;
  if(lit > 0){
    w = variables_[lit].get_weight();
    if(!mpf_cmp_d(w.get_mpf_t(), 2)){
      partial_solution_ /= 1;
    }
    else{
      partial_solution_ /= w;
    }
  }
  else{
    w = variables_[-1*lit].get_weight();
    if(!mpf_cmp_d(w.get_mpf_t(), 2)){
      partial_solution_ /= 1;
    }
    else{
      partial_solution_ /= 1 - w;
    }
  }
}
void ComponentManager::include_solution(int lit, StackLevel &top){
  mpf_class w = 0;
  if(lit > 0){
    w = variables_[lit].get_weight();
    if(!mpf_cmp_d(w.get_mpf_t(), 2)){
      top.includeSolution(1);
    }
    else{
      top.includeSolution(w);
    }
  }
  else{
    w = variables_[-1*lit].get_weight();
    if(!mpf_cmp_d(w.get_mpf_t(), 2)){
      top.includeSolution(1);
    }
    else{
      top.includeSolution(1 - w);
    }
  }
}

void ComponentManager::include_partial_solution(int lit){
  // cout << "include partial solution "<< lit <<"partial solution before "<<partial_solution_<<
  // " partial solution after ";
  mpf_class w = 0;
  if(lit > 0){
    w = variables_[lit].get_weight();
    if(!mpf_cmp_d(w.get_mpf_t(), 2)){
      partial_solution_ *= 1;
    }
    else{
      partial_solution_ *= w;
    }
  }
  else{
    w = variables_[-1*lit].get_weight();
    if(!mpf_cmp_d(w.get_mpf_t(), 2)){
      partial_solution_ *= 1;
    }
    else{
      partial_solution_ *= 1 - w;
    }
  }
  // cout << partial_solution_ << endl;
  // if (lit == 1 || lit == -1){
  //   cout <<"Weight of 1 "<< w << endl; 
  // }
}

void ComponentManager::recordRemainingCompsFor(StackLevel &top) {
    // cout <<"[ComponentManager::recordRemainingCompsFor]"<<endl;
   Component & super_comp = superComponentOf(top);
   unsigned new_comps_start_ofs = component_stack_.size();

   ana_.setupAnalysisContext(top, super_comp);
   for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++)
     if (ana_.isUnseenAndActive(*vt) &&
         ana_.exploreRemainingCompOf(*vt)){
       Component *p_new_comp = ana_.makeComponentFromArcheType();
       CacheableComponent *packed_comp = NULL;
       if (config_.useIsomorphicComponentCaching){
        //  cout << "New Component"<< endl;
        //  ana_.getArchetype().current_comp_for_caching_.printHash();
         if (ana_.getArchetype().current_comp_for_caching_.getuseisomorphism() ==false){
          packed_comp = new CacheableComponent(keyforHash,
         ana_.getArchetype().current_comp_for_caching_.getHash());
         }
         else{
          //  cout << "Num variables "<< ana_.getArchetype().current_comp_for_caching_.num_variables()<< endl;
           packed_comp = new CacheableComponent(ana_.getArchetype().current_comp_for_caching_.getblissHash(),
           ana_.getArchetype().current_comp_for_caching_.num_variables());
         }
       }
       
       else if (config_.usecachetencoding){
         packed_comp = new CacheableComponent(keyforHash,
         ana_.getArchetype().current_comp_for_caching_.getHash());
        //  cout << "New Component"<< endl;
        //  ana_.getArchetype().current_comp_for_caching_.printHash();
       }
      // //  else 
      else if (config_.usepcc){
        packed_comp = new CacheableComponent(keyforHash,
        ana_.getArchetype().current_comp_for_caching_);
       }
       else{
         packed_comp = new CacheableComponent(ana_.getArchetype().current_comp_for_caching_);
       }
         if (!cache_.manageNewComponent(top, *packed_comp)){
            component_stack_.push_back(p_new_comp);
            p_new_comp->set_id(cache_.storeAsEntry(*packed_comp, super_comp.id()));
         }
         else {
           if(config_.allowactivitydecrease){
            statistics_.numcachedec_++;
				    if(statistics_.numcachedec_ % 128 == 0){
					    increasecachescores();
			    	}
            for (vector<VariableIndex>::const_iterator it = p_new_comp->varsBegin(); 
            *it != varsSENTINEL; it++){
              cachescore_[*it] -= 1;
            }
           }
           delete packed_comp;
           delete p_new_comp;
         }
     }

   top.set_unprocessed_components_end(component_stack_.size());
   sortComponentStackRange(new_comps_start_ofs, component_stack_.size());
}

#endif /* COMPONENT_MANAGEMENT_H_ */
