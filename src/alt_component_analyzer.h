/*
 * alt_component_analyzer.h
 *
 *  Created on: Mar 5, 2013
 *      Author: mthurley
 */

#ifndef ALT_COMPONENT_ANALYZER_H_
#define ALT_COMPONENT_ANALYZER_H_



#include "statistics.h"
#include "component_types/component.h"
#include "component_types/base_packed_component.h"
#include "component_types/component_archetype.h"



#include <vector>
#include <cmath>
#include <gmpxx.h>
#include "containers.h"
#include "stack.h"

using namespace std;

class AltComponentAnalyzer {
public:
	AltComponentAnalyzer(DataAndStatistics &statistics,
        LiteralIndexedVector<TriValue> & lit_values) :
        statistics_(statistics), literal_values_(lit_values) {
  }

  unsigned scoreOf(VariableIndex v) {
    return var_frequency_scores_[v];
  }

  ComponentArchetype &current_archetype(){
    return archetype_;
  }

  void initialize(LiteralIndexedVector<Literal> & literals,
      vector<LiteralID> &lit_pool);


  bool isUnseenAndActive(VariableIndex v){
    assert(v <= max_variable_id_);
    return archetype_.var_unseen_in_sup_comp(v);
  }

  // manages the literal whenever it occurs in component analysis
  // returns true iff the underlying variable was unseen before
  //
  bool manageSearchOccurrenceOf(LiteralID lit){
      if(archetype_.var_unseen_in_sup_comp(lit.var())){
        search_stack_.push_back(lit.var());
        archetype_.setVar_seen(lit.var());
        return true;
      }
      return false;
    }
  bool manageSearchOccurrenceAndScoreOf(LiteralID lit){
    var_frequency_scores_[lit.var()]+= isActive(lit);
    return manageSearchOccurrenceOf(lit);
  }

  void setSeenAndStoreInSearchStack(VariableIndex v){
    assert(isActive(v));
    search_stack_.push_back(v);
        archetype_.setVar_seen(v);
  }


  void setupAnalysisContext(StackLevel &top, Component & super_comp){
     archetype_.reInitialize(top,super_comp);

     for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++)
       if (isActive(*vt)) {
         archetype_.setVar_in_sup_comp_unseen(*vt);
         var_frequency_scores_[*vt] = 0;
       }

     for (auto itCl = super_comp.clsBegin(); *itCl != clsSENTINEL; itCl++)
         archetype_.setClause_in_sup_comp_unseen(*itCl);
  }

  // returns true, iff the component found is non-trivial
  bool exploreRemainingCompOf(VariableIndex v) {
    assert(archetype_.var_unseen_in_sup_comp(v));
    recordComponentOf(v);

    if (search_stack_.size() == 1) {
      archetype_.stack_level().includeSolution(2);
      archetype_.setVar_in_other_comp(v);
      return false;
    }
    return true;
  }


  inline Component *makeComponentFromArcheType(){
    return archetype_.makeComponentFromState(search_stack_.size());
  }

  unsigned max_clause_id(){
     return max_clause_id_;
  }
  unsigned max_variable_id(){
    return max_variable_id_;
  }

  ComponentArchetype &getArchetype(){
    return archetype_;
  }

private:
  DataAndStatistics &statistics_;

  // the id of the last clause
  // note that clause ID is the clause number,
  // different from the offset of the clause in the literal pool
  unsigned max_clause_id_ = 0;
  unsigned max_variable_id_ = 0;

  // this contains clause offsets of the clauses
  // where each variable occurs in;
  vector<ClauseOfs> variable_occurrence_lists_pool_;

  // this is a new idea,
  // for every variable we have a list
  // 0 binarylinks 0 occurrences 0
  // this should give better cache behaviour,
  // because all links of one variable (binary and nonbinray) are found
  // in one contiguous chunk of memory
  vector<unsigned> unified_variable_links_lists_pool_;


  vector<unsigned> variable_link_list_offsets_;

  LiteralIndexedVector<TriValue> & literal_values_;

  vector<unsigned> var_frequency_scores_;

  ComponentArchetype  archetype_;

  vector<VariableIndex> search_stack_;

  bool isResolved(const LiteralID lit) {
    return literal_values_[lit] == F_TRI;
  }

  bool isSatisfied(const LiteralID lit) {
    return literal_values_[lit] == T_TRI;
  }
  bool isActive(const LiteralID lit) {
      return literal_values_[lit] == X_TRI;
  }

  bool isActive(const VariableIndex v) {
    return literal_values_[LiteralID(v, true)] == X_TRI;
  }

  unsigned *beginOfLinkList(VariableIndex v) {
    return &unified_variable_links_lists_pool_[variable_link_list_offsets_[v]];
  }

  // stores all information about the component of var
  // in variables_seen_, clauses_seen_ and
  // component_search_stack
  // we have an isolated variable iff
  // after execution component_search_stack.size()==1
  void recordComponentOf(const VariableIndex var);


  void getClause(vector<unsigned> &tmp,
   		       vector<LiteralID>::iterator & it_start_of_cl,
   		       LiteralID & omitLit){
  	  tmp.clear();
  	  for (auto it_lit = it_start_of_cl;*it_lit != SENTINEL_LIT; it_lit++) {
   		  if(it_lit->var() != omitLit.var())
   			 tmp.push_back(it_lit->raw());
   	  }
     }


  void searchClause(VariableIndex vt, ClauseIndex clID, LiteralID * pstart_cls){
    auto itVEnd = search_stack_.end();
    bool all_lits_active = true;
    for (auto itL = pstart_cls; *itL != SENTINEL_LIT; itL++) {
      assert(itL->var() <= max_variable_id_);
      if(!archetype_.var_nil(itL->var()))
        manageSearchOccurrenceAndScoreOf(*itL);
      else {
        assert(!isActive(*itL));
        all_lits_active = false;
        if (isResolved(*itL))
          continue;
        //BEGIN accidentally entered a satisfied clause: undo the search process
        while (search_stack_.end() != itVEnd) {
          assert(search_stack_.back() <= max_variable_id_);
          archetype_.setVar_in_sup_comp_unseen(search_stack_.back());
          search_stack_.pop_back();
        }
        archetype_.setClause_nil(clID);
        while(*itL != SENTINEL_LIT)
          if(isActive(*(--itL)))
            var_frequency_scores_[itL->var()]--;
        //END accidentally entered a satisfied clause: undo the search process
        break;
      }
    }

    if (!archetype_.clause_nil(clID)){
      var_frequency_scores_[vt]++;
      archetype_.setClause_seen(clID,all_lits_active);
    }
  }


//  void searchThreeClause(VariableIndex vt, ClauseIndex clID, LiteralID * pstart_cls){
//      auto itVEnd = search_stack_.end();
//      bool all_lits_active = true;
//      // LiteralID * pstart_cls = reinterpret_cast<LiteralID *>(p + 1 + *(p+1));
//      for (auto itL = pstart_cls; itL != pstart_cls+2; itL++) {
//        assert(itL->var() <= max_variable_id_);
//        if(archetype_.var_nil(itL->var())){
//          assert(!isActive(*itL));
//          all_lits_active = false;
//          if (isResolved(*itL))
//            continue;
//          //BEGIN accidentally entered a satisfied clause: undo the search process
//          while (search_stack_.end() != itVEnd) {
//            assert(search_stack_.back() <= max_variable_id_);
//            archetype_.setVar_in_sup_comp_unseen(search_stack_.back());
//            search_stack_.pop_back();
//          }
//          archetype_.setClause_nil(clID);
//          while(itL != pstart_cls - 1)
//            if(isActive(*(--itL)))
//              var_frequency_scores_[itL->var()]--;
//          //END accidentally entered a satisfied clause: undo the search process
//          break;
//        } else {
//          assert(isActive(*itL));
//          var_frequency_scores_[itL->var()]++;
//          if(isUnseenAndActive(itL->var()))
//            setSeenAndStoreInSearchStack(itL->var());
//        }
//      }
//
//      if (!archetype_.clause_nil(clID)){
//        var_frequency_scores_[vt]++;
//        archetype_.setClause_seen(clID,all_lits_active);
//      }
//    }
};


#endif /* ALT_COMPONENT_ANALYZER_H_ */
