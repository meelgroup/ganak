/*
 * new_component_analyzer.h
 *
 *  Created on: Mar 1, 2013
 *      Author: mthurley
 */

#ifndef NEW_COMPONENT_ANALYZER_H_
#define NEW_COMPONENT_ANALYZER_H_





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


struct CAClauseHeader {
  unsigned clause_id = 0;
  LiteralID lit_A;
  LiteralID lit_B;

  static unsigned overheadInLits() {
    return (sizeof(CAClauseHeader) - 2 * sizeof(LiteralID)) / sizeof(LiteralID);
  }
};


class NewComponentAnalyzer;


class NewComponentAnalyzer {
public:
	NewComponentAnalyzer(DataAndStatistics &statistics,
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

  //begin DEBUG
  void test_checkArchetypeRepForClause(unsigned *pcl_ofs){
      ClauseIndex clID = getClauseID(*pcl_ofs);
      bool all_a = true;
      for (auto itL = beginOfClause(*pcl_ofs); *itL != SENTINEL_LIT; itL++) {
        if(!isActive(*itL))
          all_a = false;
      }
      assert(all_a == archetype_.clause_all_lits_active(clID));
  }
  //end DEBUG

private:
  DataAndStatistics &statistics_;

  // the id of the last clause
  // note that clause ID is the clause number,
  // different from the offset of the clause in the literal pool
  unsigned max_clause_id_ = 0;
  unsigned max_variable_id_ = 0;

  // pool of clauses as lists of LiteralIDs
  // Note that with a clause begin position p we have
  // the clause ID at position p-1
  vector<LiteralID> literal_pool_;

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

  vector<unsigned> map_clause_id_to_ofs_;
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

  bool isSatisfiedByFirstTwoLits(ClauseOfs cl_ofs) {
      return isSatisfied(getHeaderOf(cl_ofs).lit_A) || isSatisfied(getHeaderOf(cl_ofs).lit_B);
    }

  bool isActive(const VariableIndex v) {
    return literal_values_[LiteralID(v, true)] == X_TRI;
  }

  unsigned getClauseID(ClauseOfs cl_ofs) {
    return reinterpret_cast<CAClauseHeader *>(&literal_pool_[cl_ofs
        - CAClauseHeader::overheadInLits()])->clause_id;
  }

  CAClauseHeader &getHeaderOf(ClauseOfs cl_ofs) {
    return *reinterpret_cast<CAClauseHeader *>(&literal_pool_[cl_ofs
        - CAClauseHeader::overheadInLits()]);
  }

  unsigned *beginOfLinkList(VariableIndex v) {
    return &unified_variable_links_lists_pool_[variable_link_list_offsets_[v]];
  }

  vector<LiteralID>::iterator beginOfClause(ClauseOfs cl_ofs) {
    return literal_pool_.begin() + cl_ofs;
  }

  // stores all information about the component of var
  // in variables_seen_, clauses_seen_ and
  // component_search_stack
  // we have an isolated variable iff
  // after execution component_search_stack.size()==1
  void recordComponentOf(const VariableIndex var);

  void pushLitsInto(vector<unsigned> &target,
		       const vector<LiteralID> &lit_pool,
		       unsigned start_of_cl,
		       LiteralID & omitLit){
	  for (auto it_lit = lit_pool.begin() + start_of_cl;
			  (*it_lit != SENTINEL_LIT); it_lit++) {
		  if(it_lit->var() != omitLit.var())
			  target.push_back(it_lit->raw());
	  }
	  target.push_back(SENTINEL_LIT.raw());
  }

};



#endif /* NEW_COMPONENT_ANALYZER_H_ */
