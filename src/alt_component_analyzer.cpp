/*
 * alt_component_analyzer.cpp
 *
 *  Created on: Mar 5, 2013
 *      Author: mthurley
 */

#include "alt_component_analyzer.h"

// Builds occurrence lists and sets things up
void ComponentAnalyzer::initialize(
    LiteralIndexedVector<LitWatchList> & litWatchList, // binary clauses
    vector<Lit> &lit_pool) // longer-than-2-long clauses
{
  max_variable_id_ = litWatchList.end_lit().var() - 1;
  search_stack_.reserve(max_variable_id_ + 1);
  var_frequency_scores_.resize(max_variable_id_ + 1, 0);


  // maps var -> [cl_id, var1, var2, cl_id, var1, var2 ...]
  vector<vector<unsigned>>  occ_ternary_clauses(max_variable_id_ + 1);

  // maps var -> [var1..varn, SENTINEL_LIT, var1...varn, SENTINEL_LIT, ...]
  vector<vector<unsigned>>  occ_long_clauses(max_variable_id_ + 1);

  // maps var -> [cl_id, offset in occ_long_clauses, cl_id, offset in ...]
  vector<vector<ClauseOfs>> occs(max_variable_id_ + 1);

  print_debug(COLBLBACK "Building occurrence list in ComponentAnalyzer::initialize");

  vector<unsigned> tmp;
  max_clause_id_ = 0;
  auto it_curr_cl_st = lit_pool.begin();

  // lit_pool contains all non-binary clauses
  for (auto it_lit = lit_pool.begin(); it_lit != lit_pool.end(); it_lit++) {
    // Builds the occurrence list for 3-long and long clauses
    // it_curr_cl_st is the starting point of the clause
    // for each lit in the clause, it adds the clause to the occurrence list

    if (*it_lit == SENTINEL_LIT) { //End of this clause
      if (it_lit + 1 == lit_pool.end()) break;
      max_clause_id_++;
      it_lit += ClauseHeader::overheadInLits();
      it_curr_cl_st = it_lit + 1; // Point to next clause
    } else {
      const unsigned var = it_lit->var();
      assert(var <= max_variable_id_);
      getClause(tmp, it_curr_cl_st, *it_lit);
      assert(tmp.size() > 1);

      if(tmp.size() == 2) {
        // Ternary clause (but "tmp" is missing *it_lit, so it' of size 2)
        occ_ternary_clauses[var].push_back(max_clause_id_);
        occ_ternary_clauses[var].insert(
            occ_ternary_clauses[var].end(),
            tmp.begin(), tmp.end());
      } else {
        // Long clauses
        occs[var].push_back(max_clause_id_);
        occs[var].push_back(occ_long_clauses[var].size());
        occ_long_clauses[var].insert(occ_long_clauses[var].end(),
            tmp.begin(), tmp.end());
        occ_long_clauses[var].push_back(SENTINEL_LIT.raw());
      }
    }
  }

  archetype_.initSeen(max_variable_id_, max_clause_id_);

  // the unified link list
  // This is an array that contains, flattened:
  // [  [vars of binary clauses],
  //    [cl_ids and lits] of tri clauses]
  //    [cl_id, offset in occs+offset in unified_variable_links_lists_pool]
  //    [the occ_long_clauses] ]
  unified_variable_links_lists_pool_.clear();

  // a map into unified_variable_Links_lists_pool.
  // maps var -> starting point in unified_variable_Links_lists_pool
  variable_link_list_offsets_.clear();
  variable_link_list_offsets_.resize(max_variable_id_ + 1, 0);

  // now fill unified link list, for each variable
  for (unsigned v = 1; v < max_variable_id_ + 1; v++) {
    variable_link_list_offsets_[v] = unified_variable_links_lists_pool_.size();

    // data for binary clauses
    for (const auto& l: litWatchList[Lit(v, false)].binary_links_)
      if (l != SENTINEL_LIT)
        unified_variable_links_lists_pool_.push_back(l.var());

    for (const auto& l: litWatchList[Lit(v, true)].binary_links_)
      if (l != SENTINEL_LIT)
        unified_variable_links_lists_pool_.push_back(l.var());

    unified_variable_links_lists_pool_.push_back(0);

    // data for ternary clauses
    unified_variable_links_lists_pool_.insert(
        unified_variable_links_lists_pool_.end(),
        occ_ternary_clauses[v].begin(),
        occ_ternary_clauses[v].end());

    unified_variable_links_lists_pool_.push_back(0);

    // data for long clauses
    for(auto it = occs[v].begin(); it != occs[v].end(); it+=2) { // +2 because [cl_id, offset]
      unified_variable_links_lists_pool_.push_back(*it); //cl_id
      unified_variable_links_lists_pool_.push_back(*(it + 1) + (occs[v].end() - it));
    }

    unified_variable_links_lists_pool_.push_back(0);

    unified_variable_links_lists_pool_.insert(
        unified_variable_links_lists_pool_.end(),
        occ_long_clauses[v].begin(),
        occ_long_clauses[v].end());
  }
}

// Check which component a variable is in
void ComponentAnalyzer::recordComponentOf(const VariableIndex var) {
  search_stack_.clear();
  setSeenAndStoreInSearchStack(var);

  for (auto vt = search_stack_.begin(); vt != search_stack_.end(); vt++) {
    assert(isUnknown(*vt));

    //traverse binary clauses
    unsigned const* p = beginOfLinkList(*vt);
    for (; *p; p++) {
      if(manageSearchOccurrenceOf(Lit(*p,true))){
        var_frequency_scores_[*p]++;
        var_frequency_scores_[*vt]++;
      }
    }

    //traverse ternary clauses
    for (p++; *p ; p+=3) {
      if(archetype_.clause_unseen_in_sup_comp(*p)){
        const Lit litA = *(Lit*)(p + 1);
        const Lit litB = *(Lit*)(p + 2);
        if(isTrue(litA)|| isTrue(litB))
          archetype_.setClause_nil(*p);
        else {
          var_frequency_scores_[*vt]++;
          manageSearchOccurrenceAndScoreOf(litA);
          manageSearchOccurrenceAndScoreOf(litB);
          archetype_.setClause_seen(
              *p
              ,isUnknown(litA) && isUnknown(litB));
        }
      }
    }

    // traverse long clauses
    for (p++; *p ; p +=2)
      if(archetype_.clause_unseen_in_sup_comp(*p))
        searchClause(*vt,*p, (Lit const*)(p + 1 + *(p+1)));
  }
}
