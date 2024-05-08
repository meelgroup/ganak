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

#include "comp_analyzer.hpp"
#include "common.hpp"
#include "counter.hpp"
#include "clauseallocator.hpp"
#include "structures.hpp"

// There is exactly ONE of these
CompAnalyzer::CompAnalyzer(
        const LiteralIndexedVector<TriValue> & lit_values,
        const uint32_t& _indep_support_end,
        Counter* _solver) :
        conf(_solver->get_conf()),
        values(lit_values),
        indep_support_end(_indep_support_end),
        solver(_solver)
{}


// Builds occ lists and sets things up, Done exactly ONCE for a whole counting runkk
// this sets up unified_var_links_lists_pool and variable_link_list_offsets_
void CompAnalyzer::initialize(
    const LiteralIndexedVector<LitWatchList> & watches, // binary clauses
    const ClauseAllocator* alloc, const vector<ClauseOfs>& long_irred_cls) // longer-than-2-long clauses
{
  max_var = watches.end_lit().var() - 1;
  comp_vars.reserve(max_var + 1);

  // maps var -> [cl_id, var1, var2, cl_id, var1, var2 ...]
  vector<vector<uint32_t>>  occ_ternary_clauses(max_var + 1);

  // maps var -> [var1..varn, SENTINEL_LIT, var1...varn, SENTINEL_LIT, ...]
  vector<vector<uint32_t>>  occ_long_clauses(max_var + 1);

  // maps var -> [cl_id, offset in occ_long_clauses, cl_id, offset in ...]
  vector<vector<ClauseOfs>> occs(max_var + 1);

  debug_print(COLBLBACK "Building occ list in CompAnalyzer::initialize...");

  vector<uint32_t> tmp;
  max_clid = 1;
  assert(idx_to_cl_map.empty());
  assert(idx_to_cl_data.empty());
  idx_to_cl_map.push_back(0);
  // lit_pool contains all non-binary clauses
  for (const auto& off: long_irred_cls) {
    // Builds the occ list for 3-long and long clauses
    // it_curr_cl_st is the starting point of the clause
    // for each lit in the clause, it adds the clause to the occ list
    const Clause& cl = *alloc->ptr(off);
    assert(cl.size() > 2);

    for(const auto& l: cl) {
      const uint32_t var = l.var();
      assert(var <= max_var);
      get_cl(tmp, cl, l);
      assert(tmp.size() > 1);

      if(tmp.size() == 2) {
        // Ternary clause (but "tmp" is missing *it_lit, so it' of size 2)
        occ_ternary_clauses[var].push_back(max_clid);
        occ_ternary_clauses[var].insert(occ_ternary_clauses[var].end(), tmp.begin(), tmp.end());
      } else {
        // Long clauses
        occs[var].push_back(max_clid);
        occs[var].push_back(occ_long_clauses[var].size());
        occ_long_clauses[var].insert(occ_long_clauses[var].end(), tmp.begin(), tmp.end());
        occ_long_clauses[var].push_back(SENTINEL_LIT.raw());
      }
    }
    // Fill idx
    const uint32_t at = idx_to_cl_data.size();
    for(const auto& l: cl) idx_to_cl_data.push_back(l);
    idx_to_cl_data.push_back(SENTINEL_LIT);
    assert(idx_to_cl_map.size() == max_clid);
    idx_to_cl_map.push_back(at);
    max_clid++;
  }
  debug_print(COLBLBACK "Built occ list in CompAnalyzer::initialize.");

  archetype.init_seen(max_var, max_clid);

  debug_print(COLBLBACK "Building unified link list in CompAnalyzer::initialize...");
  // the unified link list
  // This is an array that contains, flattened:
  // [  [vars of binary clauses],
  //    [cl_ids and lits] of tri clauses]
  //    [cl_id, offset in occs+offset in unified_var_links_lists_pool]
  //    [the occ_long_clauses] ]
  unified_var_links_lists_pool.clear();

  // a map into unified_var_Links_lists_pool.
  // maps var -> starting point in unified_var_links_lists_pool
  variable_link_list_offsets.clear();
  variable_link_list_offsets.resize(max_var + 1, 0);

  // now fill unified link list, for each variable
  for (uint32_t v = 1; v < max_var + 1; v++) {
    variable_link_list_offsets[v] = unified_var_links_lists_pool.size();

    // data for binary clauses
    for(uint32_t i = 0; i < 2; i++) {
      for (const auto& bincl: watches[Lit(v, i)].binaries) {
        if (bincl.irred()) unified_var_links_lists_pool.push_back(bincl.lit().var());
      }
    }

    // data for ternary clauses
    unified_var_links_lists_pool.push_back(0);
    unified_var_links_lists_pool.insert(
        unified_var_links_lists_pool.end(),
        occ_ternary_clauses[v].begin(),
        occ_ternary_clauses[v].end());

    // data for long clauses
    unified_var_links_lists_pool.push_back(0);
    for(auto it = occs[v].begin(); it != occs[v].end(); it+=2) { // +2 because [cl_id, offset]
      unified_var_links_lists_pool.push_back(*it); //cl_id
      unified_var_links_lists_pool.push_back(*(it + 1) + (occs[v].end() - it));
    }

    unified_var_links_lists_pool.push_back(0);
    unified_var_links_lists_pool.insert(
        unified_var_links_lists_pool.end(),
        occ_long_clauses[v].begin(),
        occ_long_clauses[v].end());
  }
  debug_print(COLBLBACK "Built unified link list in CompAnalyzer::initialize.");
}

// returns true, iff the comp found is non-trivial
bool CompAnalyzer::explore_comp(const uint32_t v) {
  SLOW_DEBUG_DO(assert(archetype.var_unseen_in_sup_comp(v)));
  record_comp(v); // sets up the component that "v" is in

  // comp only contains one variable
  if (comp_vars.size() == 1) {
    debug_print("explore remaining with single var, v is: " <<  v);
    if (v >= indep_support_end) {
      archetype.stack_level().includeSolution(1);
      /* CHECK_COUNT_DO(assert(solver->check_count(true, v) == 1)); */
    } else {
      archetype.stack_level().includeSolution(2);
      /* CHECK_COUNT_DO(assert(solver->check_count(true, v) == 2)); */
    }
    archetype.set_var_in_peer_comp(v);
    return false;
  }
  return true;
}

// Create a component based on variable provided
void CompAnalyzer::record_comp(const uint32_t var) {
  comp_vars.clear();
  setSeenAndStoreInSearchStack(var);

  debug_print(COLWHT "We are NOW going through all binary/tri/long clauses "
      "recursively and put into search_stack_ all the variables that are connected to var: " << var);

  // manageSearchOccurrenceOf and search_clause
  // will push into search_stack_ which will make this
  // a recursive search for all clauses & variables that this variable is connected to
  for (auto vt = comp_vars.begin(); vt != comp_vars.end(); vt++) {
    const auto v = *vt;
    SLOW_DEBUG_DO(assert(is_unknown(v)));

    //traverse binary clauses
    uint32_t const* p = begin_cls_of_var(v);
    for (; *p; p++) {
      // NOTE: This below gives 10% slowdown(!) just to count the number of binary cls
      /* BUDDY_DO(if (solver->val(*p) == X_TRI) archetype.num_bin_cls++); */
      manageSearchOccurrenceOf(*p);
    }

    //traverse ternary clauses
    for (p++; *p ; p+=3) {
      auto clid = *p;
      const Lit a = *(Lit*)(p + 1);
      const Lit b = *(Lit*)(p + 2);
      if (archetype.clause_unseen_in_sup_comp(*p)){
        /* cout << "Tern cl. (-?" << v << ") " << litA << " " << litB << endl; */
        if(is_true(a)|| is_true(b)) {
          archetype.clear_cl(clid);
        } else {
          manageSearchOccurrenceOf(a.var());
          manageSearchOccurrenceOf(b.var());
          archetype.set_clause_seen(clid ,is_unknown(a) && is_unknown(b));
        }
      }
    }

    // traverse long clauses
    for (p++; *p ; p +=2)
      if (archetype.clause_unseen_in_sup_comp(*p))
        search_clause(v,*p, (Lit const*)(p + 1 + *(p+1)));
  }

  debug_print(COLWHT "-> Went through all bin/tri/long and now comp_vars is "
      << comp_vars.size() << " long");
}
