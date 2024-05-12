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
#include "cryptominisat5/solvertypesmini.h"
#include "structures.hpp"

using std::make_pair;

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

  archetype.init_data(max_var, max_clid);

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
  map<uint32_t, Lit> best_alters;
  for (uint32_t v = 1; v < max_var + 1; v++) {
    vector<uint32_t> lits_here(2*(max_var+1), 0);
    variable_link_list_offsets[v] = unified_var_links_lists_pool.size();

    // data for binary clauses
    for(uint32_t i = 0; i < 2; i++) {
      for (const auto& bincl: watches[Lit(v, i)].binaries) {
        if (bincl.irred()) {
          unified_var_links_lists_pool.push_back(bincl.lit().var());
          lits_here[bincl.lit().raw()]++;
        }
      }
    }

    // data for ternary clauses
    unified_var_links_lists_pool.push_back(0);
    for(uint32_t i = 0; i < occ_ternary_clauses[v].size();) {
      auto cl_id = occ_ternary_clauses[v][i++];
      unified_var_links_lists_pool.push_back(cl_id);
      Lit l;
      l = Lit::toLit(occ_ternary_clauses[v][i++]);
      unified_var_links_lists_pool.push_back(l.raw());
      lits_here[l.raw()]++;
      l = Lit::toLit(occ_ternary_clauses[v][i++]);
      unified_var_links_lists_pool.push_back(l.raw());
      lits_here[l.raw()]++;
    }

    // data for long clauses
    unified_var_links_lists_pool.push_back(0);
    for(auto it = occs[v].begin(); it != occs[v].end(); it+=2) { // +2 because [cl_id, offset]
      auto cl_id = *it;
      unified_var_links_lists_pool.push_back(cl_id);
      auto offs = *(it + 1) + (occs[v].end() - it);
      /* cout << "clid" << cl_id << " offs " << offs << endl; */
      unified_var_links_lists_pool.push_back(offs);
    }

    unified_var_links_lists_pool.push_back(0);
    for(const auto& raw: occ_long_clauses[v]) {
      Lit l = Lit::toLit(raw);
      unified_var_links_lists_pool.push_back(l.raw());
      if (l != SENTINEL_LIT) lits_here[l.raw()]++;
    }

    Lit best = Lit(1, 0);;
    uint32_t best_num = 0;
    for(uint32_t i = 2; i < 2*(max_var+1); i++) {
      if (best == SENTINEL_LIT || best_num <= lits_here[i]) {
        best = Lit::toLit(i);
        best_num = lits_here[i];
      }
    }
    best_alters[v] = best;
  }

  /* cout << "ALT!!!!!!!!!!!!!!!" << endl; */
  variable_link_list_offsets_alt.clear();
  variable_link_list_offsets_alt.resize(max_var +1);
  solver->vivif_setup();
  for (uint32_t v = 1; v < max_var + 1; v++) {
    const Lit true_l = best_alters[v];
    solver->v_new_lev();
    solver->v_enqueue(true_l);
    bool ret = solver->v_propagate();
    assert(ret);
    variable_link_list_offsets_alt[v] = make_pair(true_l, unified_var_links_lists_pool.size());

    // data for binary clauses
    for(uint32_t i = 0; i < 2; i++) {
      for (const auto& bincl: watches[Lit(v, i)].binaries) {
        if (bincl.irred() && bincl.lit() != true_l) {
          unified_var_links_lists_pool.push_back(bincl.lit().var());
        }
      }
    }

    // data for ternary clauses
    unified_var_links_lists_pool.push_back(0);
    for(uint32_t i = 0; i < occ_ternary_clauses[v].size();) {
      if ((Lit::toLit(occ_ternary_clauses[v][i+1]) == true_l) ||
          Lit::toLit(occ_ternary_clauses[v][i+2]) == true_l) {
        i += 3;
        continue;
      }
      auto cl_id = occ_ternary_clauses[v][i++];
      unified_var_links_lists_pool.push_back(cl_id);
      Lit l;
      l = Lit::toLit(occ_ternary_clauses[v][i++]);
      unified_var_links_lists_pool.push_back(l.raw());
      l = Lit::toLit(occ_ternary_clauses[v][i++]);
      unified_var_links_lists_pool.push_back(l.raw());
    }

    // data for long clauses
    vector<pair<uint32_t, uint32_t>> cl_ids;
    unified_var_links_lists_pool.push_back(0);
    for(auto it = occs[v].begin(); it != occs[v].end(); it+=2) { // +2 because [cl_id, offset]
      auto cl_id = *it;
      const Clause& cl = *alloc->ptr(long_irred_cls[cl_id-1]);
      bool sat = false;
      for(const auto& l: cl) if (solver->v_val(l) == T_TRI) sat = true;
      if (sat) continue;
      cl_ids.push_back(make_pair(cl_id, unified_var_links_lists_pool.size()+1));
      unified_var_links_lists_pool.push_back(cl_id);
      unified_var_links_lists_pool.push_back(SENTINEL_LIT.raw()); //placeholder
    }

    unified_var_links_lists_pool.push_back(0);
    for(const auto& cl_id: cl_ids) {
      unified_var_links_lists_pool[cl_id.second] = unified_var_links_lists_pool.size()-cl_id.second;
      const Clause& cl = *alloc->ptr(long_irred_cls[cl_id.first-1]);
      for(const auto& l: cl) {
        if (l.var() == v) continue;
        if (solver->v_val(l) == F_TRI) continue;
        unified_var_links_lists_pool.push_back(l.raw());
      }
      unified_var_links_lists_pool.push_back(SENTINEL_LIT.raw());
    }
    solver->v_backtrack();
  }

  debug_print(COLBLBACK "Built unified link list in CompAnalyzer::initialize.");
}

// returns true, iff the comp found is non-trivial
bool CompAnalyzer::explore_comp(const uint32_t v) {
  SLOW_DEBUG_DO(assert(archetype.var_unvisited_in_sup_comp(v)));
  record_comp(v); // sets up the component that "v" is in

  if (comp_vars.size() == 1) {
    debug_print("in " <<  __FUNCTION__ << " with single var: " <<  v);
    if (v >= indep_support_end) archetype.stack_level().includeSolution(1);
    else archetype.stack_level().includeSolution(2);
    archetype.set_var_in_peer_comp(v);
    return false;
  }
  return true;
}

// Create a component based on variable provided
void CompAnalyzer::record_comp(const uint32_t var) {
  SLOW_DEBUG_DO(assert(is_unknown(var)));
  comp_vars.clear();
  comp_vars.push_back(var);
  archetype.set_var_visited(var);

  debug_print(COLWHT "We are NOW going through all binary/tri/long clauses "
      "recursively and put into search_stack_ all the variables that are connected to var: " << var);

  // manage_occ_of and search_clause
  // will push into search_stack_ which will make this
  // a recursive search for all clauses & variables that this variable is connected to
  for (auto vt = comp_vars.begin(); vt != comp_vars.end(); vt++) {
    const auto v = *vt;
    SLOW_DEBUG_DO(assert(is_unknown(v)));

    //traverse binary clauses
    uint32_t const* p;
    if (is_true(variable_link_list_offsets_alt[v].first)) p = begin_cls_of_var_alt(v);
    else p = begin_cls_of_var(v);
    for (; *p; p++) {
      // NOTE: This below gives 10% slowdown(!) just to count the number of binary cls
      /* BUDDY_DO(if (solver->val(*p) == X_TRI) archetype.num_bin_cls++); */
      manage_occ_of(*p);
    }

    //traverse ternary clauses
    for (p++; *p ; p+=3) {
      auto clid = *p;
      const Lit a = *(Lit*)(p + 1);
      const Lit b = *(Lit*)(p + 2);
      if (archetype.clause_unvisited_in_sup_comp(*p)){
        /* cout << "Tern cl. (-?" << v << ") " << litA << " " << litB << endl; */
        if(is_true(a)|| is_true(b)) {
          archetype.clear_cl(clid);
        } else {
          manage_occ_of(a.var());
          manage_occ_of(b.var());
          archetype.set_clause_visited(clid ,is_unknown(a) && is_unknown(b));
        }
      }
    }

    // traverse long clauses
    for (p++; *p ; p +=2)
      if (archetype.clause_unvisited_in_sup_comp(*p)) {
        search_clause(*p, (Lit const*)(p + 1 + *(p+1)));
      }
  }

  debug_print(COLWHT "-> Went through all bin/tri/long and now comp_vars is "
      << comp_vars.size() << " long");
}
