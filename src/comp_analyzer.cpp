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
#include "mpreal.h"

using std::make_pair;

template class CompAnalyzer<mpz_class>;
template class CompAnalyzer<mpfr::mpreal>;

// Builds occ lists and sets things up, Done exactly ONCE for a whole counting runkk
// this sets up unif_occ and unif_occ_offs
template<typename T>
void CompAnalyzer<T>::initialize(
    const LiteralIndexedVector<LitWatchList> & watches, // binary clauses
    const ClauseAllocator<T>* alloc, const vector<ClauseOfs>& long_irred_cls) // longer-than-2-long clauses
{
  max_var = watches.end_lit().var() - 1;
  comp_vars.reserve(max_var + 1);
  VAR_FREQ_DO(var_freq_scores.resize(max_var + 1, 0));

  // maps var -> [cl_id, var1, var2, cl_id, var1, var2 ...]
  vector<vector<uint32_t>>  occ_ternary_clauses(max_var + 1);

  // maps var -> [var1..varn, SENTINEL_LIT, var1...varn, SENTINEL_LIT, ...]
  vector<vector<uint32_t>>  occ_long_clauses(max_var + 1);

  // maps var -> [cl_id, offset in occ_long_clauses, cl_id, offset in ...]
  vector<vector<ClauseOfs>> occs(max_var + 1);

  debug_print(COLBLBACK "Building occ list in CompAnalyzer<T>::initialize...");

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
  debug_print(COLBLBACK "Built occ list in CompAnalyzer<T>::initialize.");

  archetype.init_data(max_var, max_clid);

  debug_print(COLBLBACK "Building unified link list in CompAnalyzer<T>::initialize...");
  // the unified link list
  // This is an array that contains, flattened:
  // [  [vars of binary clauses],
  //    [cl_ids and lits] of tri clauses]
  //    [cl_id, offset in occs+offset in unif_occ]
  //    [the occ_long_clauses] ]
  unif_occ.clear();

  // a map into unif_occ.
  // maps var -> starting point in unif_occ
  unif_occ_offs.clear();
  unif_occ_offs.resize(max_var + 1, 0);

  for (uint32_t v = 1; v < max_var + 1; v++) {
    vector<uint32_t> lits_here(2*(max_var+1), 0);
    unif_occ_offs[v] = unif_occ.size();

    // data for binary clauses
    for(uint32_t i = 0; i < 2; i++) {
      for (const auto& bincl: watches[Lit(v, i)].binaries) {
        if (bincl.irred()) {
          unif_occ.push_back(bincl.lit().var());
        }
      }
    }

    // data for ternary clauses
    unif_occ.push_back(0);
    for(uint32_t i = 0; i < occ_ternary_clauses[v].size();) {
      auto cl_id = occ_ternary_clauses[v][i++];
      unif_occ.push_back(cl_id);
      Lit l;
      l = Lit::toLit(occ_ternary_clauses[v][i++]);
      unif_occ.push_back(l.raw());
      lits_here[l.raw()]++;
      l = Lit::toLit(occ_ternary_clauses[v][i++]);
      unif_occ.push_back(l.raw());
      lits_here[l.raw()]++;
    }

    // data for long clauses
    unif_occ.push_back(0);
    for(auto it = occs[v].begin(); it != occs[v].end(); it+=2) { // +2 because [cl_id, offset]
      auto cl_id = *it;
      unif_occ.push_back(cl_id);
      auto offs = *(it + 1) + (occs[v].end() - it);
      /* cout << "clid" << cl_id << " offs " << offs << endl; */
      unif_occ.push_back(offs);
    }

    unif_occ.push_back(0);
    for(const auto& raw: occ_long_clauses[v]) {
      Lit l = Lit::toLit(raw);
      unif_occ.push_back(l.raw());
      if (l != SENTINEL_LIT) lits_here[l.raw()]+=2;
    }
  }

  debug_print(COLBLBACK "Built unified link list in CompAnalyzer<T>::initialize.");
}

// returns true, iff the comp found is non-trivial
template<typename T>
bool CompAnalyzer<T>::explore_comp(const uint32_t v) {
  SLOW_DEBUG_DO(assert(archetype.var_unvisited_in_sup_comp(v)));
  record_comp(v); // sets up the component that "v" is in

  if (comp_vars.size() == 1) {
    debug_print("in " <<  __FUNCTION__ << " with single var: " <<  v);
    if (v >= indep_support_end) archetype.stack_level().include_solution(1);
    else {
      if (weighted()) archetype.stack_level().include_solution(counter->get_weight(v));
      else archetype.stack_level().include_solution(2);
    }
    archetype.set_var_in_peer_comp(v);
    return false;
  }
  return true;
}

// Create a component based on variable provided
template<typename T>
void CompAnalyzer<T>::record_comp(const uint32_t var) {
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
    uint32_t const* p = begin_cls_of_var(v);
    for (; *p; p++) {
      // NOTE: This below gives 10% slowdown(!) just to count the number of binary cls
      /* BUDDY_DO(if (counter->val(*p) == X_TRI) archetype.num_bin_cls++); */
      if (manage_occ_of(*p)) {
        if (!conf.do_check_unkn_bin || is_unknown(*p)) {
          VAR_FREQ_DO(bump_freq_score(*p); bump_freq_score(v));
        }
      }
    }

    //traverse ternary clauses
    for (p++; *p ; p+=3) {
      const auto clid = *p;
      const Lit a = *(Lit*)(p + 1);
      const Lit b = *(Lit*)(p + 2);
      if (archetype.clause_unvisited_in_sup_comp(*p)){
        /* cout << "Tern cl. (-?" << v << ") " << litA << " " << litB << endl; */
        if(is_true(a)|| is_true(b)) {
          archetype.clear_cl(clid);
        } else {
          VAR_FREQ_DO(bump_freq_score(v));
          manage_occ_and_score_of(a.var());
          manage_occ_and_score_of(b.var());
          archetype.set_clause_visited(clid ,is_unknown(a) && is_unknown(b));
        }
      }
    }

    // traverse long clauses
    for (p++; *p ; p +=2)
      if (archetype.clause_unvisited_in_sup_comp(*p)) {
        search_clause(v, *p, (Lit const*)(p + 1 + *(p+1)));
      }

  }

  debug_print(COLWHT "-> Went through all bin/tri/long and now comp_vars is "
      << comp_vars.size() << " long");
}
