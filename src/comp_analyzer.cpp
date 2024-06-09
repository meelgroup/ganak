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
// this sets up unif_occ
template<typename T>
void CompAnalyzer<T>::initialize(
    const LiteralIndexedVector<LitWatchList> & watches, // binary clauses
    const ClauseAllocator<T>* alloc, const vector<ClauseOfs>& long_irred_cls) // longer-than-2-long clauses
{
  max_var = watches.end_lit().var() - 1;
  comp_vars.reserve(max_var + 1);
  VAR_FREQ_DO(var_freq_scores.resize(max_var + 1, 0));

  vector<vector<ClauseOfs>> occs(max_var + 1);

  debug_print(COLBLBACK "Building occ list in CompAnalyzer<T>::initialize...");

  max_clid = 1;
  unif_occ.clear();
  unif_occ.resize(max_var + 1);
  long_clauses_data.clear();
  vector<uint32_t> tmp;
  for (const auto& off: long_irred_cls) {
    const Clause& cl = *alloc->ptr(off);
    assert(cl.size() > 2);
    uint32_t long_cl_off = long_clauses_data.size();
    for(const auto&l: cl) long_clauses_data.push_back(l);
    long_clauses_data.push_back(SENTINEL_LIT);

    for(const auto& l: cl) {
      const uint32_t var = l.var();
      assert(var <= max_var);
      unif_occ[var].push_back(max_clid);
      unif_occ[var].push_back(long_cl_off);
    }
    max_clid++;
  }
  debug_print(COLBLBACK "Built occ list in CompAnalyzer<T>::initialize.");

  archetype.init_data(max_var, max_clid);

  debug_print(COLBLBACK "Building unified link list in CompAnalyzer<T>::initialize...");


  // data for binary clauses
  unif_occ_bin.clear();
  unif_occ_bin.resize(max_var+1);
  vector<uint32_t> tmp2;
  for (uint32_t v = 1; v < max_var + 1; v++) {
    tmp2.clear();
    for(uint32_t i = 0; i < 2; i++) {
      for (const auto& bincl: watches[Lit(v, i)].binaries) {
        if (bincl.irred()) tmp2.push_back(bincl.lit().var());
      }
    }
    unif_occ_bin[v] = tmp2;
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
    for(const auto& p: unif_occ_bin[v]) {
      if (manage_occ_of(p)) {
        if (!conf.do_check_unkn_bin || is_unknown(p)) {
          VAR_FREQ_DO(bump_freq_score(p); bump_freq_score(v));
        }
      }
    }

    // traverse long clauses
    for (uint32_t i = 0; i < unif_occ[v].size(); i+=2) {
      const uint32_t cl_id = unif_occ[v][i];
      const ClauseOfs cl_off = unif_occ[v][i+1];
      if (archetype.clause_unvisited_in_sup_comp(cl_id)) {
        search_clause(v, cl_id, long_clauses_data.data()+cl_off);
      }
    }
  }

  debug_print(COLWHT "-> Went through all bin/tri/long and now comp_vars is "
      << comp_vars.size() << " long");
}
