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
#include <climits>
#include <cstdint>

using std::make_pair;

template class CompAnalyzer<mpz_class>;
template class CompAnalyzer<mpfr::mpreal>;

inline std::ostream& operator<<(std::ostream& os, const ClData& d)
{
  os << "[id: " << d.id << " off: " << d.off << "]";
  /* os << "id: " << d.id; */
  return os;
}


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
  vector<vec<ClData>> unif_occ;
  unif_occ.clear();
  unif_occ.resize(max_var + 1);
  long_clauses_data.clear();
  long_clauses_data.push_back(SENTINEL_LIT); // MUST start with a sentinel!
  vector<uint32_t> tmp;
  for (const auto& off: long_irred_cls) {
    const Clause& cl = *alloc->ptr(off);
    assert(cl.size() > 2);
    uint32_t long_cl_off = long_clauses_data.size();
    if (cl.size() > 3) {
      Lit blk_lit = cl[cl.size()/2];
      for(const auto&l: cl) long_clauses_data.push_back(l);
      long_clauses_data.push_back(SENTINEL_LIT);

      for(const auto& l: cl) {
        const uint32_t var = l.var();
        assert(var <= max_var);
        ClData d;
        d.tri = false;
        d.id = max_clid;
        d.off = long_cl_off;
        d.blk_lit = blk_lit;
        unif_occ[var].push_back(d);
      }
    } else {
      for(const auto& l: cl) {
        uint32_t at = 0;
        Lit lits[2];
        for(const auto&l2: cl) if (l.var() != l2.var()) lits[at++] = l2;
        assert(at == 2);
        ClData d;
        d.tri = true;
        d.id = max_clid;
        d.blk_lit = lits[0];
        d.off = lits[1].raw();
        unif_occ[l.var()].push_back(d);
      }
    }
    max_clid++;
  }

  long_sz_declevs.resize(1);
  long_sz_declevs[0].resize(max_var+1, MemData());
  for(uint32_t var = 1; var < max_var+1; var++) {
    long_sz_declevs[0][var] = MemData(unif_occ[var].size());
    /* std::sort(unif_occ[var].begin(), unif_occ[var].end()); */
  }

  if (true) {
    // fill holder
    uint32_t total_sz = 0;
    for(const auto& u: unif_occ) total_sz += u.size()*(sizeof(ClData)/sizeof(uint32_t)) + 2;
    uint32_t* data = new uint32_t[total_sz];
    holder.data  = data;

    uint32_t* data_start = holder.data + unif_occ.size()*2;
    for(uint32_t v = 0; v < unif_occ.size(); v++) {
      const auto& u = unif_occ[v];

      holder.data[v*2+1] = u.size();
      uint32_t offs = data_start - holder.data;
      holder.data[v*2] = offs;
      assert(offs < total_sz);
      memcpy(data_start, u.data, u.size()*sizeof(ClData));
      data_start += u.size()*(sizeof(ClData)/sizeof(uint32_t));
    }
    assert(data_start == data + total_sz);

    for(uint32_t v = 0; v < unif_occ.size(); v++) {
      assert(unif_occ[v].size() == holder.size(v));
      for(uint32_t i = 0; i < unif_occ[v].size(); i++) {
        assert(unif_occ[v][i] == holder.begin(v)[i]);
      }
    }
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
    unif_occ_bin[v].clear();
    unif_occ_bin[v].resize(tmp2.size());
    for(uint32_t i = 0; i < tmp2.size(); i++) unif_occ_bin[v][i] = tmp2[i];
  }
  last_seen.resize(max_var+1, 0);

  debug_print(COLBLBACK "Built unified link list in CompAnalyzer<T>::initialize.");
}

// returns true, iff the comp found is non-trivial
template<typename T>
bool CompAnalyzer<T>::explore_comp(const uint32_t v, int32_t dec_lev) {
  SLOW_DEBUG_DO(assert(archetype.var_unvisited_in_sup_comp(v)));
  record_comp(v, dec_lev); // sets up the component that "v" is in

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
void CompAnalyzer<T>::record_comp(const uint32_t var, int32_t declev) {
  SLOW_DEBUG_DO(assert(is_unknown(var)));
  comp_vars.clear();
  comp_vars.push_back(var);
  archetype.set_var_visited(var);

  debug_print(COLWHT "We are NOW going through all binary/tri/long clauses "
      "recursively and put into search_stack_ all the variables that are connected to var: " << var);

  if (declev >= (int)long_sz_declevs.size()) {
    long_sz_declevs.resize(declev+1);
    long_sz_declevs[declev].resize(max_var+1, MemData());
  }



  for (auto vt = comp_vars.begin(); vt != comp_vars.end(); vt++) {
    const auto v = *vt;
    SLOW_DEBUG_DO(assert(is_unknown(v)));

    /* cout << "cur declev is: " << declev << endl; */
    /* for(int32_t k = last_seen[v]; k >= std::min(var_data[v].dirty_lev, declev); k--) { */
    int32_t k = std::min(var_data[v].dirty_lev, declev);
    if (last_seen[v] >= k) {
      int32_t d = std::max(k, 0);
      /* unif_occ[v].resize(long_sz_declevs[d][v].sz); */
      holder.resize(v, long_sz_declevs[d][v].sz);
      /* std::sort(unif_occ[v].begin(), unif_occ[v].end()); */
      /* cout << "resetting size of occ[v " << v << "] to " << unif_occ[v].size() << " stamp:" << stamp << " d: " << d << endl; */
      /* else cout << "not resetting v " << v << " stamp does not match. stamp:" << stamp << " d: " << d << endl; */
    }
    var_data[v].dirty_lev = INT_MAX;
    /* for(uint32_t v = 1; v < max_var+1; v++) { */
      /* std::sort(unif_occ[v].begin(), unif_occ[v].end()); */
      /* cout << "Now v " << v << " contents are: "; */
      /* for(const auto& d: unif_occ[v]) cout <<  d << " , "; */
      /* cout << endl; */
    /* } */

    //traverse binary clauses
    for(const auto& p: unif_occ_bin[v]) {
      if (manage_occ_of(p)) {
        if (is_unknown(p)) { bump_freq_score(p); bump_freq_score(v); }
      }
    }

    // traverse long clauses
    /* cout << "going through v " << v << endl; */
    if (declev != 0) {
      long_sz_declevs[declev][v] = MemData(holder.size(v));
      /* cout << "Remembering v " << v << " size: " << unif_occ[v].size() << " stamp: " << stamp << " declev: " << declev << endl; */
    }

    last_seen[v] = declev;
    uint32_t i = 0;
    while (i < holder.size(v)) {
      ClData& d = holder.begin(v)[i];
      if (archetype.clause_sat(d.id)) goto satl2;
      if (archetype.clause_unvisited_in_sup_comp(d.id)) {
        bool sat = false;
        if (d.tri) {
          Lit l1 = d.get_lit1();
          Lit l2 = d.get_lit2();
          sat = is_true(l1) || is_true(l2);
          if (!sat) {
            if (is_unknown(l1) && !archetype.var_nil(l1.var())) {
              bump_freq_score(l1.var());
              manage_occ_of(l1.var());
            }
            if (is_unknown(l2) && !archetype.var_nil(l1.var())) {
              bump_freq_score(l2.var());
              manage_occ_of(l2.var());
            }
            bump_freq_score(v);
            archetype.set_clause_visited(d.id);
          } else {
            goto satl;
          }
        } else {
          sat = is_true(d.blk_lit);
          if (!sat) sat = search_clause(d, long_clauses_data.data()+d.off);
          if (sat) goto satl;
        }
        i++;
      } else i++;
      continue;

      satl:
      archetype.set_clause_sat(d.id);
      satl2:
      ClData tmp = holder.begin(v)[i];
      holder.begin(v)[i] = holder.back(v);
      holder.back(v) = tmp;
      holder.pop_back(v);
      /* cout << "shrinking size of occ[v " << v << "] to " << unif_occ[v].size() << endl; */
    }
  }

  debug_print(COLWHT "-> Went through all bin/tri/long and now comp_vars is "
      << comp_vars.size() << " long");
}
