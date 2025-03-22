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
#include <climits>
#include <cstdint>
#include <cstdint>
#include <memory>

using namespace GanakInt;

inline std::ostream& operator<<(std::ostream& os, const ClData& d)
{
  os << "[id: " << d.id << " off: " << d.off << "]";
  /* os << "id: " << d.id; */
  return os;
}

// Builds occ lists and sets things up, Done exactly ONCE for a whole counting runkk
// this sets up unif_occ
void CompAnalyzer::initialize(
    const LiteralIndexedVector<LitWatchList> & watches, // binary clauses
    const ClauseAllocator* alloc, const vector<ClauseOfs>& _long_irred_cls) // longer-than-2-long clauses
{
  max_var = watches.end_lit().var() - 1;
  comp_vars.reserve(max_var + 1);
  var_freq_scores.resize(max_var + 1, 0);
  vector<vector<ClauseOfs>> occs(max_var + 1);
  const uint32_t n = max_var+1;

  debug_print(COLBLBACK "Building occ list in CompAnalyzer::initialize...");

  auto mysorter = [&] (ClauseOfs a1, ClauseOfs b1) {
    const Clause& a = *alloc->ptr(a1);
    const Clause& b = *alloc->ptr(b1);
    return a.size() < b.size();
  };
  auto long_irred_cls = _long_irred_cls;
  std::sort(long_irred_cls.begin(), long_irred_cls.end(), mysorter);

  max_clid = 1;
  max_tri_clid = 1;
  vector<vector<ClData>> unif_occ_long;
  unif_occ_long.clear();
  unif_occ_long.resize(n);
  long_clauses_data.clear();
  long_clauses_data.push_back(SENTINEL_LIT); // MUST start with a sentinel!
  vector<uint32_t> tmp;
  for (const auto& off: long_irred_cls) {
    const Clause& cl = *alloc->ptr(off);
    assert(cl.size() > 2);
    const uint32_t long_cl_off = long_clauses_data.size();
    if (cl.size() > 3) {
      Lit blk_lit = cl[cl.size()/2];
      for(const auto&l: cl) long_clauses_data.push_back(l);
      long_clauses_data.push_back(SENTINEL_LIT);

      for(const auto& l: cl) {
        const uint32_t var = l.var();
        assert(var < n);
        ClData d;
        d.id = max_clid;
        d.blk_lit = blk_lit;
        d.off = long_cl_off;
        unif_occ_long[var].push_back(d);
      }
    } else {
      assert(cl.size() == 3);
      for(const auto& l: cl) {
        uint32_t at = 0;
        Lit lits[2];
        for(const auto&l2: cl) if (l.var() != l2.var()) lits[at++] = l2;
        assert(at == 2);
        ClData d;
        d.id = max_clid;
        d.blk_lit = lits[0];
        d.off = lits[1].raw();
        unif_occ_long[l.var()].push_back(d);
      }
      assert(max_tri_clid == max_clid && "it's sorted by clause size!");
      max_tri_clid++;
    }
    /* cout << "cl id: " << max_clid << " size: " << cl.size() << " cl: "; */
    /* for(const auto& l: cl) cout << " " << l; */
    /* cout << endl; */
    max_clid++;
  }
  /* cout << "max clid: " << max_clid << " max_tri_clid: " << max_tri_clid << endl;; */
  debug_print(COLBLBACK "Built occ list in CompAnalyzer::initialize.");

  archetype.init_data(max_var, max_clid);
  debug_print(COLBLBACK "Building unified link list in CompAnalyzer::initialize...");


  // data for binary clauses
  vector<vector<uint32_t>> unif_occ_bin;
  unif_occ_bin.clear();
  unif_occ_bin.resize(n);
  vector<uint32_t> tmp2;
  for (uint32_t v = 1; v < n; v++) {
    tmp2.clear();
    for(uint32_t i = 0; i < 2; i++) {
      for (const auto& bincl: watches[Lit(v, i)].binaries) {
        if (bincl.irred()) tmp2.push_back(bincl.lit().var());
      }
    }
    // No duplicates, please -- note we converted to VARs so it maybe unique in LIT but not in VAR
    set<uint32_t> s( tmp2.begin(), tmp2.end());
    tmp2.assign( s.begin(), s.end() );

    unif_occ_bin[v].clear();
    unif_occ_bin[v].resize(tmp2.size());
    for(uint32_t i = 0; i < tmp2.size(); i++) unif_occ_bin[v][i] = tmp2[i];
  }

  if (true) {
    // fill holder
    assert(unif_occ_bin.size() == unif_occ_long.size());
    assert(unif_occ_bin.size() == n);

    uint32_t total_sz = 0;
    for(const auto& u: unif_occ_long) total_sz += u.size()*(sizeof(ClData)/sizeof(uint32_t));
    for(const auto& u: unif_occ_bin) total_sz += u.size();
    total_sz += hstride*n;
    uint32_t* data = new uint32_t[total_sz];
    holder.data = data;
    uint32_t* data_start = holder.data + n*hstride;

    for(uint32_t v = 0; v < n; v++) {
      // fill bins
      const auto& u_bins = unif_occ_bin[v];
      holder.size_bin(v) = u_bins.size();
      holder.orig_size_bin(v) = u_bins.size();
      uint32_t offs = data_start - holder.data;
      holder.data[v*hstride+holder.offset] = offs;
      assert(offs <= total_sz);
      memcpy(data_start, u_bins.data(), u_bins.size()*sizeof(uint32_t));
      data_start += u_bins.size();

      // fill longs
      const auto& u_longs = unif_occ_long[v];
      holder.orig_size_long(v) = u_longs.size();
      holder.size_long(v) = u_longs.size();
      offs = data_start - holder.data;
      holder.data[v*hstride+holder.offset+3] = offs;
      assert(offs <= total_sz);
      memcpy(data_start, u_longs.data(), u_longs.size()*sizeof(ClData));
      data_start += u_longs.size()*(sizeof(ClData)/sizeof(uint32_t));
      holder.tstamp(v) = 0;
      holder.lev(v) = 0;
    }
    assert(data_start == data + total_sz);

    // check bins
    for(uint32_t v = 0; v < unif_occ_bin.size(); v++) {
      assert(unif_occ_bin[v].size() == holder.size_bin(v));
      for(uint32_t i = 0; i < unif_occ_bin[v].size(); i++) {
        assert(unif_occ_bin[v][i] == holder.begin_bin(v)[i]);
      }
    }

    // check longs
    for(uint32_t v = 0; v < unif_occ_long.size(); v++) {
      assert(unif_occ_long[v].size() == holder.size_long(v));
      for(uint32_t i = 0; i < unif_occ_long[v].size(); i++) {
        assert(unif_occ_long[v][i] == holder.begin_long(v)[i]);
      }
    }
  }

  debug_print(COLBLBACK "Built unified link list in CompAnalyzer::initialize.");
}

// returns true, iff the comp found is non-trivial
bool CompAnalyzer::explore_comp(const uint32_t v, const uint32_t sup_comp_long_cls, const uint32_t sup_comp_bin_cls) {
  SLOW_DEBUG_DO(assert(archetype.var_unvisited_in_sup_comp(v)));
  record_comp(v, sup_comp_long_cls, sup_comp_bin_cls); // sets up the component that "v" is in

  if (comp_vars.size() == 1) {
    debug_print("in " <<  __FUNCTION__ << " with single var: " <<  v);
    if (v >= indep_support_end) archetype.stack_level().include_solution(counter->get_fg()->one());
    else {
      if (counter->weighted()) archetype.stack_level().include_solution(counter->get_weight(v));
      else archetype.stack_level().include_solution(counter->get_two());
    }
    archetype.set_var_in_peer_comp(v);
    return false;
  }
  return true;
}

// Each variable knows the level it was visited at, and the stimestamp at the time
// Each level knows the HIGHEST stamp it has been seen
// When checking a var, we go to the level, see the stamp, if it's larger than the stamp of the var,
// we need to reset the size

// Create a component based on variable provided
void CompAnalyzer::record_comp(const uint32_t var, const uint32_t sup_comp_long_cls, const uint32_t sup_comp_bin_cls) {
  SLOW_DEBUG_DO(assert(is_unknown(var)));
  comp_vars.clear();
  comp_vars.push_back(var);
  archetype.set_var_visited(var);
  counter->set_stamp();

  debug_print(COLWHT "We are NOW going through all binary/tri/long clauses "
      "recursively and put into search_stack_ all the variables that are connected to var: " << var);
  if (comps_recorded % (128*1024*1024) == 0 && reset_comps+non_reset_comps > 0) {
    verb_print(1, "[record_comp] Tot recorded " << (comps_recorded)/1000 << "K non-reset ratio: "
        << (double)(non_reset_comps)/(double)(reset_comps+non_reset_comps));
  }
  comps_recorded++;

  for (auto vt = comp_vars.begin(); vt != comp_vars.end(); vt++) {
    const auto v = *vt;
    SLOW_DEBUG_DO(assert(is_unknown(v)));

#ifdef ANALYZE_DEBUG
    cout << "-----------------------" << endl;
    cout << "record v: " << v << " start" << endl;
    cout << "holder.lev(v): " << holder.lev(v) << endl;
    cout << "holder.tstamp(v): " << holder.tstamp(v) << endl;
    cout << "counter->dec_level(): " << counter->dec_level() << endl;
    cout << "counter->get_tstamp(holder.lev(v)): " << counter->get_tstamp(holder.lev(v)) << endl;
    counter->print_trail();
#endif

    if (holder.tstamp(v) < counter->get_tstamp(holder.lev(v))) {
      /* holder.size_bin(v) = holder.orig_size_bin(v); */
      holder.size_long(v) = holder.orig_size_long(v);
      reset_comps++;
    } else {
      non_reset_comps++;
    }
    holder.lev(v) = counter->dec_level();
    holder.tstamp(v) = counter->get_tstamp();
#ifdef ANALYZE_DEBUG
    cout << "AFTER holder.lev(v): " << holder.lev(v) << endl;
    cout << "AFTER holder.tstamp(v): " << holder.tstamp(v) << endl;
#endif

    if (sup_comp_bin_cls == archetype.num_bin_cls) {
      // we have seen all bin clauses
      if (sup_comp_long_cls == archetype.num_long_cls) {
        // we have seen all bin and long clauses
        break;
      }
    } else {
      //traverse binary clauses
      auto bins = holder.begin_bin(v);
      auto bins_end = bins + holder.size_bin(v);
      while(bins != bins_end) {
        uint32_t v2 = *(bins++);
        // v2 must be true or unknown, because if it's false, this variable would be TRUE, and that' not the case
        manage_occ_of(v2);
        if (is_unknown(v2)) {
          archetype.num_bin_cls++;
          bump_freq_score(v2);
          bump_freq_score(v);
        } else {
          // it's true
          /* bins--; */
          /* bins_end--; */
          /* std::swap(*bins, *bins_end); */
          /* holder.size_bin(v)--; */
        }
      }
    }
    SLOW_DEBUG_DO(assert(archetype.num_bin_cls <= sup_comp_bin_cls));

    if (sup_comp_long_cls == archetype.num_long_cls) {
      // we have seen all long clauses
      continue;
    }

#ifdef SLOW_DEBUG
    {
      // checks that longs between size and orig_size are all satisfied
      /* cout << "holder.size_long(v): " << holder.size_long(v) << endl; */
      /* cout << "holder.orig_size_long(v): " << holder.orig_size_long(v) << endl; */
      ClData* longs = holder.begin_long(v);
      ClData* longs_end2 = longs+holder.size_long(v);
      ClData* longs_end3 = longs+holder.orig_size_long(v);
      while (longs_end2 != longs_end3) {
        const ClData& d = *longs_end2;
        const Lit* start = long_clauses_data.data()+d.off;
        if (d.id < max_tri_clid) {
          const Lit l1 = d.get_lit1();
          const Lit l2 = d.get_lit2();
          assert(is_true(l1) || is_true(l2));
        } else {
          bool sat = false;
          for (auto it_l = start; *it_l != SENTINEL_LIT; it_l++) {
            if (is_true(*it_l)) sat = true;
          }
          if (!sat) {
            cout << "long clause id: " << d.id << " not satisfied: ";
            for (auto it_l = start; *it_l != SENTINEL_LIT; it_l++) {
              cout << *it_l << " ";
            }
            cout << endl;
            assert(sat);
          }
        }
        longs_end2++;
      }
    }
#endif
    ClData* longs = holder.begin_long(v);
    ClData* longs_end = longs+holder.size_long(v);
    while (longs != longs_end) {
      SLOW_DEBUG_DO(assert(archetype.num_long_cls <= sup_comp_long_cls));
      ClData& d = *longs;
      longs++;
      bool sat = false;
      if (d.id < max_tri_clid) {
        // traverse ternary clauses
        if (archetype.clause_unvisited_in_sup_comp(d.id)) {
          const Lit l1 = d.get_lit1();
          const Lit l2 = d.get_lit2();
          if (is_true(l1) || is_true(l2)) {
            archetype.clear_cl(d.id);
            sat = true;
            /* goto end_sat; */
          } else {
            bump_freq_score(v);
            manage_occ_and_score_of(l1.var());
            manage_occ_and_score_of(l2.var());
            archetype.set_clause_visited(d.id);
          }
        } else continue;
      } else {
        if (archetype.clause_unvisited_in_sup_comp(d.id)) {
          if (is_true(d.blk_lit)) {
            archetype.clear_cl(d.id);
            sat = true;
            goto end_sat;
          }
          const Lit* start = long_clauses_data.data()+d.off;
          sat = search_clause(d, start);
          if (sat) goto end_sat;
        } else continue;
      }
      if (!sat) archetype.num_long_cls++;
      continue;

end_sat:;
      longs--;
      longs_end--;
#ifdef ANALYZE_DEBUG
      const ClData& d2 = *longs;
      cout << "SAT clause id: " << d2.id << " cl:";
      for (auto it_l = long_clauses_data.data()+d2.off; *it_l != SENTINEL_LIT; it_l++) {
        cout << *it_l << " ";
      }
      cout << endl;
#endif
      std::swap(*longs, *longs_end);
      holder.size_long(v)--;
    }
    /* cout << "AFTER holder.size_long(v): " << holder.size_long(v) << endl; */
    /* cout << "AFTER holder.orig_size_long(v): " << holder.orig_size_long(v) << endl; */
  }
  debug_print(COLWHT "-> Went through all bin/tri/long and now comp_vars is "
      << comp_vars.size() << " long");
}

// There is exactly ONE of these
CompAnalyzer::CompAnalyzer(
    const LiteralIndexedVector<TriValue> & lit_values,
    Counter* _counter) :
      values(lit_values),
      conf(_counter->get_conf()),
      indep_support_end(_counter->get_indep_support_end()),
      counter(_counter)
{}
