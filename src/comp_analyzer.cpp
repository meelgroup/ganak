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
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include "chibihash64.h"

using namespace GanakInt;

std::ostream& operator<<(std::ostream& os, const ClData& d)
{
  os << "[id: " << d.id << " off: " << d.off << "]";
  /* os << "id: " << d.id; */
  return os;
}

// Builds occ lists and sets things up, Done exactly ONCE for a whole counting run
// this sets up unif_occ
void CompAnalyzer::initialize(
    const LiteralIndexedVector<LitWatchList> & watches, // binary clauses
    ClauseAllocator const* alloc, const vector<ClauseOfs>& _long_irred_cls) // longer-than-2-long clauses
{
  watches_ = &watches;
  max_var = watches.end_lit().var() - 1;
  comp_vars.reserve(max_var + 1);
  var_freq_scores.resize(max_var + 1, 0);
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
  vector<vector<ClData>> unif_occ_long(n);
  long_clauses_data.clear();
  long_clauses_data.push_back(SENTINEL_LIT); // MUST start with a sentinel!
  for (const auto& off: long_irred_cls) {
    const Clause& cl = *alloc->ptr(off);
    assert(cl.size() > 2);
    const uint32_t long_cl_off = long_clauses_data.size();
    if (cl.size() > 3) {
      Lit const blk_lit = cl[cl.size()/2];
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
  vector<vector<uint32_t>> unif_occ_bin(n);
  vector<uint32_t> tmp2;
  for (uint32_t v = 1; v < n; v++) {
    tmp2.clear();
    for(bool const sign : {false, true}) {
      for (const auto& bincl: watches[Lit(v, sign)].binaries) {
        if (bincl.irred()) tmp2.push_back(bincl.lit().var());
      }
    }
    // No duplicates, please -- note we converted to VARs so it maybe unique in LIT but not in VAR
    std::sort(tmp2.begin(), tmp2.end());
    tmp2.erase(std::unique(tmp2.begin(), tmp2.end()), tmp2.end());

    unif_occ_bin[v] = tmp2;
  }

  // fill holder
  assert(unif_occ_bin.size() == unif_occ_long.size());
  assert(unif_occ_bin.size() == n);

  size_t const total_sz = hstride * n
    + std::accumulate(unif_occ_long.begin(), unif_occ_long.end(), size_t{0},
        [](size_t acc, const auto& u) { return acc + u.size() * (sizeof(ClData)/sizeof(uint32_t)); })
    + std::accumulate(unif_occ_bin.begin(), unif_occ_bin.end(), size_t{0},
        [](size_t acc, const auto& u) { return acc + u.size(); });
  holder.data = std::make_unique<uint32_t[]>(total_sz);
  uint32_t* const data = holder.data.get();
  uint32_t* data_start = data + n*hstride;

  for(uint32_t v = 0; v < n; v++) {
    // fill bins
    const auto& u_bins = unif_occ_bin[v];
    holder.size_bin(v) = u_bins.size();
    holder.orig_size_bin(v) = u_bins.size();
    uint32_t offs = data_start - data;
    holder.data[v*hstride+holder.offset] = offs;
    assert(offs <= total_sz);
    if (!u_bins.empty()) {
      memcpy(data_start, u_bins.data(), u_bins.size()*sizeof(uint32_t));
      data_start += u_bins.size();
    }

    // fill longs
    const auto& u_longs = unif_occ_long[v];
    holder.orig_size_long(v) = u_longs.size();
    holder.size_long(v) = u_longs.size();
    offs = data_start - data;
    holder.data[v*hstride+holder.offset+3] = offs;
    assert(offs <= total_sz);
    if (!u_longs.empty()) {
      memcpy(data_start, u_longs.data(), u_longs.size()*sizeof(ClData));
      data_start += u_longs.size()*(sizeof(ClData)/sizeof(uint32_t));
    }
    holder.set_tstamp(v, 0);
    holder.set_lev(v, 0);
  }
  assert(data_start == data + total_sz);

  // check bins
  for(uint32_t v = 0; v < unif_occ_bin.size(); v++) {
    assert(unif_occ_bin[v].size() == holder.size_bin(v));
    assert(std::equal(unif_occ_bin[v].begin(), unif_occ_bin[v].end(), holder.begin_bin(v)));
  }

  // check longs
  for(uint32_t v = 0; v < unif_occ_long.size(); v++) {
    assert(unif_occ_long[v].size() == holder.size_long(v));
    assert(std::equal(unif_occ_long[v].begin(), unif_occ_long[v].end(), holder.begin_long(v)));
  }

  debug_print(COLBLBACK "Built unified link list in CompAnalyzer::initialize.");
}

// returns true, iff the comp found is non-trivial
bool CompAnalyzer::explore_comp(const uint32_t v, const uint32_t sup_comp_long_cls, const uint32_t sup_comp_bin_cls) {
  SLOW_DEBUG_DO(assert(archetype.var_unvisited_in_sup_comp(v)));
  record_comp(v, sup_comp_long_cls, sup_comp_bin_cls); // sets up the component that "v" is in

  if (comp_vars.size() == 1) {
    debug_print("in " <<  __FUNCTION__ << " with single var: " <<  v);
    if (v >= counter->get_indep_support_end()) {
      SLOW_DEBUG_DO(
        if (v < counter->get_opt_indep_support_end()) {
            counter->check_trail(true, true);
            counter->check_opt_sampling_determined();
            debug_print("This is a VERY interesting phenomenon."
               << " We MUST be in a situation where we are UNSAT, but the solver hasn't yet determined this"
               << " We simply multiply by one. It'll be all undone anyway, as unsat MUST be detected later");
            counter->check_current_state_unsat();
        }
      );
      archetype.stack_level().include_solution(counter->get_fg()->one());
    } else {
      if (counter->weighted()) archetype.stack_level().include_solution(counter->get_weight(v));
      else archetype.stack_level().include_solution(counter->get_two());
    }
    archetype.set_var_clear(v);
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

  debug_print(COLWHT "We are NOW going through all binary/tri/long clauses "
      "recursively and put into search_stack_ all the variables that are connected to var: " << var);
  stats.comps_recorded++;

  for (uint32_t i = 0; i < comp_vars.size(); i++) {
    const auto v = comp_vars[i];
    SLOW_DEBUG_DO(assert(is_unknown(v)));
    analyze_verb(
      debug_print("-----------------------");
      debug_print("record v: " << v << " start");
      debug_print("holder.lev(v): " << holder.lev(v));
      debug_print("holder.tstamp(v): " << holder.tstamp(v));
      debug_print("counter->dec_level(): " << counter->dec_level());
      debug_print("counter->get_tstamp(holder.lev(v)): " << counter->get_tstamp(holder.lev(v)));
      counter->print_trail());

    bool reset = false;
    analyze_verb(debug_print("v: " << v << " holder.lev(v): " << holder.lev(v)
      << " holder.tstamp(v): " << holder.tstamp(v)
      << " counter->dec_lev(): " << counter->dec_level()
      << " counter->get_tstamp(holder.lev(v))): " << counter->get_tstamp(holder.lev(v))));
    if (holder.tstamp(v) < counter->get_tstamp(holder.lev(v))) {
      /* holder.size_bin(v) = holder.orig_size_bin(v); */
      holder.size_long(v) = holder.orig_size_long(v);
      stats.comps_reset++;
      reset = true;
      analyze_verb(debug_print("analyze RESET"));
    } else {
      analyze_verb(debug_print("analyze NORESET"));
      stats.comps_non_reset++;
    }
    bool const update = (counter->last_dec_candidates > conf.analyze_cand_update) || reset;
    if (update) {
      holder.set_tstamp(v, counter->get_tstamp());
      holder.set_lev(v, counter->dec_level());
      analyze_verb(debug_print("analyze tstamp UPDATED. v: " << v << " holder.lev(v): " << holder.lev(v)
          << " holder.tstamp(v): " << holder.tstamp(v)));
    } else {
      analyze_verb(debug_print("analyze tstamp REMAIN. v: " << v << " holder.lev(v): " << holder.lev(v)
          << " holder.tstamp(v): " << holder.tstamp(v)));
    }

    SLOW_DEBUG_DO(
      // checks that bins between size and orig_size are all satisfied
      uint32_t* bins = holder.begin_bin(v);
      uint32_t* bins_end2 = bins+holder.size_bin(v);
      uint32_t* bins_end3 = bins+holder.orig_size_bin(v);
      while (bins_end2 != bins_end3) {
        const uint32_t v2 = *bins_end2;
        if (is_unknown(v2)) {
          cerr << "ERROR: bin clause var: " << v2 << " unknown, but we thought it's set (and in the bin, true)!" << endl;
          release_assert(false);
        }
        bins_end2++;
      }
      if (reset) assert(holder.size_bin(v) == holder.orig_size_bin(v));
    );

    if (sup_comp_bin_cls != archetype.num_bin_cls) {
      // we have not seen all binary clauses
      // traverse binary clauses
      uint32_t* bins = holder.begin_bin(v);
      uint32_t const* bins_end = bins + holder.size_bin(v);
      while(bins != bins_end) {
        uint32_t const v2 = *(bins++);
        // v2 must be true or unknown, because if it's false, this variable would be TRUE, and that' not the case
        manage_occ_of(v2);
        if (is_unknown(v2)) {
          archetype.num_bin_cls++;
          bump_freq_score(v2);
          bump_freq_score(v);
        } else {
          /* if (update) { */
          /*   // it's satisfied */
          /*   bins--; */
          /*   bins_end--; */
          /*   std::swap(*bins, *bins_end); */
          /*   holder.size_bin(v)--; */
          /*   analyze_verb(verb_debug("analyze remove bin, var: "<< *bins_end)); */
          /* } */
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
            archetype.set_cl_clear(d.id);
            sat = true;
            goto end_sat;
          } else {
            bump_freq_score(v);
            manage_occ_and_score_of(l1);
            manage_occ_and_score_of(l2);
            archetype.set_clause_visited(d.id);
          }
        } else continue;
      } else {
        if (archetype.clause_unvisited_in_sup_comp(d.id)) {
          if (is_true(d.blk_lit)) {
            archetype.set_cl_clear(d.id);
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
      if (update) {
        longs--;
        longs_end--;
        analyze_verb(
          const ClData& d2 = *longs;
          cout << "analyze remove SAT clause id: " << d2.id << " cl:";
          for (auto it_l = long_clauses_data.data()+d2.off; *it_l != SENTINEL_LIT; it_l++) {
            cout << *it_l << " ";
          } cout << endl);
        std::swap(*longs, *longs_end);
        holder.size_long(v)--;
      }
    }
    /* cout << "AFTER holder.size_long(v): " << holder.size_long(v) << endl; */
    /* cout << "AFTER holder.orig_size_long(v): " << holder.orig_size_long(v) << endl; */
  }
  debug_print(COLWHT "-> Went through all bin/tri/long and now comp_vars is "
      << comp_vars.size() << " long");
}

// Computes Weisfeiler-Lehman canonical component information for cache lookup.
//
// Approach:
//   1. Build a clause->variable mapping restricted to this component's vars/clauses.
//   2. Compute an initial "color" for each variable: (long_degree, binary_degree, is_indep).
//   3. Run one round of WL refinement: each variable's new color incorporates the sorted
//      multiset of its long-clause neighbours' initial colors.
//   4. Sort variables by WL color (original var_id as tiebreaker) to get canonical order.
//   5. Express every clause (long + binary) in terms of canonical variable indices and sort.
//   6. Hash the resulting canonical clause list to produce a structure-invariant cache key.
//
// The result is invariant to variable/clause ID renaming: two structurally isomorphic
// components will produce identical sorted_canon_clauses and the same hash, enabling a
// cache hit even if they involve completely different variable numberings.
CanonInfo CompAnalyzer::compute_canon_info(const Comp& comp,
                                           uint64_t hash_seed,
                                           uint32_t threshold) const {
  CanonInfo info;
  const uint32_t n = comp.nVars();
  if (threshold == 0 || n > threshold || n == 0) return info;

  // --- Build membership lookups ---

  // comp clause set for O(1) membership test
  std::unordered_set<uint32_t> comp_clause_set;
  comp_clause_set.reserve(comp.num_long_cls() * 2);
  for (auto it = comp.cls_begin(); *it != sentinel; ++it) comp_clause_set.insert(*it);

  // var_id -> position in comp (0..n-1)
  std::unordered_map<uint32_t, uint32_t> var_to_pos;
  var_to_pos.reserve(n * 2);
  for (uint32_t i = 0; i < n; ++i) var_to_pos[comp.vars_begin()[i]] = i;

  // --- Compute degree-based initial color per variable (position) ---

  // long_deg[i]  = number of comp clauses containing position-i variable
  // bin_deg[i]   = number of comp variables that position-i variable shares a binary clause with
  vector<uint32_t> long_deg(n, 0);
  vector<uint32_t> bin_deg(n, 0);

  // clause_id -> list of comp-variable positions that appear in it (for WL neighbor graph)
  std::unordered_map<uint32_t, vector<uint32_t>> clause_to_pos;
  clause_to_pos.reserve(comp.num_long_cls() * 2);

  // clause_id -> all Lits with signs (for polarity-aware canonical form)
  // For ternary clauses: accumulated across all 3 variable visits.
  // For long clauses: read from long_clauses_data on first encounter.
  std::unordered_map<uint32_t, vector<Lit>> clause_all_lits;
  clause_all_lits.reserve(comp.num_long_cls() * 2);

  for (uint32_t i = 0; i < n; ++i) {
    const uint32_t v = comp.vars_begin()[i];

    // Long clauses: iterate over all long clauses v appears in.
    const ClData* longs     = holder.begin_long(v);
    const ClData* longs_end = longs + holder.orig_size_long(v);
    for (const ClData* d = longs; d != longs_end; ++d) {
      if (!comp_clause_set.count(d->id)) continue;
      clause_to_pos[d->id].push_back(i);
      ++long_deg[i];

      auto& lits = clause_all_lits[d->id];
      if (d->id < max_tri_clid) {
        // Ternary clause: contribute the 2 "other" literals from this visit.
        // After all 3 vars are visited, lits will contain all 3 signed literals.
        const Lit l1 = d->get_lit1();
        const Lit l2 = d->get_lit2();
        if (std::find(lits.begin(), lits.end(), l1) == lits.end()) lits.push_back(l1);
        if (std::find(lits.begin(), lits.end(), l2) == lits.end()) lits.push_back(l2);
      } else {
        // Long clause: read all signed literals from the literal pool on first encounter.
        if (lits.empty()) {
          const Lit* start = long_clauses_data.data() + d->off;
          for (const Lit* it_l = start; *it_l != SENTINEL_LIT; ++it_l) lits.push_back(*it_l);
        }
      }
    }

    // Binary degree (polarity-blind, used for WL initial color only).
    const uint32_t* bins = holder.begin_bin(v);
    for (uint32_t j = 0; j < holder.orig_size_bin(v); ++j) {
      if (var_to_pos.count(bins[j])) ++bin_deg[i];
    }
  }

  // --- Initial WL color per position ---
  using Color3 = std::tuple<uint32_t, uint32_t, uint32_t>;
  vector<Color3> init_color(n);
  for (uint32_t i = 0; i < n; ++i) {
    const bool is_indep = (comp.vars_begin()[i] < indep_support_end);
    init_color[i] = {long_deg[i], bin_deg[i], static_cast<uint32_t>(is_indep)};
  }
  VERBOSE_DEBUG_DO(
    cout << "WL canon: nVars=" << n
         << " nLongCls=" << clause_all_lits.size()
         << " initial colors (var longdeg bindeg isindep):";
    for (uint32_t i = 0; i < n; ++i)
      cout << " [" << comp.vars_begin()[i] << " "
           << long_deg[i] << " " << bin_deg[i] << " "
           << (comp.vars_begin()[i] < indep_support_end ? 1 : 0) << "]";
    cout << endl;
  );

  // --- Build long-clause neighbour adjacency (for WL round) ---
  // cl_neighbors[i] = positions of variables that share a long clause with position i
  vector<vector<uint32_t>> cl_neighbors(n);
  for (auto& [cl_id, positions] : clause_to_pos) {
    for (uint32_t u : positions)
      for (uint32_t w : positions)
        if (u != w) cl_neighbors[u].push_back(w);
  }

  // --- One round of WL refinement ---
  // wl1[i] = hash of (init_color[i], sorted multiset of init_colors of clause-neighbours)
  vector<uint64_t> wl1(n);
  for (uint32_t i = 0; i < n; ++i) {
    vector<Color3> ncolors;
    ncolors.reserve(cl_neighbors[i].size());
    for (uint32_t j : cl_neighbors[i]) ncolors.push_back(init_color[j]);
    sort(ncolors.begin(), ncolors.end());

    // Mix into a 64-bit hash using simple polynomial mixing
    uint64_t h = (static_cast<uint64_t>(get<0>(init_color[i])) * 2654435761ULL)
               ^ (static_cast<uint64_t>(get<1>(init_color[i])) * 40503ULL)
               ^ (static_cast<uint64_t>(get<2>(init_color[i])) * 2246822519ULL);
    for (const auto& [ld, bd, indp] : ncolors) {
      h ^= (h >> 16) * 0x45d9f3bULL;
      h += (static_cast<uint64_t>(ld) * 2654435761ULL)
         ^ (static_cast<uint64_t>(bd) * 40503ULL)
         ^ (static_cast<uint64_t>(indp));
    }
    wl1[i] = h;
  }

  VERBOSE_DEBUG_DO(
    cout << "WL canon: wl1 colors (var wl1hash):";
    for (uint32_t i = 0; i < n; ++i)
      cout << " [" << comp.vars_begin()[i] << " 0x" << std::hex << wl1[i] << std::dec << "]";
    cout << endl;
  );

  // --- Sort variables by (wl1, init_color, var_id) to get canonical order ---
  vector<uint32_t> perm(n);
  iota(perm.begin(), perm.end(), 0);
  sort(perm.begin(), perm.end(), [&](uint32_t a, uint32_t b) {
    if (wl1[a] != wl1[b]) return wl1[a] < wl1[b];
    if (init_color[a] != init_color[b]) return init_color[a] < init_color[b];
    return comp.vars_begin()[a] < comp.vars_begin()[b]; // stable tiebreak
  });

  // canon_vars[i] = original var_id at canonical position i
  info.canon_vars.resize(n);
  for (uint32_t i = 0; i < n; ++i) info.canon_vars[i] = comp.vars_begin()[perm[i]];

  // canon_is_indep[i] = 1 if canonical position i is in the independent support, else 0.
  // Must be included in the hash/equality data so that two structurally isomorphic
  // components that differ only in their indep-support assignments are not confused.
  info.canon_is_indep.resize(n);
  for (uint32_t i = 0; i < n; ++i)
    info.canon_is_indep[i] = static_cast<uint32_t>(comp.vars_begin()[perm[i]] < indep_support_end);

  // orig_pos -> canonical index
  vector<uint32_t> orig_to_canon(n);
  for (uint32_t i = 0; i < n; ++i) orig_to_canon[perm[i]] = i;

  // --- Build polarity-aware canonical clause representations ---
  // Canonical literal index: 2 * canon_pos + (uint32_t)lit.sign()
  //   (sign()=false → negative literal, sign()=true → positive literal)
  //
  // Long/ternary: from clause_all_lits, filter to in-component (unknown) vars.
  info.sorted_canon_clauses.reserve(clause_all_lits.size() + n);
  for (auto& [cl_id, lits] : clause_all_lits) {
    vector<uint32_t> cv;
    cv.reserve(lits.size());
    for (const Lit l : lits) {
      auto it = var_to_pos.find(l.var());
      if (it == var_to_pos.end()) continue; // satisfied/false lit, skip
      cv.push_back(2 * orig_to_canon[it->second] + static_cast<uint32_t>(l.sign()));
    }
    if (cv.size() < 2) continue; // should not happen for in-comp clauses
    sort(cv.begin(), cv.end());
    info.sorted_canon_clauses.push_back(std::move(cv));
  }

  // Binary clauses: use watches_ directly to get full literal polarities.
  // watches_[Lit(v,s)] contains binary clauses of the form (Lit(v,s) ∨ bincl.lit()).
  // Deduplicate via a 64-bit key (packed canonical lit-index pair, lo<<32|hi).
  std::unordered_set<uint64_t> seen_bin;
  seen_bin.reserve(n * 4);
  for (uint32_t pos_i = 0; pos_i < n; ++pos_i) {
    const uint32_t v = comp.vars_begin()[pos_i];
    for (const bool s : {false, true}) {
      for (const auto& bincl : (*watches_)[Lit(v, s)].binaries) {
        if (!bincl.irred()) continue;
        const Lit other = bincl.lit();
        auto it = var_to_pos.find(other.var());
        if (it == var_to_pos.end()) continue; // other end not in comp
        const uint32_t cli = 2 * orig_to_canon[pos_i] + static_cast<uint32_t>(s);
        const uint32_t clj = 2 * orig_to_canon[it->second] + static_cast<uint32_t>(other.sign());
        const uint32_t lo = std::min(cli, clj);
        const uint32_t hi = std::max(cli, clj);
        const uint64_t key = (static_cast<uint64_t>(lo) << 32) | hi;
        if (seen_bin.insert(key).second)
          info.sorted_canon_clauses.push_back({lo, hi});
      }
    }
  }

  // --- Sort all canonical clauses lexicographically ---
  sort(info.sorted_canon_clauses.begin(), info.sorted_canon_clauses.end());

  // --- Compute structural hash of (nVars, canonical clauses, is_indep profile) ---
  // Encoding: [n, n_total_clauses, for each clause: (size, v0, v1, ...), is_indep[0..n-1]]
  // The is_indep profile distinguishes components whose clause structures are isomorphic
  // but differ in which canonical positions belong to the independent (projection) support.
  vector<uint32_t> hdata;
  hdata.reserve(2 + info.sorted_canon_clauses.size() * 4 + n);
  hdata.push_back(n);
  hdata.push_back(static_cast<uint32_t>(info.sorted_canon_clauses.size()));
  for (const auto& cv : info.sorted_canon_clauses) {
    hdata.push_back(static_cast<uint32_t>(cv.size()));
    for (uint32_t idx : cv) hdata.push_back(idx);
  }
  for (uint32_t i = 0; i < n; ++i) hdata.push_back(info.canon_is_indep[i]);
  info.hash = chibihash64(hdata.data(), hdata.size() * sizeof(uint32_t), hash_seed);

  VERBOSE_DEBUG_DO(
    cout << "WL canon: final hash=0x" << std::hex << info.hash << std::dec
         << " nclauses=" << info.sorted_canon_clauses.size()
         << " canonical clauses:";
    for (const auto& cv : info.sorted_canon_clauses) {
      cout << " (";
      for (uint32_t idx = 0; idx < cv.size(); ++idx) {
        if (idx) cout << ",";
        cout << cv[idx];
      }
      cout << ")";
    }
    cout << endl;
  );

  info.valid = true;
  return info;
}

// There is exactly ONE of these
CompAnalyzer::CompAnalyzer(
    const LiteralIndexedVector<TriValue> & lit_values,
    Counter* _counter) :
      stats(_counter->get_stats()),
      values(lit_values),
      conf(_counter->get_conf()),
      indep_support_end(_counter->get_indep_support_end()),
      counter(_counter)
{}
