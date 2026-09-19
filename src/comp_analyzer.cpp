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
      counter->get_compiler().free_var(v);
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
        // assigned vars are never marked unvisited-in-sup-comp, so manage_occ_of()
        // is a no-op for them
        if (is_unknown(v2)) {
          manage_occ_of(v2);
          archetype.num_bin_cls++;
          bump_freq_score(v2, 1+conf.freq_short_bonus);
          bump_freq_score(v, 1+conf.freq_short_bonus);
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
            // v is unknown. If one of the others is false, it is a binary now
            const uint32_t by = (is_unknown(l1) && is_unknown(l2)) ? 1 : 1+conf.freq_short_bonus;
            bump_freq_score(v, by);
            manage_occ_and_score_of(l1, by);
            manage_occ_and_score_of(l2, by);
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

// Iterative Tarjan over the unknown vars of the comp. Binary and ternary
// clauses are var-var edges, a longer clause is a node of its own so that it
// costs its length, not its length squared. Reads the FULL occurrence lists
// (orig_size) and checks satisfiedness itself, so it does not depend on, or
// disturb, the trimmed lists explore_comp() maintains.
uint32_t CompAnalyzer::compute_cut_gains(const Comp& comp) {
  const uint32_t num_nodes = max_var + max_clid + 2;
  if (cut_epoch_of.size() < num_nodes) {
    cut_epoch_of.assign(num_nodes, 0);
    cut_disc.resize(num_nodes);
    cut_low.resize(num_nodes);
    cut_subsz.resize(num_nodes);
    cut_gain.assign(max_var+1, 0);
    cut_sep_sum.resize(max_var+1);
    cut_sep_max.resize(max_var+1);
  }
  cut_epoch++;
  if (cut_epoch == 0) {
    std::fill(cut_epoch_of.begin(), cut_epoch_of.end(), 0);
    cut_epoch = 1;
  }

  uint32_t tot_unknown = 0;
  all_vars_in_comp(comp, vt) {
    const uint32_t root = *vt;
    if (!is_unknown(root)) continue;
    if (cut_epoch_of[root] == cut_epoch) continue;

    uint32_t timer = 0;
    uint32_t root_children = 0;
    cut_tree_vars.clear();
    cut_stack.clear();
    auto enter = [&](const uint32_t node, const uint32_t parent, Lit const* lit_at) {
      cut_epoch_of[node] = cut_epoch;
      cut_disc[node] = cut_low[node] = ++timer;
      if (node <= max_var) {
        cut_subsz[node] = 1;
        cut_sep_sum[node] = 0;
        cut_sep_max[node] = 0;
        cut_gain[node] = 0;
        cut_tree_vars.push_back(node);
      } else cut_subsz[node] = 0;
      cut_stack.push_back(CutFrame{node, parent, 0, 0, lit_at});
    };
    // returns true if w was entered, i.e. the caller must stop and descend
    auto visit = [&](const uint32_t u, const uint32_t parent, const uint32_t w,
        Lit const* lit_at) -> bool {
      if (cut_epoch_of[w] != cut_epoch) { enter(w, u, lit_at); return true; }
      if (w != parent) cut_low[u] = std::min(cut_low[u], cut_disc[w]);
      return false;
    };

    enter(root, 0, nullptr);
    while (!cut_stack.empty()) {
      // NOTE: enter() may reallocate cut_stack, so index, don't hold a ref across it
      const size_t at = cut_stack.size()-1;
      const uint32_t u = cut_stack[at].node;
      const uint32_t parent = cut_stack[at].parent;
      bool descended = false;

      if (u <= max_var) {
        uint32_t const* bins = holder.begin_bin(u);
        const uint32_t nbins = holder.orig_size_bin(u);
        while (!descended && cut_stack[at].i_bin < nbins) {
          const uint32_t v2 = bins[cut_stack[at].i_bin++];
          if (is_unknown(v2)) descended = visit(u, parent, v2, nullptr);
        }
        ClData const* longs = holder.begin_long(u);
        const uint32_t nlongs = holder.orig_size_long(u);
        while (!descended && cut_stack[at].i_long < nlongs) {
          const ClData& d = longs[cut_stack[at].i_long++];
          if (d.id < max_tri_clid) {
            const Lit l1 = d.get_lit1();
            const Lit l2 = d.get_lit2();
            if (is_true(l1) || is_true(l2)) continue;
            // We may descend into l1 and so never look at l2 from here. That
            // is fine: l2 has this clause too, and l1-l2 are connected by it
            if (is_unknown(l1)) descended = visit(u, parent, l1.var(), nullptr);
            if (!descended && is_unknown(l2)) descended = visit(u, parent, l2.var(), nullptr);
            else if (descended && is_unknown(l2)) cut_stack[at].i_long--; // come back for l2
          } else {
            if (is_true(d.blk_lit)) continue;
            const uint32_t cnode = max_var + 1 + d.id;
            if (cut_epoch_of[cnode] != cut_epoch) {
              Lit const* start = long_clauses_data.data()+d.off;
              bool sat = false;
              for (auto it_l = start; *it_l != SENTINEL_LIT; it_l++)
                if (is_true(*it_l)) {sat = true; break;}
              if (sat) {
                // mark it seen with a disc that can never lower anyone's low
                cut_epoch_of[cnode] = cut_epoch;
                cut_disc[cnode] = std::numeric_limits<uint32_t>::max();
                continue;
              }
              enter(cnode, u, start);
              descended = true;
            } else if (cnode != parent) cut_low[u] = std::min(cut_low[u], cut_disc[cnode]);
          }
        }
      } else {
        while (!descended && *cut_stack[at].lit_at != SENTINEL_LIT) {
          const Lit l = *(cut_stack[at].lit_at++);
          if (is_unknown(l)) descended = visit(u, parent, l.var(), nullptr);
        }
      }
      if (descended) continue;

      // u is finished
      cut_stack.pop_back();
      if (cut_stack.empty()) break;
      const uint32_t p = cut_stack.back().node;
      cut_low[p] = std::min(cut_low[p], cut_low[u]);
      cut_subsz[p] += cut_subsz[u];
      if (p <= max_var) {
        if (p == root) root_children++;
        if (p == root || cut_low[u] >= cut_disc[p]) {
          cut_sep_sum[p] += cut_subsz[u];
          cut_sep_max[p] = std::max(cut_sep_max[p], cut_subsz[u]);
        }
      }
    }

    const uint32_t tree_sz = cut_tree_vars.size();
    tot_unknown += tree_sz;
    for(const auto& v: cut_tree_vars) {
      if (cut_sep_sum[v] == 0) continue;
      if (v == root && root_children < 2) continue;
      const uint32_t rest = tree_sz - 1 - cut_sep_sum[v]; // 0 for the root
      const uint32_t largest = std::max(rest, cut_sep_max[v]);
      cut_gain[v] = tree_sz - 1 - largest;
    }
  }
  return tot_unknown;
}

void CompAnalyzer::check_cut_gains(const Comp& comp) {
  vector<uint32_t> vars;
  all_vars_in_comp(comp, vt) if (is_unknown(*vt)) vars.push_back(*vt);
  vector<char> seen(max_var+1, 0);
  vector<uint32_t> todo;

  // flood fill from "from", never entering "removed" (0 = nothing removed)
  auto fill = [&](const uint32_t from, const uint32_t removed) -> uint32_t {
    uint32_t sz = 0;
    todo.clear();
    todo.push_back(from);
    seen[from] = 1;
    auto go = [&](const uint32_t w) {
      if (w != removed && is_unknown(w) && !seen[w]) { seen[w] = 1; todo.push_back(w); }
    };
    while (!todo.empty()) {
      const uint32_t u = todo.back();
      todo.pop_back();
      sz++;
      uint32_t const* bins = holder.begin_bin(u);
      for(uint32_t i = 0; i < holder.orig_size_bin(u); i++) go(bins[i]);
      ClData const* longs = holder.begin_long(u);
      for(uint32_t i = 0; i < holder.orig_size_long(u); i++) {
        const ClData& d = longs[i];
        if (d.id < max_tri_clid) {
          if (is_true(d.get_lit1()) || is_true(d.get_lit2())) continue;
          go(d.get_lit1().var());
          go(d.get_lit2().var());
        } else {
          Lit const* start = long_clauses_data.data()+d.off;
          bool sat = false;
          for (auto it_l = start; *it_l != SENTINEL_LIT; it_l++) if (is_true(*it_l)) sat = true;
          if (sat) continue;
          for (auto it_l = start; *it_l != SENTINEL_LIT; it_l++) go(it_l->var());
        }
      }
    }
    return sz;
  };

  for(const auto& rem: vars) {
    // the connected piece rem lives in
    for(const auto& v: vars) seen[v] = 0;
    const uint32_t tree_sz = fill(rem, 0);
    vector<uint32_t> tree;
    for(const auto& v: vars) if (seen[v]) tree.push_back(v);
    for(const auto& v: vars) seen[v] = 0;
    uint32_t largest = 0;
    uint32_t pieces = 0;
    for(const auto& v: tree) if (v != rem && !seen[v]) {
      largest = std::max(largest, fill(v, rem));
      pieces++;
    }
    const uint32_t expect = pieces >= 2 ? tree_sz - 1 - largest : 0;
    if (expect != cut_gain[rem]) {
      cout << "ERROR: cut gain of var " << rem << " is " << cut_gain[rem]
        << " but brute force says " << expect << " pieces: " << pieces
        << " tree_sz: " << tree_sz << " largest: " << largest << endl;
      release_assert(false);
    }
  }
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
