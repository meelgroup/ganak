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

  // --synthesis (wDNNF): set up sharing of pure output vars. We build signed
  // occurrence data below so residual polarity purity can be computed per round.
  share_mode = (conf.synthesis && !conf.compile_fname.empty());
  // Shareability cutoff is opt_indep_support_end, NOT indep_support_end. Reason:
  // main-DPLL's find_best_branch picks vars `< opt_indep_support_end`, so vars in
  // [indep_support_end, opt_indep_support_end) ARE branched in the main search.
  // If those vars were shareable, synth_forced_lit would override the polarity
  // heuristic on every such decision (forcing orig_polarity, which defaults to
  // +v) -- a serious search-quality regression that empirically dominated the
  // wDNNF runtime cost (see git history). Vars >= opt_indep_support_end are
  // ONLY set by the SAT oracle (witness recording) or by propagation, so pinning
  // them inside the SAT loop is free of main-DPLL heuristic damage. Sharing
  // soundness is unchanged -- the SAT-loop pin (counter.cpp synth_forced_lit
  // call from the SAT branch) still enforces sibling agreement for the vars
  // that remain shareable. The vars we're DROPPING from the shareable set are
  // now treated exactly as in faithful mode: ordinary branch vars.
  indep_end = counter->get_opt_indep_support_end();

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
  // --synthesis: full signed literals per clause id (ids run 1..#long_cls).
  if (share_mode) { cls_lits.clear(); cls_lits.resize(long_irred_cls.size() + 2); }
  for (const auto& off: long_irred_cls) {
    const Clause& cl = *alloc->ptr(off);
    assert(cl.size() > 2);
    if (share_mode) cls_lits[max_clid].assign(cl.begin(), cl.end());
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

  // --synthesis (wDNNF): build signed binary occ lists and per-round scratch.
  if (share_mode) {
    claimed_share.assign(max_var + 1, 0);
    shareable.assign(max_var + 1, 0);
    seen_pos.assign(max_var + 1, 0);
    seen_neg.assign(max_var + 1, 0);
    bin_pos.assign(n, {});
    bin_neg.assign(n, {});
    for (uint32_t v = 1; v < n; v++) {
      for (const auto& bincl: watches[Lit(v, true)].binaries)
        if (bincl.irred()) bin_pos[v].push_back(bincl.lit());
      for (const auto& bincl: watches[Lit(v, false)].binaries)
        if (bincl.irred()) bin_neg[v].push_back(bincl.lit());
    }
    // ORIGINAL polarity: scan every irredundant clause once and record which
    // polarities of v occur. State-independent: never recomputed. Synthesis uses
    // this to pin a deterministic polarity when the residual is ambiguous (esp.
    // FREE = no active occurrence, which would otherwise default-pin differently
    // in different sibling sub-comps -- the wDNNF bug).
    orig_polarity.assign(n, 0);
    for (uint32_t v = 1; v < n; v++) {
      if (!bin_pos[v].empty()) orig_polarity[v] |= 1;
      if (!bin_neg[v].empty()) orig_polarity[v] |= 2;
    }
    for (uint32_t cid = 1; cid < cls_lits.size(); cid++) {
      for (const Lit l : cls_lits[cid]) {
        orig_polarity[l.var()] |= (l.sign() ? 1 : 2);
      }
    }
    // If NO output var has orig_polarity < 3 (i.e. every output var appears in
    // both polarities in the original CNF), compute_shareable_vars can NEVER
    // mark anything as shareable -- the orig_polarity==3 filter (the soundness
    // gate for "no globally-consistent pin polarity") rules them all out
    // regardless of the residual. In that case the whole share_mode pipeline is
    // dead weight: O(super-comp) work per analysis round and a SAT-loop pin
    // check per decision, all producing nothing. Disable share_mode up front.
    bool any_candidate = false;
    for (uint32_t v = indep_end; v <= max_var; v++) {
      if (orig_polarity[v] != 3) { any_candidate = true; break; }
    }
    if (!any_candidate) {
      share_mode = false;
      verb_print(1, "[compile-synthesis] no output var is originally pure"
        " -- share_mode disabled (faithful-mode performance)");
    } else {
      verb_print(1, "[compile-synthesis] wDNNF share-and-branch over pure outputs >= " << indep_end);
    }
  }
}

// SLOW_DEBUG: current residual polarity of v (see header). Recomputed fresh from
// `values`, so it reflects the live assignment regardless of shareable[] staleness.
int CompAnalyzer::residual_polarity(uint32_t v) {
  int res = 0;
  for (const Lit o : bin_pos[v]) if (!is_true(o)) { res |= 1; break; }
  for (const Lit o : bin_neg[v]) if (!is_true(o)) { res |= 2; break; }
  const ClData* longs = holder.begin_long(v);
  const uint32_t nl = holder.orig_size_long(v);
  for (uint32_t i = 0; i < nl; i++) {
    const auto& lits = cls_lits[longs[i].id];
    bool sat = false;
    for (const Lit l : lits) if (is_true(l)) { sat = true; break; }
    if (sat) continue;
    for (const Lit l : lits) if (l.var() == v) res |= (l.sign() ? 1 : 2);
  }
  return res;
}

// --synthesis (wDNNF): recompute shareable[] for the current super-comp.
// A var is initially shareable iff it is an unknown OUTPUT var (>= indep_end) that
// is PURE (appears in a single polarity) in the residual formula -- this is the
// weak-decomposability condition (Akshay et al. 2018): a shared var in a single
// polarity never causes a witness-extraction conflict. Then a demotion fixpoint
// un-shares one var of every *active* clause whose remaining literals are ALL
// shareable, so each clause keeps a bridging (non-shareable) var and is claimed by
// exactly one component -- nothing is dropped, keeping the circuit faithful as a
// Boolean function. (The model count stays meaningless; that's fine -- in
// --synthesis we synthesize, not count.) Inputs (< indep_end) are never shareable,
// so components stay disjoint over the inputs.
//
// SOUNDNESS NOTE (the historical ~1-2% synth-failure bug): purity is a *residual*
// property -- it depends on which clauses are still active. Two sibling sub-comps
// (post-decomposition) see DIFFERENT residuals once decisions inside one have
// flipped clauses on/off, so the same var can end up "pure positive" in one
// sibling and "pure negative" (or free) in another. synth_forced_lit then pins
// the var to opposite polarities across siblings, and the AND-of-siblings node
// loses wDNNF. The mitigation lives in synth_forced_lit and pin_polarity[].
void CompAnalyzer::compute_shareable_vars(const Comp& super_comp) {
  // Tier-2A memoization: if we ran this same analysis for the same super-comp
  // at the same global trail tstamp, shareable[] is still valid -- no decision
  // happened (tstamp unchanged) AND no other super-comp's analysis ran in
  // between (otherwise last_share_super_comp would point elsewhere). Skipping
  // saves an O(super-comp) walk and the demotion fixpoint.
  const uint64_t cur_tstamp = counter->get_tstamp();
  if (&super_comp == last_share_super_comp && cur_tstamp == last_share_tstamp) {
    stats.synth_shareable_memo_hits++;
    return;
  }
  last_share_super_comp = &super_comp;
  last_share_tstamp = cur_tstamp;
  stats.synth_shareable_calls++;
  // --- residual polarity purity ---
  all_vars_in_comp(super_comp, vt) {
    const uint32_t v = *vt;
    shareable[v] = 0; seen_pos[v] = 0; seen_neg[v] = 0;
  }
  // long/ternary clauses (signed lits via cls_lits); skip satisfied clauses.
  all_cls_in_comp(super_comp, ci) {
    const auto& lits = cls_lits[*ci];
    bool sat = false;
    for (const Lit l : lits) if (is_true(l)) { sat = true; break; }
    if (sat) continue;
    for (const Lit l : lits) {
      const uint32_t vv = l.var();
      if (vv < indep_end || !is_unknown(vv)) continue;
      if (l.sign()) seen_pos[vv] = 1; else seen_neg[vv] = 1;
    }
  }
  // binary clauses (signed via bin_pos/bin_neg); (v OR o) is active iff o not true.
  all_vars_in_comp(super_comp, vt) {
    const uint32_t v = *vt;
    if (v < indep_end || !is_unknown(v)) continue;
    for (const Lit o : bin_pos[v]) if (!is_true(o)) { seen_pos[v] = 1; break; }
    for (const Lit o : bin_neg[v]) if (!is_true(o)) { seen_neg[v] = 1; break; }
  }
  // pure output var -> shareable, but only if it's also ORIGINALLY pure (or
  // never occurs) in the formula. Originally-impure outputs (both +v and -v
  // appear somewhere in the CNF) are NEVER shared because two sibling sub-comps
  // can land in residuals where v is pure with DIFFERENT polarities (one sibling
  // has the +v clauses satisfied -> pure neg; another has -v satisfied -> pure
  // pos). No globally-consistent pin polarity exists for such a var, so sharing
  // it inevitably violates wDNNF.
  all_vars_in_comp(super_comp, vt) {
    const uint32_t v = *vt;
    if (v < indep_end || !is_unknown(v)) continue;
    if (seen_pos[v] && seen_neg[v]) continue;     // impure right now
    if (orig_polarity[v] == 3) continue;           // both polarities in original CNF
    shareable[v] = 1;
    stats.synth_shareable_marked++;
  }

  // --- demotion fixpoint: no active clause may be entirely shareable ---
  all_cls_in_comp(super_comp, ci) {
    const auto& lits = cls_lits[*ci];
    bool sat = false;
    for (const Lit l : lits) if (is_true(l)) { sat = true; break; }
    if (sat) continue;
    int demote = -1; bool all_share = true;
    for (const Lit l : lits) {
      if (is_false(l)) continue;
      const uint32_t vv = l.var();
      if (shareable[vv]) demote = (int)vv;
      else { all_share = false; break; }
    }
    if (all_share && demote >= 0) { shareable[demote] = 0; stats.synth_shareable_demoted++; }
  }
  // binaries: handle each (v OR o) once, demoting v at its smaller-index endpoint.
  all_vars_in_comp(super_comp, vt) {
    const uint32_t v = *vt;
    if (!shareable[v]) continue;
    for (const Lit o : bin_pos[v])
      if (!is_true(o) && !is_false(o) && o.var() > v && shareable[o.var()]) {
        shareable[v] = 0; stats.synth_shareable_demoted++; break;
      }
    if (!shareable[v]) continue;
    for (const Lit o : bin_neg[v])
      if (!is_true(o) && !is_false(o) && o.var() > v && shareable[o.var()]) {
        shareable[v] = 0; stats.synth_shareable_demoted++; break;
      }
  }
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
      if (counter->compiling()) counter->compile_add_free_var(v);
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
    // --synthesis (wDNNF): a shareable var (pure output, see compute_shareable_vars)
    // is a component member but does not bridge -- its other clauses aren't pulled
    // in, keeping comps disjoint over inputs. Re-claimable by later siblings (see
    // make_comp_from_archetype). The demotion fixpoint guarantees every active
    // clause still has a non-shareable var, so no clause is ever dropped.
    if (is_shareable(v)) { claimed_share[v] = 1; continue; }
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
    // SLOW_DEBUG sanity bound: outside synthesis share mode every bin is counted
    // symmetrically from both endpoints in every level, so a sub-comp's bin
    // count never exceeds its super-comp's. Under synthesis share mode the
    // count is ASYMMETRIC -- a shareable endpoint skips its bin traversal, so a
    // bin (u,v) counts once if exactly one endpoint is shareable. When v's
    // shareable status changes between the parent level (counted once) and this
    // level (counted twice if both endpoints are now non-shareable), the
    // inequality breaks. The bin SET is still a subset; this is a counting
    // artifact, not a real over-claim of bins.
    SLOW_DEBUG_DO(
      if (!share_mode) assert(archetype.num_bin_cls <= sup_comp_bin_cls);
    );

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
