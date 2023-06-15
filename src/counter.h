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

#pragma once

#include "clauseallocator.h"
#include "common.h"
#include "cryptominisat5/cryptominisat.h"
#include "primitive_types.h"
#include "statistics.h"
#include "instance.h"
#include "comp_management.h"

#include "MersenneTwister.h"
#include "comp_management.h"
#include "boundedqueue.h"
#include "TreeDecomposition.h"
#include "structures.h"
#include <deque>
#include <map>

using std::pair;
using std::map;

using std::deque;

enum retStateT
{
  EXIT,
  RESOLVED,
  PROCESS_COMPONENT,
  BACKTRACK,
  GO_AGAIN
};

struct VS {
  VS() {}
  VS(uint32_t _v, double _score1, uint32_t _score2) : v(_v), score1(_score1), score2(_score2) {}
  bool operator<(const VS& other) const {
    if (score1 != other.score1) return score1 > other.score1;
    else return score2 > other.score2;
  }
  uint32_t v;
  double score1 = 0;
  uint32_t score2 = 0;
};

class ClauseAllocator;

// There is only one counter
class Counter : public Instance {
public:
  Counter(const CounterConfiguration& conf);
  ~Counter();
  friend class ClauseAllocator;
  struct ConflictData {
    ConflictData() :
      nHighestLevel(-1),
      bOnlyOneLitFromHighest(false)
    {}

    int32_t nHighestLevel;
    bool bOnlyOneLitFromHighest;
  };
  ConflictData find_conflict_level();

  double scoreOf(VariableIndex v) {
    if (config_.branch_type == branch_t::sharptd) {
      double score = 0;
      score += comp_manager_->scoreOf(v);
      score += 10*watches_[Lit(v, false)].activity + 10*watches_[Lit(v, true)].activity;
      score += tdscore[v];
      return score;
    } else {
      assert(config_.branch_type == branch_t::old_ganak);
      return
        comp_manager_->scoreOf(v)*act_inc +
        10*watches_[Lit(v, false)].activity + 10*watches_[Lit(v, true)].activity;
    }
  }
  mpz_class outer_count(CMSat::SATSolver* solver = NULL);
  CounterConfiguration &config() { return config_; }
  DataAndStatistics &statistics() { return stats; }
  void set_target_polar(const vector<CMSat::lbool>& model);
  void set_indep_support(const set<uint32_t>& indeps);
  void add_irred_cl(const vector<Lit>& lits);
  void add_red_cl(const vector<Lit>& lits, int lbd = -1);
  void get_activities(vector<double>& acts, vector<uint8_t>& polars, double& ret_act_inc, vector<uint32_t>& comp_acts) const;
  void set_activities(const vector<double>& act, const vector<uint8_t>& polars, double act_inc, vector<uint32_t>& comp_acts);
  const DataAndStatistics& get_stats() const;
  void shuffle_activities(MTRand& mtrand2);
  void end_irred_cls();
  void get_unit_cls(vector<Lit>& units) const;
  void init_activity_scores();
  void set_next_restart(uint64_t next) { config_.next_restart = next; }
  bqueue<uint32_t> comp_size_q;
  int32_t dec_level() const { return decision_stack_.get_decision_level(); }
  void print_restart_data() const;
  double get_start_time() const { return start_time;}
  void fill_cl(const Antecedent& ante, Lit*& c, uint32_t& size, Lit p) const;

  // deal with saved uip
  enum class SavedUIPRet {prop_again, ret_false, cont};
  SavedUIPRet deal_with_saved_uips();

  // test public
  uint64_t check_count(bool include_all_dec = false, int32_t single_var = -1);

private:
  mpz_class count(vector<Lit>& largest_cube_ret);
  CMSat::SATSolver* sat_solver = NULL;
  bool isindependent = true;
  vector<double> scores;
  bqueue<uint32_t> depth_q;
  bqueue<double, double> cache_miss_rate_q;
  vector<VS> vars_scores; // for branch picking

  // Temporaries, used during recordLastUIPClause
  mutable vector<Lit> tmpLit; //used in recoredLastUIPClause
  vector<uint32_t> toClear;
  set<Lit> toSet;

  // Used during minimizeAndStoreUIPClause
  deque<Lit> tmp_clause_minim;

  // Temporaries for failedLitProbeInternal
  vector<Lit> test_lits;
  vector<uint8_t> viewed_lits;

  double start_time;
  MTRand mtrand;

  DecisionStack decision_stack_;
  vector<Lit> trail;
  uint32_t qhead = 0;
  ComponentManager* comp_manager_ = NULL;

  // the last time conflict clauses have been deleted
  uint64_t last_reduceDB_conflicts = 0;
  // the last time the conflict clause storage has been compacted
  uint64_t last_ccl_cleanup_decs_ = 0;

  void simplePreProcess();
  bool prepFailedLiteralTest();
  bool is_implied(const vector<Lit>& cp);
  void check_implied(const vector<Lit>& cl);

  SOLVER_StateT countSAT();
  bool decideLiteral();
  uint32_t find_best_branch_gpmc(bool do_indep);
  uint32_t find_best_branch(bool do_indep);
  bool clause_falsified(const vector<Lit>& cl) const;
  bool clause_asserting(const vector<Lit>& cl) const;
  template<class T> bool clause_satisfied(const T& cl) const;
  bool prop_and_probe();
  bool failed_lit_probe();
  bool failed_lit_probe_no_bprop();
  bool failed_lit_probe_with_bprop();
  bool one_lit_probe(Lit lit, bool set);
  void computeLargestCube();
  void compute_score(TreeDecomposition& tdec);
  void td_decompose();

  // this is the actual BCP algorithm
  // starts propagating all literal in trail_
  // beginning at offset start_at_trail_ofs
  bool propagate();
  template<uint32_t start = 2>
  void get_maxlev_maxind(ClauseOfs ofs, int32_t& maxlev, uint32_t& maxind);
  bool check_watchlists() const;

  void print_all_levels();
  bool restart_if_needed();
  retStateT backtrack_nonindep();
  retStateT backtrack();
  void print_dec_info() const;
  template<class T> void print_cl(const T& cl) const;
  template<class T> void v_print_cl(const T& cl) const;
  void print_conflict_info() const;
  void print_comp_stack_info() const;

  // if on the current decision level
  // a second branch can be visited, RESOLVED is returned
  // otherwise returns BACKTRACK
  retStateT resolveConflict();
  void go_back_to(int32_t backj);
  uint32_t find_lev_to_set(int32_t other_lev);
  size_t find_backtrack_level_of_learnt();
  void print_trail(bool check_entail = true, bool check_anything = true) const;
  void check_trail(bool check_entail = true) const;

  void setLiteral(const Lit lit, int32_t dec_lev, Antecedent ant = Antecedent())
  {
    assert(val(lit) == X_TRI);
    if (ant.isNull())
      print_debug("setLiteral called with a decision. Lit: " << lit << " lev: " << dec_lev);
    else print_debug("-> lit propagated: " << lit << " trail pos will be: " << trail.size());

    VERBOSE_DEBUG_DO(cout << "setting lit: " << lit << " to lev: " << dec_lev << " cur val: " << lit_val_str(lit) << " ante: " << ant << " sublev: " << trail.size() << endl);
    var(lit).decision_level = dec_lev;
    var(lit).ante = ant;
    if (!ant.isNull()) {
      var(lit).last_polarity = lit.sign();
      var(lit).set_once = true;
    }
    var(lit).sublevel = trail.size();
    qhead = std::min<uint32_t>(qhead, trail.size());
    trail.push_back(lit);
    __builtin_prefetch(watches_[lit.neg()].binary_links_.data());
    __builtin_prefetch(watches_[lit.neg()].watch_list_.data());
    if (ant.isAnt() && ant.isAClause()) {
      Clause& cl = *alloc->ptr(ant.asCl());
      if (cl.red && cl.lbd > lbd_cutoff) {
        cl.increaseScore();
        cl.update_lbd(calc_lbd(cl));
      }
    }
    lit_values_[lit] = T_TRI;
    lit_values_[lit.neg()] = F_TRI;
  }

  void checkProbabilisticHashSanity() const {
      const uint64_t t = stats.num_cache_look_ups_ + 1;
      // The +32 is because there is actually another hash, which is 32b and is used
      // by the caching subsystem. Both must match.
      if (2 * log2(t) > log2(config_.delta) + (64+32) * 0.9843) {
        // 1 - log_2(2.004)/64 = 0.9843
        cout << "ERROR: We need to change the hash range (-1)" << endl;
        exit(-1);
      }
  }

  void setConflictState(Lit litA, Lit litB)
  {
    conflLit = litA;
    confl = Antecedent(litB);
  }

  void setConflictState(Clause* cl)
  {
    if (cl->red && cl->lbd > lbd_cutoff) {
      cl->increaseScore();
      cl->update_lbd(calc_lbd(*cl));
    }
    confl = Antecedent(alloc->get_offset(cl));
    conflLit = NOT_A_LIT;
  }

  // The literals that have been set in this decision level
  vector<Lit>::const_iterator top_declevel_trail_begin() const
  {
    return trail.begin() + variables_[decision_stack_.top().var].sublevel;
  }
  vector<Lit>::iterator top_declevel_trail_begin()
  {
    return trail.begin() + variables_[decision_stack_.top().var].sublevel;
  }

  void init_decision_stack()
  {
    decision_stack_.clear();
    trail.clear();
    // initialize the stack to contain at least level zero
    decision_stack_.push_back(StackLevel(
          1, // super comp
          2)); //comp stack offset
    decision_stack_.top().on_path_to_target_ = true;

    // I guess this is needed so the system later knows it's fully counted
    // since this is only a dummy.
    decision_stack_.back().change_to_right_branch();
  }

  const Lit &top_dec_lit() const
  {
    return *top_declevel_trail_begin();
  }

  uint32_t trail_at_dl(uint32_t dl) const {
    return variables_[decision_stack_.at(dl).var].sublevel;
  }

  uint32_t trail_at_top() const {
    return variables_[decision_stack_.top().var].sublevel;
  }

  void reactivate_comps_and_backtrack_trail([[maybe_unused]] bool check_ws = true)
  {
    VERBOSE_PRINT("->reactivate and backtrack...");
    auto jt = top_declevel_trail_begin();
    auto it = jt;
    /* qhead = std::min<int32_t>(qhead, jt - trail.begin()); // TODO something is wrong here. */
    for (; it != trail.end(); it++) {
      int32_t dl = var(*it).decision_level;
      assert(dl != -1);
      if (dl < decision_stack_.get_decision_level()) {
        var(*it).sublevel = jt - trail.begin();
        *jt++ = *it;
        VERBOSE_DEBUG_DO(cout << "Backing up, setting sublevel lit " << std::setw(5) << *it << " sublev: " << var(*it).sublevel << endl);
      } else {
        VERBOSE_DEBUG_DO(cout << "Backing up, unsetting: " << *it << " lev: " << var(*it).decision_level << " ante was: " << var(*it).ante << endl);
        unSet(*it);
      }
    }
    VERY_SLOW_DEBUG_DO(if (check_ws && !check_watchlists()) {print_trail(false, false);assert(false);});
    comp_manager_->cleanRemainingComponentsOf(decision_stack_.top());
    trail.resize(jt - trail.begin());
    decision_stack_.top().resetRemainingComps();
  }

  /////////////////////////////////////////////
  //  Conflict analysis below
  /////////////////////////////////////////////

  // if the state name is CONFLICT,
  // then violated_clause contains the clause determining the conflict;
  Lit conflLit = NOT_A_LIT;
  Antecedent confl;
  // this is an array of all the clauses found
  // during the most recent conflict analysis
  // it might contain more than 2 clauses
  // but always will have:
  //      uip_clause the 1UIP clause found
  //  possible clauses in between will be other UIP clauses
  vector<Lit> uip_clause;
  vector<vector<Lit>> saved_uip_cls;

  int32_t get_confl_maxlev(const Lit p) const;
  void create_fake(Lit p, uint32_t& size, Lit*& c) const;
  void recordLastUIPCauses();
  void minimizeUIPClause();
  uint32_t abstractLevel(const uint32_t x) const;
  bool litRedundant(Lit p, uint32_t abstract_levels);
  vector<Lit> analyze_stack;
  void recursiveConfClauseMin();
  bool takeSolution();
  bool get_polarity(const uint32_t var) const;
  bool standard_polarity(const uint32_t var) const;

  // Vivification
  int64_t v_tout;
  vector<Lit> v_tmp;
  vector<Lit> v_tmp2;
  vector<Lit> v_cl;
  uint64_t last_confl_vivif = 0;
  map<ClauseOfs, pair<Lit, Lit>> ws_pos;
  void v_cl_repair(ClauseOfs off);
  void vivify_cls(vector<ClauseOfs>& cls);
  void vivify_clauses();
  bool vivify_cl(const ClauseOfs off);
  bool v_propagate();
  void v_backtrack();
  void v_unset(const Lit l);
  void v_enqueue(const Lit l);
  TriValue v_val(const Lit l) const;
  void v_fix_watch(Clause& cl, uint32_t i);
  template<class T> bool propagating_cl(T& cl) const;
  template<class T> bool conflicting_cl(T& cl) const;
  template<class T> bool should_have_propagated_earlier(const T& cl) const;
  void v_new_lev();
  template<class T> bool v_clause_satisfied(const T& cl) const;
  void vivif_backtrack();
  vector<Lit> v_trail;
  uint32_t v_qhead;
  uint32_t v_lev;
  vector<int32_t> v_levs;
  uint32_t v_backtrack_to;
  LiteralIndexedVector<TriValue> v_values;

  void print_stat_line();
  uint64_t next_print_stat_cache = 20000;
  uint64_t next_print_stat_confl = 5000;

  // indicates if we have called end_irred_cls()
  bool ended_irred_cls = false;
};

template<class T> void Counter::print_cl(const T& cl) const {
  for(uint32_t i = 0; i < cl.size(); i ++) {
    const auto l = cl[i];
    cout << "lit " << std::setw(6) << l
      << " lev: " << std::setw(4) << var(l).decision_level
      << " ante: " << std::setw(5) << std::left << var(l).ante
    << " val: " << lit_val_str(l) << endl;
  }
}

template<class T> void Counter::v_print_cl(const T& cl) const {
  for(uint32_t i = 0; i < cl.size(); i ++) {
    const auto l = cl[i];
    cout << "lit " << std::setw(6) << l
      << " lev: " << std::setw(4) << v_levs[l.var()]
    << " val: " << val_str(v_val(l)) << endl;
  }
}

template<class T> bool Counter::conflicting_cl(T& cl) const {
  for(const auto&l: cl) {
    if (val(l) == T_TRI || val(l) == X_TRI) return false;
  }
  return true;
}

template<class T> bool Counter::propagating_cl(T& cl) const {
  uint32_t unk = 0;
  for(const auto&l: cl) {
    if (val(l) == T_TRI) return false;
    if (val(l) == X_TRI) {unk++; if (unk>1) break;}
  }
  return unk == 1;
}
