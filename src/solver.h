/*
 * solver.h
 *
 *  Created on: Aug 23, 2012
 *      Author: marc
 */

#pragma once

#include "common.h"
#include "primitive_types.h"
#include "statistics.h"
#include "instance.h"
#include "comp_management.h"
#include "solver_config.h"

#include "MersenneTwister.h"
#include "comp_management.h"
#include "boundedqueue.h"
#include <deque>

using std::deque;

enum retStateT
{
  EXIT,
  RESOLVED,
  PROCESS_COMPONENT,
  BACKTRACK,
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

// There is only one counter
class Counter : public Instance
{
public:
  Counter(const CounterConfiguration& conf);
  ~Counter();

  double scoreOf(VariableIndex v) {
    return
      comp_manager_->scoreOf(v)*act_inc +
      10*watches_[Lit(v, false)].activity + 10*watches_[Lit(v, true)].activity;
  }
  mpz_class count(vector<Lit>& largest_cube_ret);
  CounterConfiguration &config() { return config_; }
  DataAndStatistics &statistics() { return stats; }
  void set_target_polar(const vector<CMSat::lbool>& model);
  void set_indep_support(const set<uint32_t>& indeps);
  void add_red_cl(const vector<Lit>& lits);
  void get_activities(vector<double>& acts, vector<uint8_t>& polars, double& ret_act_inc, vector<uint32_t>& comp_acts) const;
  void set_activities(const vector<double>& act, const vector<uint8_t>& polars, double act_inc, vector<uint32_t>& comp_acts);
  const DataAndStatistics& get_stats() const;
  void shuffle_activities(MTRand& mtrand2);
  void end_irred_cls();
  size_t get_num_irred_long_cls() const { return conflict_clauses_.size(); }
  void get_unit_cls(vector<Lit>& units) const;
  void init_activity_scores();
  void set_next_restart(uint64_t next) { config_.next_restart = next; }
  bqueue<uint32_t> comp_size_queue;
  uint64_t dec_level() const { return decision_stack_.size(); }

private:
  vector<double> scores;
  bqueue<uint32_t> depth_queue;
  bqueue<double, double> cache_miss_rate_queue;
  vector<VS> vars_scores; // for branch picking

  // Temporaries, used during recordLastUIPClause
  vector<Lit> tmp_clause; //used in recoredLastUIPClause
  vector<uint32_t> toClear;
  set<Lit> toSet;

  // Used during minimizeAndStoreUIPClause
  deque<Lit> tmp_clause_minim;

  // Temporaries for failedLitProbeInternal
  vector<Lit> test_lits;
  vector<uint8_t> viewed_lits;

  double time_start;
  CounterConfiguration config_;
  MTRand mtrand;

  DecisionStack decision_stack_;
  vector<Lit> trail;
  ComponentManager* comp_manager_ = NULL;

  // the last time conflict clauses have been deleted
  uint64_t last_ccl_deletion_decs_ = 0;
  // the last time the conflict clause storage has been compacted
  uint64_t last_ccl_cleanup_decs_ = 0;

  void simplePreProcess();
  bool prepFailedLiteralTest();

  SOLVER_StateT countSAT();
  void decideLiteral(Lit lit = NOT_A_LIT);
  uint32_t find_best_branch();
  double alternate_score(uint32_t v, bool value);
  bool prop_and_probe();
  bool failed_lit_probe();
  bool failed_lit_probe_no_bprop();
  bool failed_lit_probe_with_bprop();
  bool one_lit_probe(Lit lit, bool set);
  void computeLargestCube();

  // this is the actual BCP algorithm
  // starts propagating all literal in trail_
  // beginning at offset start_at_trail_ofs
  bool propagate(const uint32_t start_at_trail_ofs);

  void print_all_levels();
  bool restart_if_needed();
  retStateT backtrack();

  // if on the current decision level
  // a second branch can be visited, RESOLVED is returned
  // otherwise returns BACKTRACK
  retStateT resolveConflict();

  bool setLiteralIfFree(const Lit lit, const Antecedent ant = Antecedent(NOT_A_CLAUSE),
      const bool bothprop = false)
  {
    if (lit_values_[lit] != X_TRI) return false;
    if (ant == Antecedent(NOT_A_CLAUSE)) print_debug("setLiteralIfFree called with NOT_A_CLAUSE as antecedent (i.e. it's a decision). Lit: " << lit);
    else print_debug("-> lit propagated: " << lit);

    var(lit).decision_level = decision_stack_.get_decision_level();
    var(lit).ante = ant;
    if (bothprop) var(lit).bprop = true;
    if (ant != Antecedent(NOT_A_CLAUSE)) {
      var(lit).last_polarity = lit.sign();
      var(lit).set_once = true;
    }
    trail.push_back(lit);
    __builtin_prefetch(watches_[lit.neg()].binary_links_.data());
    __builtin_prefetch(watches_[lit.neg()].watch_list_.data());
    if (ant.isAClause() && ant.asCl() != NOT_A_CLAUSE)
      getHeaderOf(ant.asCl()).increaseScore();
    lit_values_[lit] = T_TRI;
    lit_values_[lit.neg()] = F_TRI;
    return true;
  }

  void printOnlineStats();
  void checkProbabilisticHashSanity() const {
      const uint64_t t = stats.num_cache_look_ups_ + 1;
      if (2 * log2(t) > log2(config_.delta) + 64 * 0.9843) {
        // 1 - log_2(2.004)/64 = 0.9843
        cout << "ERROR: We need to change the hash range (-1)" << endl;
        exit(-1);
      }
  }

  void setConflictState(Lit litA, Lit litB)
  {
    violated_clause.clear();
    violated_clause.push_back(litA);
    violated_clause.push_back(litB);
  }
  void setConflictState(ClauseOfs cl_ofs)
  {
    getHeaderOf(cl_ofs).increaseScore();
    violated_clause.clear();
    for (auto it = beginOf(cl_ofs); *it != SENTINEL_LIT; it++)
      violated_clause.push_back(*it);
  }

  // The literals that have been set in this decision level
  vector<Lit>::const_iterator top_declevel_trail_begin() const
  {
    return trail.begin() + decision_stack_.top().trail_ofs();
  }
  vector<Lit>::iterator top_declevel_trail_begin()
  {
    return trail.begin() + decision_stack_.top().trail_ofs();
  }

  void init_decision_stack()
  {
    decision_stack_.clear();
    trail.clear();
    // initialize the stack to contain at least level zero
    decision_stack_.push_back(StackLevel(
          1, // super comp
          0, // trail offset
          2)); //comp stack offset
    decision_stack_.top().on_path_to_target_ = true;

    // I guess this is needed so the system later knows it's fully counted
    // since this is only a dummy.
    decision_stack_.back().change_to_right_branch();
  }

  const Lit &top_dec_lit()
  {
    assert(decision_stack_.top().trail_ofs() < trail.size());
    return *top_declevel_trail_begin();
  }

  void reactivate_comps_and_backtrack_trail(bool print = false)
  {
    for (auto it = top_declevel_trail_begin(); it != trail.end(); it++) unSet(*it);
    comp_manager_->cleanRemainingComponentsOf(decision_stack_.top());
    trail.resize(decision_stack_.top().trail_ofs());
    if (print) cout << "Forgetting decision: "
      << std::setw(1) << (decision_stack_.top().is_right_branch() ? "-" : "")
      << std::setw(6) << decision_stack_.top().getbranchvar()
        << " count: " << decision_stack_.top().getTotalModelCount() << endl;
    decision_stack_.top().resetRemainingComps();
  }

  /////////////////////////////////////////////
  //  Conflict analysis below
  /////////////////////////////////////////////

  // if the state name is CONFLICT,
  // then violated_clause contains the clause determining the conflict;
  vector<Lit> violated_clause;
  // this is an array of all the clauses found
  // during the most recent conflict analysis
  // it might contain more than 2 clauses
  // but always will have:
  //      uip_clauses_.front() the 1UIP clause found
  //      uip_clauses_.back() the lastUIP clause found
  //  possible clauses in between will be other UIP clauses
  vector<vector<Lit>> uip_clauses_;

  // the assertion level of uip_clauses_.back()
  // or (if the decision variable did not have an antecedent
  // before) then assertionLevel_ == DL;
  int assertion_level_ = 0;

  // build conflict clauses from most recent conflict
  // as stored in state_.violated_clause
  // solver state must be CONFLICT to work;
  // this first method record only the last UIP clause
  // so as to create clause that asserts the current decision
  // literal
  void recordLastUIPCauses();
  void recordAllUIPCauses();
  void minimizeAndStoreUIPClause(Lit uipLit, vector<Lit> &tmp_clause);
  void storeUIPClause(Lit uipLit, vector<Lit> &tmp_clause);
  int getAssertionLevel() const { return assertion_level_; }
  bool takeSolution();
  bool get_polarity(const uint32_t var);

  void print_stat_line();
  uint64_t next_print_stat = 20000;

  // indicates if we have called end_irred_cls()
  bool ended_irred_cls = false;
};
