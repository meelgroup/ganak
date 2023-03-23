/*
 * solver.h
 *
 *  Created on: Aug 23, 2012
 *      Author: marc
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "common.h"
#include "primitive_types.h"
#include "statistics.h"
#include "instance.h"
#include "comp_management.h"
#include "solver_config.h"
#include "MersenneTwister.h"

#include <sys/time.h>
#include <deque>

using std::deque;

enum retStateT
{
  EXIT,
  RESOLVED,
  PROCESS_COMPONENT,
  BACKTRACK,
  RESTART
};

// There is only one solver
class Solver : public Instance
{
public:
  Solver() { mtrand.seed((uint32_t)0U);}
  void solve(const std::string &file_name);
  SolverConfiguration &config() { return config_; }
  DataAndStatistics &statistics() { return stats; }

private:
  // Temporaries, used during recordLastUIPClause
  vector<uint8_t> tmp_seen;
  vector<Lit> tmp_clause; //used in recoredLastUIPClause
  vector<uint32_t> toClear;

  // Used during minimizeAndStoreUIPClause
  deque<Lit> tmp_clause_minim;

  // Temporaries for failedLitProbeInternal
  vector<Lit> test_lits;
  LiteralIndexedVector<uint8_t> viewed_lits;

  double time_start;
  SolverConfiguration config_;
  MTRand mtrand;

  DecisionStack decision_stack_;
  vector<Lit> trail;
  ComponentManager comp_manager_ = ComponentManager(
          config_,stats, lit_values_, indep_support_);

  // the last time conflict clauses have been deleted
  uint64_t last_ccl_deletion_decs_ = 0;
  // the last time the conflict clause storage has been compacted
  uint64_t last_ccl_cleanup_decs_ = 0;

  bool simplePreProcess();
  bool prepFailedLiteralTest();
  // we assert that the formula is consistent
  // and has not been found UNSAT yet
  // hard wires all assertions in the literal stack into the formula
  // removes all set variables and essentially reinitiallizes all
  // further data
  void HardWireAndCompact();

  SOLVER_StateT countSAT();
  void decideLiteral();
  bool failedLitProbe();
  bool failedLitProbeInternal();
  void computeLargestCube();

  void decayActivitiesOf(Component &comp)
  {
    for (auto it = comp.varsBegin(); *it != varsSENTINEL; it++)
    {
      litWatchList(Lit(*it, true)).activity_score_ *= 0.5;
      litWatchList(Lit(*it, false)).activity_score_ *= 0.5;
    }
  }

  // this is the actual BCP algorithm
  // starts propagating all literal in lit_stack_
  // beginning at offset start_at_stack_ofs
  bool propagate(const uint32_t start_at_stack_ofs);

  retStateT backtrack();

  // if on the current decision level
  // a second branch can be visited, RESOLVED is returned
  // otherwise returns BACKTRACK
  retStateT resolveConflict();

  float scoreOf(VariableIndex v)
  {
    float score = 1.0; //comp_manager_.scoreOf(v);
    score += litWatchList(Lit(v, true)).activity_score_;
    score += litWatchList(Lit(v, false)).activity_score_;

    return score;
  }

  bool setLiteralIfFree(const Lit lit,
                        const Antecedent ant = Antecedent(NOT_A_CLAUSE))
  {
    if (lit_values_[lit] != X_TRI) return false;
    if (ant == Antecedent(NOT_A_CLAUSE)) print_debug("setLiteralIfFree called with NOT_A_CLAUSE as antecedent (i.e. it's a decision). Lit: " << lit);
//    else print_debug("Literal propagated: " << lit);

    var(lit).decision_level = decision_stack_.get_decision_level();
    var(lit).ante = ant;
    var(lit).polarity = lit.sign();
    var(lit).set = true;
    trail.push_back(lit);
    if (ant.isAClause() && ant.asCl() != NOT_A_CLAUSE)
      getHeaderOf(ant.asCl()).increaseScore();
    lit_values_[lit] = T_TRI;
    lit_values_[lit.neg()] = F_TRI;
    return true;
  }

  void printOnlineStats();
  void checkProbabilisticHashSanity() const {
      const uint32_t t = stats.num_cache_look_ups_ + 1;
      if (2 * log2(t) > log2(config_.delta) + 64 * config_.hashrange * 0.9843) {
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

  vector<Lit>::const_iterator TOSLiteralsBegin()
  {
    return trail.begin() + decision_stack_.top().lit_stack_ofs();
  }

  void initStack()
  {
    decision_stack_.clear();
    trail.clear();
    // initialize the stack to contain at least level zero
    decision_stack_.push_back(StackLevel(
          1, // super comp
          0, // trail offset
          2)); //comp stack offset
    decision_stack_.back().changeBranch();
  }

  const Lit &TOS_decLit()
  {
    assert(decision_stack_.top().lit_stack_ofs() < trail.size());
    return trail[decision_stack_.top().lit_stack_ofs()];
  }

  void reactivateTOS()
  {
    for (auto it = TOSLiteralsBegin(); it != trail.end(); it++) unSet(*it);
    comp_manager_.cleanRemainingComponentsOf(decision_stack_.top());
    trail.resize(decision_stack_.top().lit_stack_ofs());
    decision_stack_.top().resetRemainingComps();
  }

  bool fail_test(const Lit lit)
  {
    uint32_t sz = trail.size();
    // we increase the decLev artificially
    // s.t. after the tentative BCP call, we can learn a conflict clause
    // relative to the assignment of *jt
    decision_stack_.startFailedLitTest();
    print_debug("Fail testing lit: " << lit);
    setLiteralIfFree(lit);

    assert(!hasAntecedent(lit));

    bool bSucceeded = propagate(sz);
    if (!bSucceeded) recordAllUIPCauses();

    decision_stack_.stopFailedLitTest();

    while (trail.size() > sz)
    {
      unSet(trail.back());
      trail.pop_back();
    }
    return bSucceeded;
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
  void minimizeAndStoreUIPClause(Lit uipLit,
                                 vector<Lit> &tmp_clause,
                                 const vector<uint8_t>& seen);
  void storeUIPClause(Lit uipLit, vector<Lit> &tmp_clause);
  int getAssertionLevel() const { return assertion_level_; }
  bool takeSolution();
  bool get_polarity(const uint32_t var);
};

#endif /* SOLVER_H_ */
