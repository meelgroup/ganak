/*
 * solver.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: marc
 */
#include "solver.h"

#include <deque>
#include <algorithm>
#include <ios>
#include "common.h"
#include "cryptominisat5/solvertypesmini.h"
#include "primitive_types.h"

StopWatch::StopWatch()
{
  interval_length_.tv_sec = 60;
  gettimeofday(&last_interval_start_, NULL);
  start_time_ = stop_time_ = last_interval_start_;
}

timeval StopWatch::getElapsedTime()
{
  timeval r;
  timeval other_time = stop_time_;
  if (stop_time_.tv_sec == start_time_.tv_sec && stop_time_.tv_usec == start_time_.tv_usec)
    gettimeofday(&other_time, NULL);
  long int ad = 0;
  long int bd = 0;

  if (other_time.tv_usec < start_time_.tv_usec)
  {
    ad = 1;
    bd = 1000000;
  }
  r.tv_sec = other_time.tv_sec - ad - start_time_.tv_sec;
  r.tv_usec = other_time.tv_usec + bd - start_time_.tv_usec;
  return r;
}

void Solver::HardWireAndCompact()
{
  compactClauses();
  compactVariables();
  literal_stack_.clear();

  for (auto l = LiteralID(1, false); l != literals_.end_lit(); l.inc())
  {
    literal(l).activity_score_ = literal(l).binary_links_.size() - 1;
    literal(l).activity_score_ += occurrence_lists_[l].size();
  }
  statistics_.num_unit_clauses_ = unit_clauses_.size();
  initStack(num_variables());
  original_lit_pool_size_ = literal_pool_.size();
}

void Solver::solve(const string &file_name)
{
  srand(config_.randomseed);
  stopwatch_.start();
  createfromFile(file_name);
  if (solver.okay()) HardWireAndCompact();
  if (config_.perform_pcc) comp_manager_.getrandomseedforclhash();

  initStack(num_variables());
  if (!config_.quiet)
  {
    cout << "c Solving file " << file_name << endl;
    statistics_.printShortFormulaInfo();
    cout << "c Sampling set size: " << independent_support_.size() << endl;
    if (independent_support_.size() > 50) {
      cout << "c Sampling set is too large, not displaying" << endl;
    } else {
      cout << "c Sampling set: ";
      for (const auto& i: independent_support_) cout << ' ' << i;
      cout << endl;
    }
  }

  if (solver.okay()) {
    if (!config_.quiet) statistics_.printShortFormulaInfo();
    last_ccl_deletion_decs_ = last_ccl_cleanup_decs_ =
      statistics_.getNumDecisions();
    comp_manager_.initialize(literals_, literal_pool_, num_variables());

    statistics_.exit_state_ = countSAT();
    statistics_.set_final_solution_count_projected(decision_stack_.top().getTotalModelCount());
    statistics_.num_long_conflict_clauses_ = num_conflict_clauses();
  } else {
    cout << "c Found UNSAT during preprocessing" << endl;
    statistics_.exit_state_ = SUCCESS;
    statistics_.set_final_solution_count(0);
  }

  stopwatch_.stop();
  statistics_.time_elapsed_ = stopwatch_.getElapsedSeconds();
  comp_manager_.gatherStatistics();
  statistics_.printShort();
}

bool Solver::takeSolution() {
  //solver.set_polarity_mode(CMSat::PolarityMode::polarmode_rnd);
  //solver.set_up_for_sample_counter(100);
  CMSat::lbool ret = solver.solve();
  assert(ret != CMSat::l_Undef);
  if (ret == CMSat::l_False) {
    cout << "CMS gave UNSAT" << endl;
    return false;
  }
  for(uint32_t i = 0; i < num_variables(); i++) {
    target_polar[i+1] = solver.get_model()[i] == CMSat::l_True;
  }
  counted_bottom_component = false;
  return ret == CMSat::l_True;
}

SOLVER_StateT Solver::countSAT() {
  retStateT state = RESOLVED;

  if (!takeSolution()) return SUCCESS;
  while (true) {
    //print_debug("var top of decision stack: " << decision_stack_.top().getbranchvar());
    //NOTE: findNextRemainingComponentOf finds disjoing components!
    while (comp_manager_.findNextRemainingComponentOf(decision_stack_.top())) {
      checkProbabilisticHashSanity();
      decideLiteral();
      while (!failedLitProbe()) {
        state = resolveConflict();
        if (state == BACKTRACK) break;
      }
      if (state == BACKTRACK) break;
    }
    if (!counted_bottom_component) {
      assert(statistics_.num_decisions_ >= statistics_.last_restart_decisions);
      print_debug("Bottom component reached, decisions since restart: "
        << statistics_.num_decisions_ - statistics_.last_restart_decisions);
      counted_bottom_component = true;
    }

    state = backtrack();
    if (state == RESTART) continue;
    if (state == EXIT) return SUCCESS;
    while (state != PROCESS_COMPONENT && !failedLitProbe()) {
      state = resolveConflict();
      if (state == BACKTRACK) {
        state = backtrack();
        if (state == EXIT) return SUCCESS;
      }
    }
  }
  return SUCCESS;
}

void Solver::decideLiteral() {
  // establish another decision stack level
  decision_stack_.push_back(
    StackLevel(decision_stack_.top().currentRemainingComponent(),
               literal_stack_.size(),
               comp_manager_.component_stack_size()));

  auto it = comp_manager_.superComponentOf(decision_stack_.top()).varsBegin();
  unsigned max_score_var = *it;
  float max_score = scoreOf(*(it));
  float score;

  while (*it != varsSENTINEL &&
           independent_support_.find(*it) == independent_support_.end()) {
    it++;
  }

  if (*it != varsSENTINEL) {
    max_score_var = *it;
    max_score = scoreOf(*it);
  }

  while (*it != varsSENTINEL) {
    if (independent_support_.find(*it) != independent_support_.end()) {
      score = scoreOf(*it);
      if (score > max_score) {
        max_score = score;
        max_score_var = *it;
      }
    }
    it++;
  }

  // Decision scores
  /*if (config_.use_csvsads) {
    float cachescore = comp_manager_.cacheScoreOf(max_score_var);
    for (auto it = comp_manager_.superComponentOf(decision_stack_.top()).varsBegin();
         *it != varsSENTINEL; it++) {
      if (independent_support_.find(*it) != independent_support_.end()) {
        score = scoreOf(*it);
        if (score > max_score * config_.csvsads_param) {
          if (comp_manager_.cacheScoreOf(*it) > cachescore) {
            max_score_var = *it;
            cachescore = comp_manager_.cacheScoreOf(*it);
          }
        }
      }
    }
  }*/

  // this assert should always hold,
  // if not then there is a bug in the logic of countSAT();
  assert(max_score_var != 0);
  bool polarity;
  if (!counted_bottom_component) polarity = target_polar[max_score_var];
  else switch (config_.polarity_config) {
    case polar_default:
      polarity = literal(LiteralID(max_score_var, true)).activity_score_ > literal(LiteralID(max_score_var, false)).activity_score_;
      break;
    case polaritycache:
      polarity = literal(LiteralID(max_score_var, true)).activity_score_ >
        literal(LiteralID(max_score_var, false)).activity_score_;
      if (literal(LiteralID(max_score_var, true)).activity_score_ >
            2 * literal(LiteralID(max_score_var, false)).activity_score_) {
        polarity = true;
      } else if (literal(LiteralID(max_score_var, false)).activity_score_ >
                  2 * literal(LiteralID(max_score_var, true)).activity_score_) {
        polarity = false;
      } else if (var(max_score_var).set) {
        int random = rand() % 3;
        switch (random) {
          case 0:
            polarity = literal(LiteralID(max_score_var, true)).activity_score_ >
              literal(LiteralID(max_score_var, false)).activity_score_;
            break;
          case 1:
            polarity = var(max_score_var).polarity;
            break;
          case 2:
            polarity = var(max_score_var).polarity;
            polarity = !polarity;
            break;
        }
      }
      break;
    default:
      assert(false);
      exit(-1);
  }
  LiteralID theLit(max_score_var, polarity);
  decision_stack_.top().setbranchvariable(max_score_var);
  decision_stack_.top().setonpath(!counted_bottom_component);

  print_debug("Deciding lit: " << theLit << " dec level: " << decision_stack_.get_decision_level());
  setLiteralIfFree(theLit);
  statistics_.num_decisions_++;
  if (statistics_.num_decisions_ % 128 == 0) {
    if (config_.use_csvsads) {
      comp_manager_.increasecachescores();
    }
    decayActivities();
  }
  assert(
      decision_stack_.top().remaining_components_ofs() <= comp_manager_.component_stack_size());
  if (decision_stack_.get_decision_level() > statistics_.max_decision_level_) {
    statistics_.max_decision_level_ = decision_stack_.get_decision_level();
    if (statistics_.max_decision_level_ % 25 == 0) {
      cout << "c Max decision level :" << statistics_.max_decision_level_ << endl;
    }
  }
}

retStateT Solver::backtrack() {
  assert(
      decision_stack_.top().remaining_components_ofs() <= comp_manager_.component_stack_size());

  //Restart
  if (statistics_.getNumDecisions() > statistics_.next_restart) {
    statistics_.num_restarts++;
    statistics_.next_restart_diff*=1.4;
    statistics_.next_restart += statistics_.next_restart_diff;
    if ((statistics_.num_restarts % 10) == 5) {
      statistics_.next_restart_diff = 1000;
    }
    statistics_.last_restart_decisions = statistics_.num_decisions_;
    cout << "Restart here" << endl;
    if (counted_bottom_component) {
      //smallest cube is valid.
      vector<CMSat::Lit> cl;
      for(const auto&l: smallest_cube) cl.push_back(~CMSat::Lit(abs(l)-1, l<0));
      solver.add_clause(cl);
      cout << "cube: ";
      for(const auto&l: smallest_cube) cout << l << " ";
      cout << endl;
      counted_bottom_component = false;
    }
    do {
      if (decision_stack_.top().branch_found_unsat() || decision_stack_.top().anotherCompProcessible()) {
        comp_manager_.removeAllCachePollutionsOf(decision_stack_.top());
      }
      reactivateTOS();
      decision_stack_.pop_back();
    } while (decision_stack_.get_decision_level() > 0);
    if (!takeSolution()) return EXIT;
    return RESTART;
  }

  do {
    if (decision_stack_.top().branch_found_unsat()) {
      comp_manager_.removeAllCachePollutionsOf(decision_stack_.top());
    } else if (decision_stack_.top().anotherCompProcessible()) {
      print_debug("Processing another component at dec lev " << decision_stack_.get_decision_level()
          << " instead of bakctracking." << " Num unprocessed components: "
          << decision_stack_.top().numUnprocessedComponents());
      return PROCESS_COMPONENT;
    }

    if (!decision_stack_.top().isSecondBranch()) {
      print_debug("isSecondBranch (i.e. active branch is FALSE)"
          << " -- dec lev: " << decision_stack_.get_decision_level());
      const LiteralID aLit = TOS_decLit();
      assert(decision_stack_.get_decision_level() > 0);
      decision_stack_.top().changeBranch();
      reactivateTOS();
      print_debug("I don't understand... setLiteralIfFree called with NOT_A_CLAUSE -- it's in the literal stack of the decision stack's top element...");
      setLiteralIfFree(aLit.neg(), NOT_A_CLAUSE);
      return RESOLVED;
    } else {
      print_debug("not isSecondBranch (i.e. active branch is TRUE)"
          << " -- dec lev: " << decision_stack_.get_decision_level());
    }
    comp_manager_.cacheModelCountOf(decision_stack_.top().super_component(),
                                    decision_stack_.top().getTotalModelCount());

    // Update cache score heuristic
    if (config_.use_csvsads) {
      statistics_.numcachedec_++;
      if (statistics_.numcachedec_ % 128 == 0) {
        comp_manager_.increasecachescores();
      }
      comp_manager_.decreasecachescore(comp_manager_.superComponentOf(decision_stack_.top()));
    }

    // Backtrack from end, i.e. finished.
    if (decision_stack_.get_decision_level() == 0) {
      print_debug("Backtracking from lev 0, i.e. ending");
      break;
    }

    reactivateTOS();
    assert(decision_stack_.size() >= 2);
    (decision_stack_.end() - 2)->includeSolution(decision_stack_.top().getTotalModelCount());
    print_debug("Backtracking from level " << decision_stack_.get_decision_level()
        << " count here is: " << decision_stack_.top().getTotalModelCount());
    decision_stack_.pop_back();
    print_debug("-> Backtracked to level " << decision_stack_.get_decision_level()
        << " num unprocessed components here: " << decision_stack_.top().numUnprocessedComponents()
        << " on_path: " << decision_stack_.top().on_path_to_target_);
    if (decision_stack_.top().on_path_to_target_) {
      int at = 0;
#ifdef VERBOSE_DEBUG
      cout << "Smallest cube so far. Size: " << decision_stack_.size()-1 << " cube:";
      for(const auto& d: decision_stack_) {
        if (at > 0) cout << (target_polar[d.getbranchvar()] ? 1 : -1)*(int)d.getbranchvar() << " ";
        at++;
      }
      cout << endl;
#endif
      smallest_cube.clear();
      at = 0;
      for(const auto& d: decision_stack_) {
        if (at > 0) smallest_cube.push_back(
            (target_polar[d.getbranchvar()] ? 1 : -1)*(int)d.getbranchvar());;
        at++;
      }
    }
    // step to the next component not yet processed
    decision_stack_.top().nextUnprocessedComponent();

    assert(
        decision_stack_.top().remaining_components_ofs() < comp_manager_.component_stack_size() + 1);
  } while (true);
  return EXIT;
}

retStateT Solver::resolveConflict() {
  recordLastUIPCauses();

  if (statistics_.num_clauses_learned_ - last_ccl_deletion_decs_ >
        statistics_.clause_deletion_interval()) {
    deleteConflictClauses();
    last_ccl_deletion_decs_ = statistics_.num_clauses_learned_;
  }

  if (statistics_.num_clauses_learned_ - last_ccl_cleanup_decs_ > 100000) {
    compactConflictLiteralPool();
    last_ccl_cleanup_decs_ = statistics_.num_clauses_learned_;
  }

  statistics_.num_conflicts_++;

  assert(
      decision_stack_.top().remaining_components_ofs() <= comp_manager_.component_stack_size());

  assert(uip_clauses_.size() == 1);

  if (uip_clauses_.back().size() == 0) {
    cout << "c EMPTY CLAUSE FOUND" << endl;
  }

  decision_stack_.top().mark_branch_unsat();
  //BEGIN Backtracking
  // maybe the other branch had some solutions
  if (decision_stack_.top().isSecondBranch()) {
    if (decision_stack_.get_decision_level() == 1) {
      cout
          << "c Solved half the solution space (i.e. one branch at dec. lev 1)." << endl
          << "c --> Conflicts: " << statistics_.num_conflicts_ << endl
          << "c --> Decisions: " << statistics_.num_decisions_ << endl;
    }
    return BACKTRACK;
  }

  Antecedent ant(NOT_A_CLAUSE);
  // this has to be checked since using implicit BCP
  // and checking literals there not exhaustively
  // we cannot guarantee that uip_clauses_.back().front() == TOS_decLit().neg()
  // this is because we might have checked a literal
  // during implict BCP which has been a failed literal
  // due only to assignments made at lower decision levels
  if (uip_clauses_.back().front() == TOS_decLit().neg()) {
    assert(TOS_decLit().neg() == uip_clauses_.back()[0]);
    var(TOS_decLit().neg()).ante = addUIPConflictClause(
        uip_clauses_.back());
    ant = var(TOS_decLit()).ante;
  }
  assert(decision_stack_.get_decision_level() > 0);
  assert(decision_stack_.top().branch_found_unsat());

  // we do not have to remove pollutions here,
  // since conflicts only arise directly before
  // remaining components are stored
  // hence
  assert(
      decision_stack_.top().remaining_components_ofs() == comp_manager_.component_stack_size());

  decision_stack_.top().changeBranch();
  const LiteralID lit = TOS_decLit();
  reactivateTOS();
  if (ant == NOT_A_CLAUSE) {
    print_debug("Conflict pushes us to: " << lit << " and due to failed literal probling, we can't guarantee it's due to the 1UIP, so setting it as a decision instead");
  } else {
    print_debug("Conflict pushes us to: " << lit);
  }
  setLiteralIfFree(lit.neg(), ant);
  //END Backtracking
  return RESOLVED;
}

bool Solver::failedLitProbe() {
  // the asserted literal has been set, so we start
  // bcp on that literal
  assert(literal_stack_.size() > 0 &&
      "We could just put this in an IF, but I don't think it should be 0");
  unsigned start_ofs = literal_stack_.size() - 1;
  print_debug("--> Setting units of this component...");
  for (const auto& lit : unit_clauses_) setLiteralIfFree(lit);
  print_debug("--> Units of this component set, propagating");
  bool bSucceeded = propagate(start_ofs);
  if (config_.perform_failed_lit_test && bSucceeded) {
    bSucceeded = failedLitProbeInternal();
  }
  return bSucceeded;
}

bool Solver::propagate(unsigned start_at_stack_ofs) {
  for (unsigned int i = start_at_stack_ofs; i < literal_stack_.size(); i++) {
    const LiteralID unLit = literal_stack_[i].neg();

    //Propagate bin Clauses
    for (auto bt = literal(unLit).binary_links_.begin();
         *bt != SENTINEL_LIT; bt++) {
      if (isResolved(*bt)) {
        setConflictState(unLit, *bt);
        return false;
      }
      setLiteralIfFree(*bt, Antecedent(unLit));
    }

    //Propagate long clauses
    for (auto itcl = literal(unLit).watch_list_.rbegin();
         *itcl != SENTINEL_CL; itcl++) {
      bool isLitA = (*beginOf(*itcl) == unLit);
      auto p_watchLit = beginOf(*itcl) + 1 - isLitA;
      auto p_otherLit = beginOf(*itcl) + isLitA;

      if (isSatisfied(*p_otherLit)) {
        continue;
      }
      auto itL = beginOf(*itcl) + 2;
      while (isResolved(*itL)) {
        itL++;
      }
      // either we found a free or satisfied lit
      if (*itL != SENTINEL_LIT) {
        literal(*itL).addWatchLinkTo(*itcl);
        swap(*itL, *p_watchLit);
        *itcl = literal(unLit).watch_list_.back();
        literal(unLit).watch_list_.pop_back();
      } else {
        // or p_unLit stays resolved
        // and we have hence no free literal left
        // for p_otherLit remain poss: Active or Resolved
        if (setLiteralIfFree(*p_otherLit, Antecedent(*itcl))) { // implication
          if (isLitA) {
            swap(*p_otherLit, *p_watchLit);
          }
        } else {
          setConflictState(*itcl);
          return false;
        }
      }
    }
  }
  return true;
}

// this is IBCP 30.08
bool Solver::failedLitProbeInternal() {
  static vector<LiteralID> test_lits(num_variables());
  static LiteralIndexedVector<unsigned char> viewed_lits(
      num_variables() + 1, 0);
  unsigned stack_ofs = decision_stack_.top().literal_stack_ofs();
  unsigned num_curr_lits = 0;
  while (stack_ofs < literal_stack_.size()) {
    test_lits.clear();
    for (auto it = literal_stack_.begin() + stack_ofs;
         it != literal_stack_.end(); it++) {
      for (auto cl_ofs : occurrence_lists_[it->neg()]) {
        if (!isSatisfied(cl_ofs)) {
          for (auto lt = beginOf(cl_ofs); *lt != SENTINEL_LIT; lt++) {
            if (isActive(*lt) && !viewed_lits[lt->neg()]) {
              test_lits.push_back(lt->neg());
              viewed_lits[lt->neg()] = true;
            }
          }
        }
      }
    }
    num_curr_lits = literal_stack_.size() - stack_ofs;
    stack_ofs = literal_stack_.size();
    for (auto jt = test_lits.begin(); jt != test_lits.end(); jt++) {
      viewed_lits[*jt] = false;
    }

    // Figure out which literals to probe
    vector<float> scores;
    scores.clear();
    for (auto jt = test_lits.begin(); jt != test_lits.end(); jt++) {
      scores.push_back(literal(*jt).activity_score_);
    }
    sort(scores.begin(), scores.end());
    num_curr_lits = 10 + num_curr_lits / 20;
    float threshold = 0.0;
    if (scores.size() > num_curr_lits) {
      threshold = scores[scores.size() - num_curr_lits];
    }
    statistics_.num_failed_literal_tests_ += test_lits.size();

    // Do the probing
    print_debug("Failed literal probing START");
    for (auto lit : test_lits) {
      if (isActive(lit) && threshold <= literal(lit).activity_score_) {
        unsigned sz = literal_stack_.size();
        // we increase the decLev artificially
        // s.t. after the tentative BCP call, we can learn a conflict clause
        // relative to the assignment of *jt
        decision_stack_.startFailedLitTest();
        setLiteralIfFree(lit);

        assert(!hasAntecedent(lit));

        bool bSucceeded = propagate(sz);
        if (!bSucceeded) recordAllUIPCauses();
        decision_stack_.stopFailedLitTest();

        while (literal_stack_.size() > sz) {
          unSet(literal_stack_.back());
          literal_stack_.pop_back();
        }

        if (!bSucceeded) {
          statistics_.num_failed_literals_detected_++;
          sz = literal_stack_.size();
          for (auto it = uip_clauses_.rbegin();
               it != uip_clauses_.rend(); it++) {
            if (it->size() == 0) cout << "c EMPTY CLAUSE FOUND" << endl;
            setLiteralIfFree(it->front(),
                             addUIPConflictClause(*it));
          }
          if (!propagate(sz)) {
            print_debug("Failed literal probing END -- this component/branch is UNSAT");
            return false;
          }
        }
      }
    }
  }
  print_debug("Failed literal probing END -- no UNSAT");
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN module conflictAnalyzer
///////////////////////////////////////////////////////////////////////////////////////////////

void Solver::minimizeAndStoreUIPClause(
  LiteralID uipLit,
  vector<LiteralID> &tmp_clause, const vector<unsigned char>& seen) {

  static deque<LiteralID> clause;
  clause.clear();
  assertion_level_ = 0;
  for (auto lit : tmp_clause) {
    if (existsUnitClauseOf(lit.var())) {
      continue;
    }
    bool resolve_out = false;
    if (hasAntecedent(lit)) {
      resolve_out = true;
      if (getAntecedent(lit).isAClause()) {
        for (auto it = beginOf(getAntecedent(lit).asCl()) + 1;
             *it != SENTINEL_CL; it++) {
          if (!seen[it->var()]) {
            resolve_out = false;
            break;
          }
        }
      } else if (!seen[getAntecedent(lit).asLit().var()]) {
        resolve_out = false;
      }
    }

    if (!resolve_out) {
      // uipLit should be the sole literal of this Decision Level
      if (var(lit).decision_level >= assertion_level_) {
        assertion_level_ = var(lit).decision_level;
        clause.push_front(lit);
      } else {
        clause.push_back(lit);
      }
    }
  }

  if (uipLit.var()) {
    assert(var(uipLit).decision_level >= 0
            && (unsigned)var(uipLit).decision_level == decision_stack_.get_decision_level());
  }

  //assert(uipLit.var() != 0);
  if (uipLit.var() != 0) {
    clause.push_front(uipLit);
  }
  uip_clauses_.push_back(vector<LiteralID>(clause.begin(), clause.end()));
}

void Solver::recordLastUIPCauses() {
  // note:
  // variables of lower dl: if seen we dont work with them anymore
  // variables of this dl: if seen we incorporate their
  // antecedent and set to unseen
  tmp_seen.clear();
  tmp_seen.resize(num_variables()+1, false);
  tmp_clause.clear();
  assert(toClear.empty());

  assertion_level_ = 0;
  uip_clauses_.clear();

  unsigned lit_stack_ofs = literal_stack_.size();
  const unsigned DL = decision_stack_.get_decision_level();
  unsigned lits_at_current_dl = 0;

  for (const auto& l: violated_clause) {
    if (var(l).decision_level == 0 || existsUnitClauseOf(l.var())) {
      continue;
    }
    if (var(l).decision_level < (int)DL) {
      tmp_clause.push_back(l);
    } else {
      lits_at_current_dl++;
    }
    literal(l).increaseActivity();
    tmp_seen[l.var()] = true;
    toClear.push_back(l.var());
  }

  LiteralID curr_lit;
  while (lits_at_current_dl) {
    assert(lit_stack_ofs != 0);
    curr_lit = literal_stack_[--lit_stack_ofs];

    if (!tmp_seen[curr_lit.var()]) {
      continue;
    }
    tmp_seen[curr_lit.var()] = false;
    toClear.push_back(curr_lit.var());

    if (lits_at_current_dl-- == 1) {
      // perform UIP stuff
      if (!hasAntecedent(curr_lit)) {
        // this should be the decision literal when in first branch
        // or it is a literal decided to explore in failed literal testing
        break;
      }
    }

    assert(hasAntecedent(curr_lit));
    if (getAntecedent(curr_lit).isAClause()) {
      updateActivities(getAntecedent(curr_lit).asCl());
      assert(curr_lit == *beginOf(getAntecedent(curr_lit).asCl()));

      for (auto it = beginOf(getAntecedent(curr_lit).asCl()) + 1;
           *it != SENTINEL_CL; it++) {
        if (tmp_seen[it->var()] || (var(*it).decision_level == 0) || existsUnitClauseOf(it->var())) {
          continue;
        }
        if (var(*it).decision_level < (int)DL) {
          tmp_clause.push_back(*it);
        } else {
          lits_at_current_dl++;
        }
        tmp_seen[it->var()] = true;
        toClear.push_back(it->var());
      }
    } else {
      LiteralID alit = getAntecedent(curr_lit).asLit();
      literal(alit).increaseActivity();
      literal(curr_lit).increaseActivity();
      if (!tmp_seen[alit.var()] && !(var(alit).decision_level == 0) &&
            !existsUnitClauseOf(alit.var())) {
        if (var(alit).decision_level < (int)DL) {
          tmp_clause.push_back(alit);
        } else {
          lits_at_current_dl++;
        }
        tmp_seen[alit.var()] = true;
        toClear.push_back(alit.var());
      }
    }
    curr_lit = NOT_A_LIT;
  }
  minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, tmp_seen);
  for(const auto& v: toClear) tmp_seen[v] = false;
  toClear.clear();
}

void Solver::recordAllUIPCauses() {
  // note:
  // variables of lower dl: if seen we dont work with them anymore
  // variables of this dl: if seen we incorporate their
  // antecedent and set to unseen
  tmp_seen.clear();
  tmp_seen.resize(num_variables()+1, false);
  tmp_clause.clear();
  assert(toClear.empty());

  assertion_level_ = 0;
  uip_clauses_.clear();

  unsigned lit_stack_ofs = literal_stack_.size();
  const unsigned DL = decision_stack_.get_decision_level();
  unsigned lits_at_current_dl = 0;

  for (auto l : violated_clause) {
    if (var(l).decision_level == 0 || existsUnitClauseOf(l.var())) {
      continue;
    }
    if (var(l).decision_level < (int)DL) {
      tmp_clause.push_back(l);
    } else {
      lits_at_current_dl++;
    }
    literal(l).increaseActivity();
    tmp_seen[l.var()] = true;
    toClear.push_back(l.var());
  }
  unsigned n = 0;
  LiteralID curr_lit;
  while (lits_at_current_dl) {
    assert(lit_stack_ofs != 0);
    curr_lit = literal_stack_[--lit_stack_ofs];

    if (!tmp_seen[curr_lit.var()]) {
      continue;
    }

    tmp_seen[curr_lit.var()] = false;

    if (lits_at_current_dl-- == 1) {
      n++;
      if (!hasAntecedent(curr_lit)) {
        // this should be the decision literal when in first branch
        // or it is a literal decided to explore in failed literal testing
        //assert(stack_.TOS_decLit() == curr_lit);
        break;
      }
      // perform UIP stuff
      minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, tmp_seen);
    }

    assert(hasAntecedent(curr_lit));

    if (getAntecedent(curr_lit).isAClause()) {
      updateActivities(getAntecedent(curr_lit).asCl());
      assert(curr_lit == *beginOf(getAntecedent(curr_lit).asCl()));

      for (auto it = beginOf(getAntecedent(curr_lit).asCl()) + 1;
           *it != SENTINEL_CL; it++) {
        if (tmp_seen[it->var()] || (var(*it).decision_level == 0) ||
              existsUnitClauseOf(it->var())) {
          continue;
        }
        if (var(*it).decision_level < (int)DL) {
          tmp_clause.push_back(*it);
        } else {
          lits_at_current_dl++;
        }
        tmp_seen[it->var()] = true;
        toClear.push_back(it->var());
      }
    } else {
      LiteralID alit = getAntecedent(curr_lit).asLit();
      literal(alit).increaseActivity();
      literal(curr_lit).increaseActivity();
      if (!tmp_seen[alit.var()] && !(var(alit).decision_level == 0) &&
            !existsUnitClauseOf(alit.var())) {
        if (var(alit).decision_level < (int)DL) {
          tmp_clause.push_back(alit);
        } else {
          lits_at_current_dl++;
        }
        tmp_seen[alit.var()] = true;
        toClear.push_back(alit.var());
      }
    }
  }
  if (!hasAntecedent(curr_lit)) {
    minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, tmp_seen);
  }
  for(const auto& v: toClear) tmp_seen[v] = false;
  toClear.clear();
}

void Solver::printOnlineStats() {
  if (config_.quiet) {
    return;
  }

  cout << endl;
  cout << "time elapsed: " << stopwatch_.getElapsedSeconds() << "s" << endl;
  if (config_.verbose) {
    cout << "conflict clauses (all / bin / unit) \t";
    cout << num_conflict_clauses();
    cout << "/" << statistics_.num_binary_conflict_clauses_ << "/"
         << unit_clauses_.size() << endl;
    cout << "failed literals found by implicit BCP \t "
         << statistics_.num_failed_literals_detected_ << endl;

    cout << "implicit BCP miss rate \t "
         << statistics_.implicitBCP_miss_rate() * 100 << "%";
    cout << endl;

    comp_manager_.gatherStatistics();

    cout << "cache size " << statistics_.cache_MB_memory_usage() << "MB" << endl;
    cout << "components (stored / hits) \t\t"
         << statistics_.cached_component_count() << "/"
         << statistics_.cache_hits() << endl;
    cout << "avg. variable count (stored / hits) \t"
         << statistics_.getAvgComponentSize() << "/"
         << statistics_.getAvgCacheHitSize();
    cout << endl;
    cout << "cache miss rate " << statistics_.cache_miss_rate() * 100 << "%"
         << endl;
  }
}
