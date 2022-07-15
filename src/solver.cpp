/*
 * solver.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: marc
 */
#include "solver.h"
#include <deque>

#include <algorithm>

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

void Solver::print(vector<LiteralID> &vec)
{
  cout << "c ";
  for (auto l : vec)
    cout << l.toInt() << " ";
  cout << endl;
}

void Solver::print(vector<unsigned> &vec)
{
  cout << "c ";
  for (auto l : vec)
    cout << l << " ";
  cout << endl;
}

bool Solver::simplePreProcess()
{

  if (!config_.perform_pre_processing)
    return true;

  assert(literal_stack_.size() == 0);
  unsigned start_ofs = 0;

  //begin process unit clauses
  for (auto lit : unit_clauses_)
  {
    if (isUnitClause(lit.neg()))
    {
      return false;
    }
    setLiteralIfFree(lit);
  }
  //end process unit clauses

  bool succeeded = BCP(start_ofs);

  if (succeeded)
    succeeded &= prepFailedLiteralTest();

  if (succeeded)
    HardWireAndCompact();

  return succeeded;
}

bool Solver::prepFailedLiteralTest()
{
  unsigned last_size;
  do
  {
    last_size = literal_stack_.size();
    for (unsigned v = 1; v < variables_.size(); v++)
      if (isActive(v))
      {
        unsigned sz = literal_stack_.size();
        setLiteralIfFree(LiteralID(v, true));
        bool res = BCP(sz);
        while (literal_stack_.size() > sz)
        {
          unSet(literal_stack_.back());
          literal_stack_.pop_back();
        }

        if (!res)
        {
          sz = literal_stack_.size();
          setLiteralIfFree(LiteralID(v, false));
          if (!BCP(sz))
            return false;
        }
        else
        {

          sz = literal_stack_.size();
          setLiteralIfFree(LiteralID(v, false));
          bool resb = BCP(sz);
          while (literal_stack_.size() > sz)
          {
            unSet(literal_stack_.back());
            literal_stack_.pop_back();
          }
          if (!resb)
          {
            sz = literal_stack_.size();
            setLiteralIfFree(LiteralID(v, true));
            if (!BCP(sz))
              return false;
          }
        }
      }
  } while (literal_stack_.size() > last_size);

  return true;
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

  statistics_.num_original_binary_clauses_ = statistics_.num_binary_clauses_;
  statistics_.num_original_unit_clauses_ = statistics_.num_unit_clauses_ =
      unit_clauses_.size();
  initStack(num_variables());
  original_lit_pool_size_ = literal_pool_.size();
}

void Solver::solve(const string &file_name)
{
  stopwatch_.setTimeBound(config_.time_bound_seconds);
  srand(config_.randomseed);
  stopwatch_.start();
  statistics_.input_file_ = file_name;

  bool notfoundUNSAT = createfromFile(file_name);

  //Found Empirically
  if (statistics_.num_original_binary_clauses_ > 0.75 * statistics_.num_original_clauses_)
  {
    config_.maxdecterminate = false;
  }

  if (config_.perform_pcc)
  {
    comp_manager_.getrandomseedforclhash();
  }

  initStack(num_variables());

  if (!config_.quiet)
  {
    cout << "c Solving " << file_name << endl;
    statistics_.printShortFormulaInfo();
  }
  if (independent_support_.size() == 0)
  {
    if (!config_.quiet)
    {
      cout << "c Sampling set not present! So doing total model counting." << endl;
    }
    config_.perform_projectedmodelcounting = false;
  }
  else if (!config_.quiet)
  {
    if (!config_.perform_projectedmodelcounting)
    {
      cout << "c Warning! Sampling set is present but projected model counting"
           << " is turned off by the user so solver is not doing projected model counting." << endl;
    }
    else
    {
      cout << "c Sampling set is present, performing projected model counting " << endl;
    }
    cout << "c Sampling set size: " << independent_support_.size() << endl;
    cout << "c Sampling set: ";
    for (auto it = independent_support_.begin(); it != independent_support_.end(); ++it)
    {
      cout << ' ' << *it;
    }
    cout << endl;
  }
  if (!config_.quiet)
  {
    cout << "c " << endl;
    cout << "c Preprocessing .." << endl;
  }

  if (notfoundUNSAT) {
    notfoundUNSAT = simplePreProcess();
  }

  if (!config_.quiet) {
    cout << "c DONE" << endl;
  }

  if (notfoundUNSAT) {

    if (!config_.quiet) {
      statistics_.printShortFormulaInfo();
    }

    last_ccl_deletion_time_ = last_ccl_cleanup_time_ =
      statistics_.getTime();

    violated_clause.reserve(num_variables());

    comp_manager_.initialize(literals_, literal_pool_, num_variables());

    statistics_.exit_state_ = countSAT();

    if (statistics_.exit_state_ == CHANGEHASH) {
      cout << "ERROR: We need to change the hash range (-1)" << endl;
      exit(1);
    }
    if (config_.perform_projectedmodelcounting) {
      statistics_.set_final_solution_count_projected(stack_.top().getTotalModelCount());
    } else {
      statistics_.set_final_solution_count(stack_.top().getTotalModelCount());
    }
    statistics_.num_long_conflict_clauses_ = num_conflict_clauses();
  } else {
    statistics_.exit_state_ = SUCCESS;
    statistics_.set_final_solution_count(0.0);
    cout << endl
         << "c FOUND UNSAT DURING PREPROCESSING " << endl;
  }

  stopwatch_.stop();
  statistics_.time_elapsed_ = stopwatch_.getElapsedSeconds();

  comp_manager_.gatherStatistics();
  string writefile;
  if (config_.perform_projectedmodelcounting) {
    writefile = "out.pmc";
  } else {
    writefile = "out.mc";
  }
//  statistics_.writeToFile(writefile, config_.perform_projectedmodelcounting);
  statistics_.printShort(config_.perform_projectedmodelcounting);
}

SOLVER_StateT Solver::countSAT() {
  retStateT state = RESOLVED;

  while (true) {
    while (comp_manager_.findNextRemainingComponentOf(stack_.top())) {
      unsigned t = statistics_.num_cache_look_ups_ + 1;
      if (2 * log2(t) > log2(config_.delta) + 64 * config_.hashrange * 0.9843) { // 1 - log_2(2.004)/64 = 0.9843
        return CHANGEHASH;
      }
      decideLiteral();
      if (stopwatch_.timeBoundBroken()) {
        return TIMEOUT;
      }
      while (!bcp()) {
        state = resolveConflict();
        if (state == BACKTRACK) {
          break;
        }
      }
      if (state == BACKTRACK) {
        break;
      }
    }

    state = backtrack();
    if (state == RESTART) {
      continue;
    }
    else if (state == EXIT) {
      return SUCCESS;
    }
    while (state != PROCESS_COMPONENT && !bcp()) {
      state = resolveConflict();
      if (state == BACKTRACK) {
        state = backtrack();
        if (state == EXIT) {
          return SUCCESS;
        }
      }
    }
  }
  return SUCCESS;
}

void Solver::decideLiteral() {
  // establish another decision stack level
  stack_.push_back(
    StackLevel(stack_.top().currentRemainingComponent(),
               literal_stack_.size(),
               comp_manager_.component_stack_size()));

  auto it = comp_manager_.superComponentOf(stack_.top()).varsBegin();
  unsigned max_score_var = *it;
  float max_score = scoreOf(*(it));
  float score;
  if (config_.perform_projectedmodelcounting)
  {
    isindependent = true;
    bool isindependent_support_present = false;
    while (*it != varsSENTINEL &&
             independent_support_.find(*it) == independent_support_.end()) {
      it++;
    }
    if (*it != varsSENTINEL) {
      isindependent_support_present = true;
      max_score_var = *it;
      max_score = scoreOf(*it);
    }
    while (*it != varsSENTINEL) {
      if (independent_support_.find(*it) != independent_support_.end()) {
        isindependent_support_present = true;
        score = scoreOf(*it);
        if (score > max_score) {
          max_score = score;
          max_score_var = *it;
        }
      }
      it++;
    }
    if (isindependent_support_present && config_.use_csvsads) {
      float cachescore = comp_manager_.cacheScoreOf(max_score_var);
      for (auto it = comp_manager_.superComponentOf(stack_.top()).varsBegin();
           *it != varsSENTINEL; it++) {
        if (independent_support_.find(*it) != independent_support_.end()) {
          score = scoreOf(*it);
          if (score > max_score * config_.csvsads_param) {
            if (comp_manager_.cacheScoreOf(*it) > cachescore) {
              isindependent_support_present = true;
              max_score_var = *it;
              cachescore = comp_manager_.cacheScoreOf(*it);
            }
          }
        }
      }
      max_score = score;
    } else if (!isindependent_support_present) {
      isindependent = false;
      max_score = -1;
      score = -1;
      for (auto it = comp_manager_.superComponentOf(stack_.top()).varsBegin();
           *it != varsSENTINEL; it++) {
        score = scoreOf(*it);
        if (score > max_score) {
          max_score = score;
          max_score_var = *it;
        }
      }
      if (config_.use_csvsads) {
        float cachescore = comp_manager_.cacheScoreOf(max_score_var);
        for (auto it = comp_manager_.superComponentOf(stack_.top()).varsBegin();
             *it != varsSENTINEL; it++) {
          score = scoreOf(*it);
          if (score > max_score * config_.csvsads_param) {
            if (comp_manager_.cacheScoreOf(*it) > cachescore) {
              max_score_var = *it;
              cachescore = comp_manager_.cacheScoreOf(*it);
            }
          }
        }
        max_score = score;
      }
    }
  } else if (config_.use_edr) {
    vector<unsigned> bucket;
    double he_param = 1.0 - 0.1 * exp(-0.0001 * statistics_.num_decisions_);
    if (statistics_.num_decisions_ > 20000) {
      config_.use_edr = false;
    }
    for (auto it = comp_manager_.superComponentOf(stack_.top()).varsBegin();
         *it != varsSENTINEL; it++) {
      score = scoreOf(*it);
      if (score > max_score) {
        max_score = score;
      }
    }
    if (max_score == 0) {
      max_score = -1;
    }
    for (auto it =
             comp_manager_.superComponentOf(stack_.top()).varsBegin();
         *it != varsSENTINEL; it++) {
      score = scoreOf(*it);
      if ((1.0 * score) / max_score > he_param) {
        bucket.push_back(*it);
      }
    }
    max_score_var = bucket[rand() % bucket.size()];
  } else {
    for (auto it =
             comp_manager_.superComponentOf(stack_.top()).varsBegin();
         *it != varsSENTINEL; it++) {
      score = scoreOf(*it);
      if (score > max_score) {
        max_score = score;
        max_score_var = *it;
      }
    }
    if (config_.use_csvsads) {
      float cachescore = comp_manager_.cacheScoreOf(max_score_var);
      for (auto it = comp_manager_.superComponentOf(stack_.top()).varsBegin();
           *it != varsSENTINEL; it++) {
        score = scoreOf(*it);
        if (score > max_score * config_.csvsads_param) {
          if (comp_manager_.cacheScoreOf(*it) > cachescore) {
            max_score_var = *it;
            cachescore = comp_manager_.cacheScoreOf(*it);
          }
        }
      }
      max_score = score;
    }
  }
  // this assert should always hold,
  // if not then there is a bug in the logic of countSAT();
  assert(max_score_var != 0);
  bool polarity = true;
  switch (config_.polarity_config) {
    case polar_false:
      polarity = false;
      break;
    case polar_true:
      polarity = true;
      break;
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
  }
  LiteralID theLit(max_score_var, polarity);
  stack_.top().setbranchvariable(max_score_var);

  setLiteralIfFree(theLit);
  statistics_.num_decisions_++;
  if (config_.maxdecterminate) {
    if (statistics_.num_decisions_ > config_.maxdec &&
        statistics_.num_conflicts_ < config_.minconflicts_) {
      cout << "c Terminating solver because the number of decisions exceeds the given value of " << config_.maxdec << " and conflicts is less than " << config_.minconflicts_ << endl;
      exit(1);
    }
  }
  if (statistics_.num_decisions_ % 128 == 0) {
    if (config_.use_csvsads) {
      comp_manager_.increasecachescores();
    }
    decayActivities();
  }
  assert(
      stack_.top().remaining_components_ofs() <= comp_manager_.component_stack_size());
  if (stack_.get_decision_level() > statistics_.max_decision_level_) {
    statistics_.max_decision_level_ = stack_.get_decision_level();
    if (statistics_.max_decision_level_ % 25 == 0) {
      cout << "c Max decision level :" << statistics_.max_decision_level_ << endl;
    }
  }
}

retStateT Solver::backtrack() {
  assert(
      stack_.top().remaining_components_ofs() <= comp_manager_.component_stack_size());

  if (config_.use_lso && statistics_.num_decisions_ >= config_.lsoafterdecisions) {
    config_.use_lso = false;
    if (!config_.quiet) {
      cout << "c Doing Restart" << endl;
    }
    do {
      if (stack_.top().branch_found_unsat() || stack_.top().anotherCompProcessible()) {
        comp_manager_.removeAllCachePollutionsOf(stack_.top());
      }
      reactivateTOS();
      stack_.pop_back();
    } while (stack_.get_decision_level() > 0);
    statistics_.num_decisions_ = 0;
    return RESTART;
  }
  if (!isindependent && config_.perform_projectedmodelcounting) {
    do {
      if (stack_.top().branch_found_unsat()) {
        comp_manager_.removeAllCachePollutionsOf(stack_.top());
      } else if (stack_.top().anotherCompProcessible()) {
        return PROCESS_COMPONENT;
      }
      if (stack_.top().getBranchSols() != 0 && isindependent == false) {
        while (independent_support_.count(stack_.top().getbranchvar()) == 0) {
          if (stack_.get_decision_level() <= 0) {
            break;
          }
          reactivateTOS();
          assert(stack_.size() >= 2);
          (stack_.end() - 2)->includeSolution(stack_.top().getTotalModelCount());
          stack_.pop_back();
          // step to the next component not yet processed
          stack_.top().nextUnprocessedComponent();
          assert(
              stack_.top().remaining_components_ofs() < comp_manager_.component_stack_size() + 1);
          if (stack_.top().anotherCompProcessible()) {
            return PROCESS_COMPONENT;
          }
        }
      }
//       isindependent = true;
      if (!stack_.top().isSecondBranch()) {
        LiteralID aLit = TOS_decLit();
        assert(stack_.get_decision_level() > 0);
        stack_.top().changeBranch();
        reactivateTOS();
        setLiteralIfFree(aLit.neg(), NOT_A_CLAUSE);
        return RESOLVED;
      }
      comp_manager_.cacheModelCountOf(stack_.top().super_component(),
                                      stack_.top().getTotalModelCount());
      if (config_.use_csvsads) {
        statistics_.numcachedec_++;
        if (statistics_.numcachedec_ % 128 == 0) {
          comp_manager_.increasecachescores();
        }
        comp_manager_.decreasecachescore(comp_manager_.superComponentOf(stack_.top()));
      }
      if (stack_.get_decision_level() <= 0) {
        break;
      }
      reactivateTOS();
      assert(stack_.size() >= 2);
      (stack_.end() - 2)->includeSolution(stack_.top().getTotalModelCount());
      stack_.pop_back();
      // step to the next component not yet processed
      stack_.top().nextUnprocessedComponent();

      assert(
          stack_.top().remaining_components_ofs() < comp_manager_.component_stack_size() + 1);
    } while (stack_.get_decision_level() >= 0);
  } else {
    do {
      if (stack_.top().branch_found_unsat()) {
        comp_manager_.removeAllCachePollutionsOf(stack_.top());
      } else if (stack_.top().anotherCompProcessible()) {
        return PROCESS_COMPONENT;
      }

      if (!stack_.top().isSecondBranch()) {
        if (stack_.get_decision_level() == 1) {
          cout << "c We have solved halfed" << endl;
          config_.maxdecterminate = false;
          config_.use_lso = false;
        }
        LiteralID aLit = TOS_decLit();
        assert(stack_.get_decision_level() > 0);
        if (stack_.get_decision_level() == 1) {
          cout << "c We have solved halfed" << endl;
          config_.use_lso = false;
        }
        stack_.top().changeBranch();
        reactivateTOS();
        setLiteralIfFree(aLit.neg(), NOT_A_CLAUSE);
        return RESOLVED;
      }
      // OTHERWISE:  backtrack further
      comp_manager_.cacheModelCountOf(stack_.top().super_component(),
                                      stack_.top().getTotalModelCount());

      //Cache score should be decreased since the component is getting added to cache
      if (config_.use_csvsads) {
        statistics_.numcachedec_++;
        if (statistics_.numcachedec_ % 128 == 0) {
          comp_manager_.increasecachescores();
        }
        comp_manager_.decreasecachescore(comp_manager_.superComponentOf(stack_.top()));
      }

      if (stack_.get_decision_level() <= 0) {
        break;
      }
      reactivateTOS();

      assert(stack_.size() >= 2);
      (stack_.end() - 2)->includeSolution(stack_.top().getTotalModelCount());
      stack_.pop_back();
      // step to the next component not yet processed
      stack_.top().nextUnprocessedComponent();

      assert(
          stack_.top().remaining_components_ofs() < comp_manager_.component_stack_size() + 1);

    } while (stack_.get_decision_level() >= 0);
  }
  return EXIT;
}

retStateT Solver::resolveConflict() {
  recordLastUIPCauses();

  if (statistics_.num_clauses_learned_ - last_ccl_deletion_time_ >
        statistics_.clause_deletion_interval()) {
    deleteConflictClauses();
    last_ccl_deletion_time_ = statistics_.num_clauses_learned_;
  }

  if (statistics_.num_clauses_learned_ - last_ccl_cleanup_time_ > 100000) {
    compactConflictLiteralPool();
    last_ccl_cleanup_time_ = statistics_.num_clauses_learned_;
  }

  statistics_.num_conflicts_++;

  assert(
      stack_.top().remaining_components_ofs() <= comp_manager_.component_stack_size());

  assert(uip_clauses_.size() == 1);

  // DEBUG
  if (uip_clauses_.back().size() == 0) {
    cout << "c EMPTY CLAUSE FOUND" << endl;
  }
  // END DEBUG

  stack_.top().mark_branch_unsat();
  //BEGIN Backtracking
  // maybe the other branch had some solutions
  if (stack_.top().isSecondBranch()) {
    if (stack_.get_decision_level() == 1) {
      cout << "c We have solved halfed" << endl;
      config_.use_lso = false;
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
  //	// RRR
  //	else if(var(uip_clauses_.back().front()).decision_level
  //			< stack_.get_decision_level()
  //			&& assertion_level_ <  stack_.get_decision_level()){
  //         stack_.top().set_both_branches_unsat();
  //         return BACKTRACK;
  //	}
  //
  //
  //	// RRR
  assert(stack_.get_decision_level() > 0);
  assert(stack_.top().branch_found_unsat());

  // we do not have to remove pollutions here,
  // since conflicts only arise directly before
  // remaining components are stored
  // hence
  assert(
      stack_.top().remaining_components_ofs() == comp_manager_.component_stack_size());

  stack_.top().changeBranch();
  LiteralID lit = TOS_decLit();
  reactivateTOS();
  setLiteralIfFree(lit.neg(), ant);
  //END Backtracking
  return RESOLVED;
}

bool Solver::bcp() {
  // the asserted literal has been set, so we start
  // bcp on that literal
  unsigned start_ofs = literal_stack_.size() - 1;

  //BEGIN process unit clauses
  for (auto lit : unit_clauses_) {
    setLiteralIfFree(lit);
  }
  //END process unit clauses

  bool bSucceeded = BCP(start_ofs);

  if (config_.perform_failed_lit_test && bSucceeded) {
    bSucceeded = implicitBCP();
  }
  return bSucceeded;
}

bool Solver::BCP(unsigned start_at_stack_ofs) {
  for (unsigned int i = start_at_stack_ofs; i < literal_stack_.size(); i++) {
    LiteralID unLit = literal_stack_[i].neg();
    //BEGIN Propagate Bin Clauses
    for (auto bt = literal(unLit).binary_links_.begin();
         *bt != SENTINEL_LIT; bt++) {
      if (isResolved(*bt)) {
        setConflictState(unLit, *bt);
        return false;
      }
      setLiteralIfFree(*bt, Antecedent(unLit));
    }
    //END Propagate Bin Clauses
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
bool Solver::implicitBCP() {
  static vector<LiteralID> test_lits(num_variables());
  static LiteralIndexedVector<unsigned char> viewed_lits(num_variables() + 1,
                                                         0);

  unsigned stack_ofs = stack_.top().literal_stack_ofs();
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

    for (auto lit : test_lits) {
      if (isActive(lit) && threshold <= literal(lit).activity_score_) {
        unsigned sz = literal_stack_.size();
        // we increase the decLev artificially
        // s.t. after the tentative BCP call, we can learn a conflict clause
        // relative to the assignment of *jt
        stack_.startFailedLitTest();
        setLiteralIfFree(lit);

        assert(!hasAntecedent(lit));

        bool bSucceeded = BCP(sz);
        if (!bSucceeded) {
          recordAllUIPCauses();
        }

        stack_.stopFailedLitTest();

        while (literal_stack_.size() > sz) {
          unSet(literal_stack_.back());
          literal_stack_.pop_back();
        }

        if (!bSucceeded) {
          statistics_.num_failed_literals_detected_++;
          sz = literal_stack_.size();
          for (auto it = uip_clauses_.rbegin();
               it != uip_clauses_.rend(); it++) {
            // DEBUG
            if (it->size() == 0) {
              cout << "c EMPTY CLAUSE FOUND" << endl;
            }
            // END DEBUG
            setLiteralIfFree(it->front(),
                             addUIPConflictClause(*it));
          }
          if (!BCP(sz)) {
            return false;
          }
        }
      }
    }
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN module conflictAnalyzer
///////////////////////////////////////////////////////////////////////////////////////////////

void Solver::minimizeAndStoreUIPClause(
  LiteralID uipLit,
  vector<LiteralID> &tmp_clause, bool seen[]) {

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
    assert(var(uipLit).decision_level == stack_.get_decision_level());
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
  bool seen[num_variables() + 1];
  memset(seen, false, sizeof(bool) * (num_variables() + 1));

  static vector<LiteralID> tmp_clause;
  tmp_clause.clear();

  assertion_level_ = 0;
  uip_clauses_.clear();

  unsigned lit_stack_ofs = literal_stack_.size();
  int DL = stack_.get_decision_level();
  unsigned lits_at_current_dl = 0;

  for (auto l : violated_clause) {
    if (var(l).decision_level == 0 || existsUnitClauseOf(l.var())) {
      continue;
    }
    if (var(l).decision_level < DL) {
      tmp_clause.push_back(l);
    } else {
      lits_at_current_dl++;
    }
    literal(l).increaseActivity();
    seen[l.var()] = true;
  }

  LiteralID curr_lit;
  while (lits_at_current_dl) {
    assert(lit_stack_ofs != 0);
    curr_lit = literal_stack_[--lit_stack_ofs];

    if (!seen[curr_lit.var()]) {
      continue;
    }

    seen[curr_lit.var()] = false;

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
        if (seen[it->var()] || (var(*it).decision_level == 0) || existsUnitClauseOf(it->var())) {
          continue;
        }
        if (var(*it).decision_level < DL) {
          tmp_clause.push_back(*it);
        } else {
          lits_at_current_dl++;
        }
        seen[it->var()] = true;
      }
    } else {
      LiteralID alit = getAntecedent(curr_lit).asLit();
      literal(alit).increaseActivity();
      literal(curr_lit).increaseActivity();
      if (!seen[alit.var()] && !(var(alit).decision_level == 0) &&
            !existsUnitClauseOf(alit.var())) {
        if (var(alit).decision_level < DL) {
          tmp_clause.push_back(alit);
        } else {
          lits_at_current_dl++;
        }
        seen[alit.var()] = true;
      }
    }
    curr_lit = NOT_A_LIT;
  }
  minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, seen);

  //	if (var(curr_lit).decision_level > assertion_level_)
  //		assertion_level_ = var(curr_lit).decision_level;
}

void Solver::recordAllUIPCauses() {
  // note:
  // variables of lower dl: if seen we dont work with them anymore
  // variables of this dl: if seen we incorporate their
  // antecedent and set to unseen
  bool seen[num_variables() + 1];
  memset(seen, false, sizeof(bool) * (num_variables() + 1));

  static vector<LiteralID> tmp_clause;
  tmp_clause.clear();

  assertion_level_ = 0;
  uip_clauses_.clear();

  unsigned lit_stack_ofs = literal_stack_.size();
  int DL = stack_.get_decision_level();
  unsigned lits_at_current_dl = 0;

  for (auto l : violated_clause) {
    if (var(l).decision_level == 0 || existsUnitClauseOf(l.var())) {
      continue;
    }
    if (var(l).decision_level < DL) {
      tmp_clause.push_back(l);
    } else {
      lits_at_current_dl++;
    }
    literal(l).increaseActivity();
    seen[l.var()] = true;
  }
  unsigned n = 0;
  LiteralID curr_lit;
  while (lits_at_current_dl) {
    assert(lit_stack_ofs != 0);
    curr_lit = literal_stack_[--lit_stack_ofs];

    if (!seen[curr_lit.var()]) {
      continue;
    }

    seen[curr_lit.var()] = false;

    if (lits_at_current_dl-- == 1) {
      n++;
      if (!hasAntecedent(curr_lit)) {
        // this should be the decision literal when in first branch
        // or it is a literal decided to explore in failed literal testing
        //assert(stack_.TOS_decLit() == curr_lit);
        break;
      }
      // perform UIP stuff
      minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, seen);
    }

    assert(hasAntecedent(curr_lit));

    if (getAntecedent(curr_lit).isAClause()) {
      updateActivities(getAntecedent(curr_lit).asCl());
      assert(curr_lit == *beginOf(getAntecedent(curr_lit).asCl()));

      for (auto it = beginOf(getAntecedent(curr_lit).asCl()) + 1;
           *it != SENTINEL_CL; it++) {
        if (seen[it->var()] || (var(*it).decision_level == 0) ||
              existsUnitClauseOf(it->var())) {
          continue;
        }
        if (var(*it).decision_level < DL) {
          tmp_clause.push_back(*it);
        } else {
          lits_at_current_dl++;
        }
        seen[it->var()] = true;
      }
    } else {
      LiteralID alit = getAntecedent(curr_lit).asLit();
      literal(alit).increaseActivity();
      literal(curr_lit).increaseActivity();
      if (!seen[alit.var()] && !(var(alit).decision_level == 0) &&
            !existsUnitClauseOf(alit.var())) {
        if (var(alit).decision_level < DL) {
          tmp_clause.push_back(alit);
        } else {
          lits_at_current_dl++;
        }
        seen[alit.var()] = true;
      }
    }
  }
  if (!hasAntecedent(curr_lit)) {
    minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, seen);
  }
  //	if (var(curr_lit).decision_level > assertion_level_)
  //		assertion_level_ = var(curr_lit).decision_level;
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
