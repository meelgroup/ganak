/*
 * solver.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: marc
 */
#include "solver.h"

#include <algorithm>
#include <ios>
#include <iomanip>
#include "common.h"
#include "comp_types/comp.h"
#include "cryptominisat5/solvertypesmini.h"
#include "primitive_types.h"
#include "stack.h"
#include "structures.h"
#include "time_mem.h"

bool Solver::simplePreProcess()
{
  uint32_t start_ofs = 0;

  for (auto lit : unit_clauses_) {
    if (isUnitClause(lit.neg())) return false;
    setLiteralIfFree(lit);
  }

  bool succeeded = propagate(start_ofs);
  if (succeeded) HardWireAndCompact();
  return succeeded;
}

void Solver::HardWireAndCompact()
{
  compactClauses();
  compactVariables();
  trail.clear();
  test_lits.resize(num_variables());
  viewed_lits.resize(num_variables() + 1, 0);

  for (auto l = Lit(1, false); l != literals_.end_lit(); l.inc())
  {
    litWatchList(l).activity_score_ = litWatchList(l).binary_links_.size() - 1;
    litWatchList(l).activity_score_ += occ_lists_[l].size();
  }
  stats.num_unit_clauses_ = unit_clauses_.size();
  initStack();
  original_lit_pool_size_ = literal_pool_.size();
}

void Solver::solve(const std::string &file_name)
{
  srand(config_.randomseed);
  time_start = cpuTime();
  createfromFile(file_name);
  if (config_.perform_pcc) comp_manager_.getrandomseedforclhash();

  initStack();
  if (config_.verb) {
    cout << "c Solving file " << file_name << endl;
    stats.printShortFormulaInfo();
    if (indep_support_given) {
      cout << "c Sampling set size: " << indep_support_.size() << endl;
      if (indep_support_.size() > 50) {
        cout << "c Sampling set is too large, not displaying" << endl;
      } else {
        cout << "c Sampling set: ";
        for (const auto& i: indep_support_) cout << ' ' << i;
        cout << endl;
      }
    } else cout << "c No sampling set, doing unprojected counting" << endl;
  }

  if (satSolver.okay() && !simplePreProcess()) {
    stats.exit_state_ = SUCCESS;
    stats.set_final_solution_count(0);
  }
  if (config_.verb) cout << "c Prepocessing done" << endl;
  if (satSolver.okay()) {
    if (config_.verb) stats.printShortFormulaInfo();
    last_ccl_deletion_decs_ = last_ccl_cleanup_decs_ = stats.getNumDecisions();
    comp_manager_.initialize(literals_, literal_pool_, num_variables());

    stats.exit_state_ = countSAT();
    stats.set_final_solution_count_projected(decision_stack_.top().getTotalModelCount());
    stats.num_long_conflict_clauses_ = num_conflict_clauses();
  } else {
    cout << "c Found UNSAT during preprocessing" << endl;
    stats.exit_state_ = SUCCESS;
    stats.set_final_solution_count(0);
  }

  stats.time_elapsed_ = time_start - cpuTime();
  comp_manager_.gatherStatistics();
  stats.printShort();
}

bool Solver::takeSolution() {
  //solver.set_polarity_mode(CMSat::PolarityMode::polarmode_rnd);
  //solver.set_up_for_sample_counter(100);
  CMSat::lbool ret = satSolver.solve();
  assert(ret != CMSat::l_Undef);
  if (ret == CMSat::l_False) {
    cout << "c CMS gave UNSAT" << endl;
    return false;
  }
  assert(ret == CMSat::l_True);
  for(uint32_t i = 0; i < num_variables(); i++) {
    target_polar[i+1] = satSolver.get_model()[i] == CMSat::l_True;
  }
  counted_bottom_comp = false;
  return true;
}

SOLVER_StateT Solver::countSAT() {
  retStateT state = RESOLVED;

  counted_bottom_comp = true;
  if (config_.restart && !takeSolution()) return SUCCESS;
  while (true) {
    print_debug("var top of decision stack: " << decision_stack_.top().getbranchvar());
    //NOTE: findNextRemainingComponentOf finds disjoint comps!
    while (comp_manager_.findNextRemainingComponentOf(decision_stack_.top())) {
      checkProbabilisticHashSanity();
      decideLiteral();

#ifdef VERBOSE_DEBUG
      cout << COLORG "--- going through all levels now, printing comps --" << endl;
      uint32_t lev = 0;
      for(const auto& s: decision_stack_) {
        auto const& sup_at = s.super_comp();
        cout << COLORG "super comp of lev " << lev
          << " is at " << sup_at
          << " branch var here: " << decision_stack_.at(lev).getbranchvar()
          << " remaining comp ofs: " << decision_stack_.at(lev).remaining_comps_ofs()
          << " num unprocess comps: " << decision_stack_.at(lev).numUnprocessedComponents()
          << endl;

        const auto& c = comp_manager_.at(sup_at);
        cout << COLORG "-> Variables in comp_manager_.at(" << sup_at << ")."
          << " num: " << c->num_variables() << " vars: ";
        for(uint32_t i = 0; i < c->num_variables(); i++) {
          const auto& v = c->varsBegin();
          cout << v[i] << " ";
        }
        cout << endl;
        lev++;
      }
      cout << COLORG "--- Went through all levels now --" << endl;
#endif

      while (!failedLitProbe()) {
        state = resolveConflict();
        if (state == BACKTRACK) break;
      }
      if (state == BACKTRACK) break;
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
  print_debug("new decision level is about to be created, lev now: " << decision_stack_.get_decision_level());
  decision_stack_.push_back(
    StackLevel(decision_stack_.top().currentRemainingComponent(),
               trail.size(),
               comp_manager_.comp_stack_size()));

  // Find variable to branch on
  auto it = comp_manager_.getSuperComponentOf(decision_stack_.top()).varsBegin();
  uint32_t max_score_var = *it;
  float max_score = scoreOf(*(it));
  float score;
  while (*it != varsSENTINEL &&
           indep_support_.find(*it) == indep_support_.end()) {
    it++;
  }
  if (*it != varsSENTINEL) {
    max_score_var = *it;
    max_score = scoreOf(*it);
  }
  while (*it != varsSENTINEL) {
    if (indep_support_.find(*it) != indep_support_.end()) {
      score = scoreOf(*it);
      if (score > max_score) {
        max_score = score;
        max_score_var = *it;
      }
    }
    it++;
  }
  assert(max_score_var != 0 &&
        "this assert should always hold, if not then there is a bug in the logic of countSAT()");

  // Figure out polarity
  bool polarity;
  if (!counted_bottom_comp) polarity = target_polar[max_score_var];
  else {
    // TODO MATE: this whole thing is a huge mess as far as I'm concerned
    polarity = litWatchList(Lit(max_score_var, true)).activity_score_ >
      litWatchList(Lit(max_score_var, false)).activity_score_;
    if (litWatchList(Lit(max_score_var, true)).activity_score_ >
          2 * litWatchList(Lit(max_score_var, false)).activity_score_) {
      polarity = true;
    } else if (litWatchList(Lit(max_score_var, false)).activity_score_ >
                2 * litWatchList(Lit(max_score_var, true)).activity_score_) {
      polarity = false;
    } else if (var(max_score_var).set) {
      // TODO MATE this sounds insane, right? Random polarities??
      uint32_t random = mtrand.randInt(2) ;
      switch (random) {
        case 0:
          polarity = litWatchList(Lit(max_score_var, true)).activity_score_ >
            litWatchList(Lit(max_score_var, false)).activity_score_;
          break;
        case 1:
          polarity = var(max_score_var).polarity;
          break;
        case 2:
          polarity = !(var(max_score_var).polarity);
          break;
      }
    }
  }

  // The decision literal is now ready. Deal with it.
  const Lit lit(max_score_var, polarity);
  print_debug(COLYEL "decideLiteral() is deciding: " << lit << " dec level: "
      << decision_stack_.get_decision_level());
  decision_stack_.top().setbranchvariable(max_score_var);
  decision_stack_.top().setonpath(!counted_bottom_comp);
  setLiteralIfFree(lit);
  stats.num_decisions_++;
  if (stats.num_decisions_ % 128 == 0) {
    if (config_.use_csvsads) comp_manager_.increasecachescores();
    decayActivities();
  }
  assert( decision_stack_.top().remaining_comps_ofs() <= comp_manager_.comp_stack_size());
}

void Solver::computeLargestCube()
{
  largest_cube.clear();
  print_debug(COLWHT "-- computeLargestCube BEGIN");

  // add decisions
  print_debug_noendl(COLWHT << "dec vars: ");
  for(uint32_t i = 1; i < decision_stack_.size(); i++) {
    const StackLevel& ds = decision_stack_[i];
    const auto dec_lit = (target_polar[ds.getbranchvar()] ? 1 : -1)*(int)ds.getbranchvar();
    print_debug_noendl(dec_lit << " ");
    largest_cube.push_back(dec_lit);
  }
  print_debug_noendl(endl);

  // Show decision stack's comps
#ifdef VERBOSE_DEBUG
  for(size_t i = 0; i < decision_stack_.size(); i++) {
    const auto& ds = decision_stack_.at(i);
    const auto dec_lit = (target_polar[ds.getbranchvar()] ? 1 : -1)*(int)ds.getbranchvar();
    print_debug(COLWHT "decision_stack.at " << i
      << " decision lit: " << dec_lit
      << " num unproc comps: " << ds.numUnprocessedComponents()
      << " unproc comps end: " << ds.getUnprocessedComponentsEnd()
      << " remain comps offs: " << ds.remaining_comps_ofs());
    const auto off_start = ds.remaining_comps_ofs();
    const auto off_end = ds.getUnprocessedComponentsEnd();
    for(uint32_t i2 = off_start; i2 < off_end; i2++) {
      assert(i2 < comp_manager_.comp_stack_size());
      const auto& c = comp_manager_.at(i2);
      cout << COLWHT "-> comp at: " << std::setw(3) << i2 << " ID: " << c->id() << " -- vars : ";
      for(auto v = c->varsBegin(); *v != varsSENTINEL; v++) cout << *v << " ";
      cout << endl;
    }
  }

  // All comps
  print_debug(COLWHT "-- comp list START");
  for(uint32_t i2 = 0; i2 < comp_manager_.comp_stack_size(); i2++) {
    const auto& c = comp_manager_.at(i2);
    cout << COLWHT "comp at: " << std::setw(3) << i2 << " ID: " << c->id() << " -- vars : ";
    if (c->empty()) { cout << "EMPTY" << endl; continue; }
    for(auto v = c->varsBegin(); *v != varsSENTINEL; v++) cout << *v << " ";
    cout << endl;
  }
  print_debug(COLWHT "-- comp list END");

   cout << COLWHT "Largest cube so far. Size: " << largest_cube.size() << " cube: ";
   for(const auto& l: largest_cube) cout << l << " ";
   cout << endl;
#endif
}

retStateT Solver::backtrack() {
  assert(
      decision_stack_.top().remaining_comps_ofs() <= comp_manager_.comp_stack_size());

  //Restart
  if (config_.restart && stats.getNumDecisions() > stats.next_restart) {
    stats.num_restarts++;
    stats.next_restart_diff*=1.4;
    stats.next_restart += stats.next_restart_diff;
    if ((stats.num_restarts % 5) == 4) {
      stats.next_restart_diff = 1000;
    }
    stats.last_restart_decisions = stats.num_decisions_;
    cout << "c Restart here" << endl;
    if (counted_bottom_comp) {
      //largest cube is valid.
      vector<CMSat::Lit> cl;
      for(const auto&l: largest_cube) cl.push_back(~CMSat::Lit(abs(l)-1, l<0));
      satSolver.add_clause(cl);
      cout << "c cube: ";
      for(const auto&l: largest_cube) cout << l << " ";
      cout << endl;
      counted_bottom_comp = false;
    }
    do {
      if (decision_stack_.top().branch_found_unsat() ||
          decision_stack_.top().anotherCompProcessible()) {
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
      print_debug("Processing another comp at dec lev "
          << decision_stack_.get_decision_level()
          << " instead of bakctracking." << " Num unprocessed comps: "
          << decision_stack_.top().numUnprocessedComponents());
      return PROCESS_COMPONENT;
    }

    // We have NOT explored the other side! Let's do it now!
    if (!decision_stack_.top().isSecondBranch()) {
      print_debug("We have NOT explored the right branch (isSecondBranch==false). Let's do it now."
          << " -- dec lev: " << decision_stack_.get_decision_level());
      //top of stack decision lit
      const Lit aLit = TOS_decLit();
      assert(decision_stack_.get_decision_level() > 0);
      decision_stack_.top().changeBranch(); //flip branch
      reactivateTOS();
      print_debug("Flipping lit to: " << aLit.neg());
      setLiteralIfFree(aLit.neg(), NOT_A_CLAUSE);
      print_debug(COLORGBG "Backtrack finished -- we flipped the branch");
      return RESOLVED;
    } else {
      print_debug(COLORGBG "We have explored BOTH branches, actually BACKTRACKING."
          << " -- dec lev: " << decision_stack_.get_decision_level());
    }
    comp_manager_.cacheModelCountOf(decision_stack_.top().super_comp(),
                                    decision_stack_.top().getTotalModelCount());

    // Update cache score heuristic
    if (config_.use_csvsads) {
      stats.numcachedec_++;
      if (stats.numcachedec_ % 128 == 0) comp_manager_.increasecachescores();
      comp_manager_.decreasecachescore(comp_manager_.getSuperComponentOf(decision_stack_.top()));
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
        << " num unprocessed comps here: " << decision_stack_.top().numUnprocessedComponents()
        << " on_path: " << decision_stack_.top().on_path_to_target_);

    if (decision_stack_.top().on_path_to_target_) {
      computeLargestCube();
      if (!counted_bottom_comp) {
        assert(stats.num_decisions_ >= stats.last_restart_decisions);
        print_debug(COLCYN "Bottom comp reached, decisions since restart: "
          << stats.num_decisions_ - stats.last_restart_decisions);
        counted_bottom_comp = true;
      }
    }
    // step to the next comp not yet processed
    decision_stack_.top().nextUnprocessedComponent();

    assert(
        decision_stack_.top().remaining_comps_ofs() < comp_manager_.comp_stack_size() + 1);
  } while (true);
  return EXIT;
}

retStateT Solver::resolveConflict() {
  recordLastUIPCauses();

  if (stats.num_clauses_learned_ - last_ccl_deletion_decs_ > stats.clause_deletion_interval()) {
    deleteConflictClauses();
    last_ccl_deletion_decs_ = stats.num_clauses_learned_;
  }

  if (stats.num_clauses_learned_ - last_ccl_cleanup_decs_ > 100000) {
    compactConflictLiteralPool();
    last_ccl_cleanup_decs_ = stats.num_clauses_learned_;
  }

  stats.num_conflicts_++;
  assert(decision_stack_.top().remaining_comps_ofs() <= comp_manager_.comp_stack_size());
  assert(uip_clauses_.size() == 1);
  if (uip_clauses_.back().empty()) { cout << "c EMPTY CLAUSE FOUND" << endl; }

  decision_stack_.top().mark_branch_unsat();
  //Backtracking
  // maybe the other branch had some solutions
  if (decision_stack_.top().isSecondBranch()) {
    if (decision_stack_.get_decision_level() == 1) {
      cout
          << "c Solved half the solution space (i.e. one branch at dec. lev 1)." << endl
          << "c --> Conflicts: " << stats.num_conflicts_ << endl
          << "c --> Decisions: " << stats.num_decisions_ << endl;
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
  if (!uip_clauses_.back().empty() && uip_clauses_.back().front() == TOS_decLit().neg()) {
    assert(TOS_decLit().neg() == uip_clauses_.back()[0]);
    var(TOS_decLit().neg()).ante = addUIPConflictClause(
        uip_clauses_.back());
    ant = var(TOS_decLit()).ante;
  }
  assert(decision_stack_.get_decision_level() > 0);
  assert(decision_stack_.top().branch_found_unsat());

  // we do not have to remove pollutions here,
  // since conflicts only arise directly before
  // remaining comps are stored
  // hence
  assert( decision_stack_.top().remaining_comps_ofs() == comp_manager_.comp_stack_size());

  decision_stack_.top().changeBranch();
  const Lit lit = TOS_decLit();
  reactivateTOS();
  if (ant == NOT_A_CLAUSE) {
    print_debug("Conflict pushes us to: " << lit<< " and due to failed literal probling, we can't guarantee it's due to the 1UIP, so setting it as a decision instead");
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
  assert(trail.size() > 0 && "This is FISHY, I fixed with a bad hack, using 'was_at_zero' but was it broken before? I think we may need propagating from 0 in these cases? Or not?");

  const bool was_at_zero = trail.size() == 0;
  const uint32_t start_ofs = trail.size() - 1;
  print_debug("--> Setting units of this comp...");
  for (const auto& lit : unit_clauses_) setLiteralIfFree(lit);
  print_debug("--> Units of this comp set, propagating");

  bool bSucceeded = true;
  if (!was_at_zero) bSucceeded = propagate(start_ofs);
  if (config_.perform_failed_lit_test && bSucceeded) {
    bSucceeded = failedLitProbeInternal();
  }
  return bSucceeded;
}

bool Solver::propagate(const uint32_t start_at_stack_ofs) {
  for (auto i = start_at_stack_ofs; i < trail.size(); i++) {
    const Lit unLit = trail[i].neg();

    //Propagate bin clauses
    for (auto bt = litWatchList(unLit).binary_links_.begin(); *bt != SENTINEL_LIT; bt++) {
      if (isFalse(*bt)) {
        setConflictState(unLit, *bt);
        return false;
      }
      setLiteralIfFree(*bt, Antecedent(unLit));
    }

    //Propagate long clauses
    for (auto itcl = litWatchList(unLit).watch_list_.rbegin(); *itcl != SENTINEL_CL; itcl++) {
      bool isLitA = (*beginOf(*itcl) == unLit);
      auto p_watchLit = beginOf(*itcl) + 1 - isLitA;
      auto p_otherLit = beginOf(*itcl) + isLitA;

      if (isTrue(*p_otherLit)) continue;
      auto itL = beginOf(*itcl) + 2;
      while (isFalse(*itL)) itL++;
      // either we found a free or satisfied lit
      if (*itL != SENTINEL_LIT) {
        litWatchList(*itL).addWatchLinkTo(*itcl);
        std::swap(*itL, *p_watchLit);
        *itcl = litWatchList(unLit).watch_list_.back();
        litWatchList(unLit).watch_list_.pop_back();
      } else {
        // or p_unLit stays resolved
        // and we have hence no free literal left
        // for p_otherLit remain poss: Active or Resolved
        if (setLiteralIfFree(*p_otherLit, Antecedent(*itcl))) { // implication
          if (isLitA) std::swap(*p_otherLit, *p_watchLit);
        } else {
          setConflictState(*itcl);
          return false;
        }
      }
    }
  }
  return true;
}

bool Solver::failedLitProbeInternal() {
  print_debug(COLRED "Failed literal probing START");

  uint32_t stack_ofs = decision_stack_.top().literal_stack_ofs();
  uint32_t num_curr_lits = 0;
  while (stack_ofs < trail.size()) {
    test_lits.clear();
    for (auto it = trail.begin() + stack_ofs;
         it != trail.end(); it++) {
      for (auto cl_ofs : occ_lists_[it->neg()]) {
        if (!isSatisfied(cl_ofs)) {
          for (auto lt = beginOf(cl_ofs); *lt != SENTINEL_LIT; lt++) {
            if (isUnknown(*lt) && !viewed_lits[lt->neg()]) {
              test_lits.push_back(lt->neg());
              print_debug("-> potential lit to test: " << lt->neg());
              viewed_lits[lt->neg()] = true;
            }
          }
        }
      }
    }
    num_curr_lits = trail.size() - stack_ofs;
    stack_ofs = trail.size();
    for (auto jt = test_lits.begin(); jt != test_lits.end(); jt++) {
      viewed_lits[*jt] = false;
    }

    // Figure out which literals to probe
    vector<float> scores;
    scores.clear();
    for (auto jt = test_lits.begin(); jt != test_lits.end(); jt++) {
      scores.push_back(litWatchList(*jt).activity_score_);
    }
    sort(scores.begin(), scores.end());
    num_curr_lits = 10 + num_curr_lits / 20;
    float threshold = 0.0;
    if (scores.size() > num_curr_lits) {
      threshold = scores[scores.size() - num_curr_lits];
    }
    stats.num_failed_literal_tests_ += test_lits.size();

    // Do the probing
    for (auto lit : test_lits) {
      if (isUnknown(lit) && threshold <= litWatchList(lit).activity_score_) {
        uint32_t sz = trail.size();
        // we increase the decLev artificially
        // s.t. after the tentative BCP call, we can learn a conflict clause
        // relative to the assignment of *jt
        decision_stack_.startFailedLitTest();
        setLiteralIfFree(lit);

        assert(!hasAntecedent(lit));

        bool bSucceeded = propagate(sz);
        if (!bSucceeded) recordAllUIPCauses();
        decision_stack_.stopFailedLitTest();

        // backtracking
        while (trail.size() > sz) {
          unSet(trail.back());
          trail.pop_back();
        }

        if (!bSucceeded) {
          stats.num_failed_literals_detected_++;
          print_debug("-> failed literal detected");
          sz = trail.size();
          for (auto it = uip_clauses_.rbegin();
               it != uip_clauses_.rend(); it++) {
            if (it->size() == 0) cout << "c EMPTY CLAUSE FOUND" << endl;
            setLiteralIfFree(it->front(),
                             addUIPConflictClause(*it));
          }
          if (!propagate(sz)) {
            print_debug("Failed literal probing END -- this comp/branch is UNSAT");
            return false;
          }
        }
      }
    }
  }
  print_debug(COLRED "Failed literal probing END -- no UNSAT, gotta check this branch");
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN module conflictAnalyzer
///////////////////////////////////////////////////////////////////////////////////////////////

void Solver::minimizeAndStoreUIPClause(
  Lit uipLit,
  vector<Lit> &tmp_clause, const vector<uint8_t>& seen) {

  tmp_clause_minim.clear();
  assertion_level_ = 0;
  for (auto lit : tmp_clause) {
    if (existsUnitClauseOf(lit.var())) continue;
    bool resolve_out = false;
    if (hasAntecedent(lit)) {
      resolve_out = true;
      if (getAntecedent(lit).isAClause()) {
        for (auto it = beginOf(getAntecedent(lit).asCl()) + 1; *it != SENTINEL_CL; it++) {
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
        tmp_clause_minim.push_front(lit);
      } else {
        tmp_clause_minim.push_back(lit);
      }
    }
  }

  if (uipLit.var()) {
    assert(var(uipLit).decision_level >= 0
            && (uint32_t)var(uipLit).decision_level == decision_stack_.get_decision_level());
  }

  //assert(uipLit.var() != 0);
  if (uipLit.var() != 0) {
    tmp_clause_minim.push_front(uipLit);
  }
  uip_clauses_.push_back(vector<Lit>(tmp_clause_minim.begin(), tmp_clause_minim.end()));
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

  uint32_t trail_ofs = trail.size();
  const uint32_t DL = decision_stack_.get_decision_level();
  uint32_t lits_at_current_dl = 0;

  for (const auto& l: violated_clause) {
    if (var(l).decision_level == 0 || existsUnitClauseOf(l.var())) {
      continue;
    }
    if (var(l).decision_level < (int)DL) tmp_clause.push_back(l);
    else lits_at_current_dl++;
    litWatchList(l).increaseActivity();
    tmp_seen[l.var()] = true;
    toClear.push_back(l.var());
  }

  Lit curr_lit;
  while (lits_at_current_dl) {
    assert(trail_ofs != 0);
    curr_lit = trail[--trail_ofs];

    if (!tmp_seen[curr_lit.var()]) continue;
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
      Lit alit = getAntecedent(curr_lit).asLit();
      litWatchList(alit).increaseActivity();
      litWatchList(curr_lit).increaseActivity();
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

  uint32_t trail_ofs = trail.size();
  const uint32_t DL = decision_stack_.get_decision_level();
  uint32_t lits_at_current_dl = 0;

  for (const auto& l : violated_clause) {
    if (var(l).decision_level == 0 || existsUnitClauseOf(l.var())) continue;
    if (var(l).decision_level < (int)DL) tmp_clause.push_back(l);
    else lits_at_current_dl++;
    litWatchList(l).increaseActivity();
    tmp_seen[l.var()] = true;
    toClear.push_back(l.var());
  }
  Lit curr_lit;
  while (lits_at_current_dl) {
    assert(trail_ofs != 0);
    curr_lit = trail[--trail_ofs];
    if (!tmp_seen[curr_lit.var()]) continue;
    tmp_seen[curr_lit.var()] = false;

    if (lits_at_current_dl-- == 1) {
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
      Lit alit = getAntecedent(curr_lit).asLit();
      litWatchList(alit).increaseActivity();
      litWatchList(curr_lit).increaseActivity();
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
  if (config_.verb == 0) return;

  cout << "c " << endl;
  cout << "c time elapsed: " << time_start - cpuTime() << "s" << endl;
  if (config_.verb >= 2) {
    cout << "conflict clauses (all / bin / unit) \t";
    cout << num_conflict_clauses();
    cout << "/" << stats.num_binary_conflict_clauses_ << "/"
         << unit_clauses_.size() << endl;
    cout << "failed literals found by implicit BCP \t "
         << stats.num_failed_literals_detected_ << endl;

    cout << "implicit BCP miss rate \t "
         << stats.implicitBCP_miss_rate() * 100 << "%";
    cout << endl;

    comp_manager_.gatherStatistics();

    cout << "cache size " << stats.cache_MB_memory_usage() << "MB" << endl;
    cout << "comps (stored / hits) \t\t"
         << stats.cached_comp_count() << "/"
         << stats.cache_hits() << endl;
    cout << "avg. variable count (stored / hits) \t"
         << stats.getAvgComponentSize() << "/"
         << stats.getAvgCacheHitSize();
    cout << endl;
    cout << "cache miss rate " << stats.cache_miss_rate() * 100 << "%"
         << endl;
  }
}
