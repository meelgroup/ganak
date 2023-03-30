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

void Solver::simplePreProcess()
{
  for (auto lit : unit_clauses_) {
    assert(!isUnitClause(lit.neg()) && "Formula is not UNSAT, we ran CMS before");;
    setLiteralIfFree(lit);
  }

  bool succeeded = propagate(0);
  release_assert(succeeded && "We ran CMS before, so it cannot be UNSAT");
  viewed_lits.resize(nVars() + 1, 0);
  stats.num_unit_irred_clauses_ = unit_clauses_.size();
  irred_lit_pool_size_ = lit_pool_.size();
  init_decision_stack();
}

void Solver::set_indep_support(const set<uint32_t> &indeps)
{
  indep_support_ = indeps;
}

void Solver::end_irred_cls()
{
  release_assert(!ended_irred_cls && "ERROR *must not* call end_irred_cls() twice");
  stats.next_restart = config_.first_restart;
  stats.maximum_cache_size_bytes_ = config_.maximum_cache_size_bytes_;
  if (config_.do_pcc) comp_manager_.getrandomseedforclhash();
  init_decision_stack();
  simplePreProcess();
  ended_irred_cls = true;
}

void Solver::add_red_cl(const vector<Lit>& lits) {
  release_assert(ended_irred_cls && "ERROR *must* call end_irred_cls() before add_red_cl()");
  assert(lits.size() <= 2); // TODO longer clauses -- but then have to add activity
  addUIPConflictClause(lits);
}

void Solver::get_unit_cls(vector<Lit>& units) const
{
  assert(units.empty());
  units = unit_clauses_;
}

void Solver::get_bin_red_cls(vector<Lit>& bins) const
{
  for(size_t i = 2; i < watches_.size(); i++) {
    Lit l(i/2, i%2);
    const auto ws = watches_[l];
    for(auto i2 = ws.last_irred_bin; i2 < ws.binary_links_.size(); i2++) {
      const auto l2 = ws.binary_links_[i2];
      if (l2 == SENTINEL_LIT) continue;
      if (l2 > l) continue; //don't add it twice
      bins.push_back(l);
      bins.push_back(ws.binary_links_[i2]);
      bins.push_back(SENTINEL_LIT);
    }
  }
}

mpz_class Solver::solve(vector<Lit>& largest_cube_ret)
{
  release_assert(ended_irred_cls && "ERROR *must* call end_irred_cls() before solve()");
  time_start = cpuTime();

  if (config_.verb) {
      cout << "c Sampling set size: " << indep_support_.size() << endl;
      if (indep_support_.size() > 50) {
        cout << "c Sampling set is too large, not displaying" << endl;
      } else {
        cout << "c Sampling set: ";
        for (const auto& i: indep_support_) cout << ' ' << i;
        cout << endl;
      }
  }
  if (config_.verb) stats.printShortFormulaInfo();
  last_ccl_deletion_decs_ = last_ccl_cleanup_decs_ = stats.getNumDecisions();
  comp_manager_.initialize(watches_, lit_pool_);

  const auto exit_state = countSAT();
  stats.num_long_red_clauses_ = num_conflict_clauses();
  comp_manager_.gatherStatistics();
  if (config_.verb) stats.printShort(this, &comp_manager_.get_cache());
  if (exit_state == RESTART) {
    largest_cube_ret = largest_cube;
    return largest_cube_val;
  } else {
    assert(exit_state == SUCCESS);
    largest_cube_ret.clear();
    return decision_stack_.top().getTotalModelCount();
  }
}

void Solver::set_target_polar(const vector<CMSat::lbool>& model) {
  assert(target_polar.size() > nVars());
  for(uint32_t i = 0; i < nVars(); i++) {
    target_polar[i+1] = model[i] == CMSat::l_True;
  }
  counted_bottom_comp = false;
}

void Solver::print_all_levels() {
  cout << COLORG "--- going through all decision levels now, printing comps --" << endl;
  uint32_t dec_lev = 0;
  for(const auto& s: decision_stack_) {
    auto const& sup_at = s.super_comp();
    cout << COLORG "super comp of dec_lev " << dec_lev
      << " is at comp_stack_ position: " << sup_at
      << " branch var here: " << decision_stack_.at(dec_lev).getbranchvar()
      << " unproc'd comp end: " << decision_stack_.at(dec_lev).getUnprocessedComponentsEnd()
      << " remaining comp ofs: " << decision_stack_.at(dec_lev).remaining_comps_ofs()
      << " num unproc'd comps: " << decision_stack_.at(dec_lev).numUnprocessedComponents()
      << " count: " << decision_stack_.at(dec_lev).getTotalModelCount()
      << endl;

    const auto& c = comp_manager_.at(sup_at);
    cout << COLORG "-> Variables in comp_manager_.at(" << sup_at << ")."
      << " num vars: " << c->nVars() << " vars: ";
    for(uint32_t i = 0; i < c->nVars(); i++) cout << c->varsBegin()[i] << " ";
    cout << endl;
    dec_lev++;
  }
  cout << COLORG "--- Went through all levels now --" << COLDEF << endl;
}

SOLVER_StateT Solver::countSAT() {
  retStateT state = RESOLVED;

  while (true) {
    print_debug("var top of decision stack: " << decision_stack_.top().getbranchvar());
    // NOTE: findNextRemainingComponentOf finds disjoint comps
    // we then solve them all with the decideLiteral & calling findNext.. again
    while (comp_manager_.findNextRemainingComponentOf(decision_stack_.top())) {
      // It's a component. It will ONLY fall into smaller pieces if we decide on a literal
      checkProbabilisticHashSanity();
      decideLiteral();
      VERBOSE_DEBUG_DO(print_all_levels());

      while (!prop_and_probe()) {
        state = resolveConflict();
        if (state == BACKTRACK) break;
      }
      if (state == BACKTRACK) break;
    }
    // we are here because there is no next component, or we had to backtrack

    if (restart_if_needed()) {return RESTART;}
    state = backtrack();
    if (state == EXIT) return SUCCESS;
    while (state != PROCESS_COMPONENT && !prop_and_probe()) {
      state = resolveConflict();
      if (state == BACKTRACK) {
        state = backtrack();
        if (state == EXIT) return SUCCESS;
      }
    }
  }
  return SUCCESS;
}

bool Solver::get_polarity(const uint32_t v)
{
  bool polarity;
  if (config_.do_restart && decision_stack_.top().on_path_to_target_) polarity = target_polar[v];
  else {
    if (var(Lit(v, false)).set_once) {
      polarity = var(Lit(v, false)).last_polarity;
    } else {
      polarity = false;
      /* polarity = litWatchList(Lit(v, true)).activity_score_ > */
      /*   litWatchList(Lit(v, false)).activity_score_; */
    }
  }
  return polarity;
}

void Solver::decideLiteral() {
  print_debug("new decision level is about to be created, lev now: " << decision_stack_.get_decision_level() << " on path: " << decision_stack_.top().on_path_to_target_ << " branch: " << decision_stack_.top().is_right_branch());
  bool on_path = true;
  if (decision_stack_.size() != 1)
    on_path = decision_stack_.top().on_path_to_target_ && !decision_stack_.top().is_right_branch();
  decision_stack_.push_back(
    StackLevel(decision_stack_.top().currentRemainingComponent(),
               trail.size(),
               comp_manager_.comp_stack_size()));
  decision_stack_.top().on_path_to_target_ = on_path;

  // Find variable to branch on
  auto it = comp_manager_.getSuperComponentOf(decision_stack_.top()).varsBegin();
  uint32_t max_score_var = *it;
  double max_score = scoreOf(*(it));

  // Find one variable that's OK to use
  while (*it != varsSENTINEL && indep_support_.find(*it) == indep_support_.end()) {
    it++;
  }
  if (*it != varsSENTINEL) {
    max_score_var = *it;
    max_score = scoreOf(*it);
  }

  // Find best variable to use
  while (*it != varsSENTINEL) {
    //TODO MATE this is expensive I think
    if (indep_support_.find(*it) != indep_support_.end()) {
      const double score = scoreOf(*it);
      if (score > max_score) {
        max_score = score;
        max_score_var = *it;
      }
    }
    it++;
  }
  assert(max_score_var != 0 &&
        "this assert should always hold, if not then there is a bug in the logic of countSAT()");

  // The decision literal is now ready. Deal with it.
  const Lit lit(max_score_var, get_polarity(max_score_var));
  print_debug(COLYEL "decideLiteral() is deciding: " << lit << " dec level: "
      << decision_stack_.get_decision_level());
  decision_stack_.top().setbranchvariable(max_score_var);
  setLiteralIfFree(lit);
  stats.num_decisions_++;
  if (stats.num_decisions_ % 128 == 0) decayActivities();
  assert( decision_stack_.top().remaining_comps_ofs() <= comp_manager_.comp_stack_size());
}

void Solver::computeLargestCube()
{
  assert(config_.do_restart);
  largest_cube.clear();
  print_debug(COLWHT "-- computeLargestCube BEGIN");
  print_debug_noendl(COLWHT "Decisions in the cube: ");

  // add decisions, components, and counts
  largest_cube_val = decision_stack_.top().getTotalModelCount();
  VERBOSE_DEBUG_DO(bool error = false;);
  for(uint32_t i = 0; i < decision_stack_.size()-1; i++) {
    const StackLevel& dec = decision_stack_[i];
    const Lit dec_lit = trail[dec.trail_ofs()];
    // Add decision
    if (i > 0) {
      const auto dec_lit2 = (target_polar[dec.getbranchvar()] ? 1 : -1)*(int)dec.getbranchvar();
      if (dec_lit2 != dec_lit.toInt()) {
        cout << "(ERROR with dec_lit: " << dec_lit << " dec_lit2: " << dec_lit2 << ") ";
        VERBOSE_DEBUG_DO(error = true;);
      }
      largest_cube.push_back(dec_lit.neg());
      print_debug_noendl(dec_lit.neg() << " ");
    }
    if (dec.getTotalModelCount() > 0) largest_cube_val *= dec.getTotalModelCount();

    const auto off_start = dec.remaining_comps_ofs();
    const auto off_end = dec.getUnprocessedComponentsEnd();
    // add all but the last component (it's the one we just counted)
    for(uint32_t i2 = off_start; i2 < off_end-1; i2++) {
      const auto& c = comp_manager_.at(i2);
      for(auto v = c->varsBegin(); *v != varsSENTINEL; v++)
        largest_cube.push_back(Lit(*v, !target_polar[*v]));
    }
  }
  print_debug_noendl(endl);

#ifdef VERBOSE_DEBUG
  // Show decision stack's comps
  for(size_t i = 0; i < decision_stack_.size(); i++) {
    const auto& dst = decision_stack_.at(i);
    const auto dec_lit = (target_polar[dst.getbranchvar()] ? 1 : -1)*(int)dst.getbranchvar();
    /* const auto dec_lit2 = trail[ds.trail_ofs()]; */
    /* cout << "dec_lit2: " << dec_lit2 << endl; */
    print_debug(COLWHT "decision_stack.at(" << i << "):"
      << " decision lit: " << dec_lit
      << " num unproc comps: " << dst.numUnprocessedComponents()
      << " unproc comps end: " << dst.getUnprocessedComponentsEnd()
      << " remain comps offs: " << dst.remaining_comps_ofs()
      << " count here: " << dst.getTotalModelCount()
      << " on path: " << dst.on_path_to_target_
      << " branch: " << dst.is_right_branch());
    const auto off_start = dst.remaining_comps_ofs();
    const auto off_end = dst.getUnprocessedComponentsEnd();
    for(uint32_t i2 = off_start; i2 < off_end; i2++) {
      assert(i2 < comp_manager_.comp_stack_size());
      const auto& c = comp_manager_.at(i2);
      cout << COLWHT "-> comp at: " << std::setw(3) << i2 << " ID: " << c->id() << " -- vars : ";
      for(auto v = c->varsBegin(); *v != varsSENTINEL; v++) cout << *v << " ";
      cout << COLDEF << endl;
    }
  }

  // All comps
  print_debug(COLWHT "== comp list START");
  for(uint32_t i2 = 0; i2 < comp_manager_.comp_stack_size(); i2++) {
    const auto& c = comp_manager_.at(i2);
    cout << COLWHT "== comp at: " << std::setw(3) << i2 << " ID: " << c->id() << " -- vars : ";
    if (c->empty()) { cout << "EMPTY" << endl; continue; }
    for(auto v = c->varsBegin(); *v != varsSENTINEL; v++) cout << *v << " ";
    cout << endl;
  }
  print_debug(COLWHT "== comp list END");

  cout << COLWHT "Largest cube so far. Size: " << largest_cube.size() << " cube: ";
  for(const auto& l: largest_cube) cout << l << " ";
  cout << endl;
  print_debug(COLWHT "cube's SOLE count: " << decision_stack_.top().getTotalModelCount());
  print_debug(COLWHT "cube's RECORDED count: " << largest_cube_val);
  assert(!error);
#endif
}

bool Solver::restart_if_needed() {
  if (config_.do_restart && stats.getNumDecisions() > stats.next_restart &&
      // don't restart if we are about to exit (i.e. empty largest cube)
      !largest_cube.empty()) {
    return true;
  }
  return false;
}

retStateT Solver::backtrack() {
  assert(decision_stack_.top().remaining_comps_ofs() <= comp_manager_.comp_stack_size());

  do {
    if (decision_stack_.top().branch_found_unsat()) {
      comp_manager_.removeAllCachePollutionsOf(decision_stack_.top());
    } else if (decision_stack_.top().anotherCompProcessible()) {
      print_debug("Processing another comp at dec lev "
          << decision_stack_.get_decision_level()
          << " instead of backtracking." << " Num unprocessed comps: "
          << decision_stack_.top().numUnprocessedComponents()
          << " so far the count: " << decision_stack_.top().getTotalModelCount());
      return PROCESS_COMPONENT;
    }

    // We have NOT explored the other side! Let's do it now!
    if (!decision_stack_.top().is_right_branch()) {
      print_debug("We have NOT explored the right branch (isSecondBranch==false). Let's do it!"
          << " -- dec lev: " << decision_stack_.get_decision_level());
      const Lit aLit = top_dec_lit();
      assert(decision_stack_.get_decision_level() > 0);
      decision_stack_.top().change_to_right_branch();
      reactivate_comps_and_backtrack_trail();
      print_debug("Flipping lit to: " << aLit.neg());
      setLiteralIfFree(aLit.neg(), NOT_A_CLAUSE);
      print_debug(COLORGBG "Backtrack finished -- we flipped the branch");
      return RESOLVED;
    }
    print_debug(COLORGBG "We have explored BOTH branches, actually BACKTRACKING."
        << " -- dec lev: " << decision_stack_.get_decision_level());
    comp_manager_.cacheModelCountOf(decision_stack_.top().super_comp(),
                                    decision_stack_.top().getTotalModelCount());

    // Backtrack from end, i.e. finished.
    if (decision_stack_.get_decision_level() == 0) {
      print_debug("Backtracking from lev 0, i.e. ending");
      break;
    }

    if (decision_stack_.top().on_path_to_target_) {
      if (!counted_bottom_comp) counted_bottom_comp = true;
      if (config_.do_restart && counted_bottom_comp) computeLargestCube();
    }

    reactivate_comps_and_backtrack_trail();
    assert(decision_stack_.size() >= 2);
#ifdef VERBOSE_DEBUG
    const auto parent_count_before = (decision_stack_.end() - 2)->getTotalModelCount();
#endif
    (decision_stack_.end() - 2)->includeSolution(decision_stack_.top().getTotalModelCount());
    print_debug("Backtracking from level " << decision_stack_.get_decision_level()
        << " count here is: " << decision_stack_.top().getTotalModelCount());
    decision_stack_.pop_back();
    auto& dst = decision_stack_.top();
    print_debug("-> Backtracked to level " << decision_stack_.get_decision_level()
        // NOTE: -1 here because we have JUST processed the child
        //     ->> (see below nextUnprocessedComponent() call)
        << " num unprocessed comps here: " << dst.numUnprocessedComponents()-1
        << " current count here: " << dst.getTotalModelCount()
        << " branch: " << dst.is_right_branch()
        << " before including child it was: " <<  parent_count_before
        << " on_path: " << dst.on_path_to_target_);

    // step to the next comp not yet processed
    dst.nextUnprocessedComponent();

    assert(dst.remaining_comps_ofs() < comp_manager_.comp_stack_size() + 1);
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
  act_inc *= 1.0/0.95;

  if (decision_stack_.top().is_right_branch()) {
    // Backtracking since finished with this AND the other branch.
    return BACKTRACK;
  }

  Antecedent ant(NOT_A_CLAUSE);
  // this has to be checked since using implicit BCP
  // and checking literals there not exhaustively
  // we cannot guarantee that uip_clauses_.back().front() == TOS_decLit().neg()
  // this is because we might have checked a literal
  // during implict BCP which has been a failed literal
  // due only to assignments made at lower decision levels
  if (!uip_clauses_.back().empty() && uip_clauses_.back().front() == top_dec_lit().neg()) {
    assert(top_dec_lit().neg() == uip_clauses_.back()[0]);
    var(top_dec_lit().neg()).ante = addUIPConflictClause( uip_clauses_.back());
    ant = var(top_dec_lit()).ante;
  }
  assert(decision_stack_.get_decision_level() > 0);
  assert(decision_stack_.top().branch_found_unsat());

  // we do not have to remove pollutions here,
  // since conflicts only arise directly before
  // remaining comps are stored hence
  assert( decision_stack_.top().remaining_comps_ofs() == comp_manager_.comp_stack_size());

  decision_stack_.top().change_to_right_branch();
  const Lit lit = top_dec_lit();
  reactivate_comps_and_backtrack_trail();
  if (ant == NOT_A_CLAUSE) {
    print_debug("Conflict pushes us to: " << lit<< " and due to failed literal probling, we can't guarantee it's due to the 1UIP, so setting it as a decision instead");
  } else {
    print_debug("Conflict pushes us to: " << lit);
  }
  setLiteralIfFree(lit.neg(), ant);
  return RESOLVED;
}

bool Solver::prop_and_probe() {
  // the asserted literal has been set, so we start
  // bcp on that literal
  assert(trail.size() > 0 && "Mate added this, but it seems OK");

  const uint32_t start_ofs = trail.size() - 1;
  print_debug("--> Setting units of this comp...");
  for (const auto& lit : unit_clauses_) setLiteralIfFree(lit);
  print_debug("--> Units of this comp set, propagating");

  bool bSucceeded = propagate(start_ofs);
  if (config_.do_failed_lit_probe && bSucceeded) {
    bSucceeded = failedLitProbeInternal();
  }
  return bSucceeded;
}

bool Solver::propagate(const uint32_t start_at_trail_ofs) {
  for (auto i = start_at_trail_ofs; i < trail.size(); i++) {
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

void Solver::get_activities(vector<double>& acts, vector<uint8_t>& polars,
    double& ret_act_inc) const
{
  acts.clear();
  for(const auto& v: variables_) {
    acts.push_back(v.activity);
  }
  polars.clear();
  for(const auto& v: variables_) {
    polars.push_back(v.last_polarity);
  }
  ret_act_inc = act_inc;
}

void Solver::shuffle_activities()
{
  for(auto& v: variables_) {
    v.activity += mtrand.randDblExc(1000);
  }
}

void Solver::set_activities(const vector<double>& act, const vector<uint8_t>& polars,
    double ret_act_inc)
{
  size_t i = 0;
  for(auto& v: variables_) {
    v.activity = act[i];
    i++;
  }

  i = 0;
  for(auto& v: variables_) {
    v.set_once = true;
    v.last_polarity = polars[i];
    i++;
  }

  act_inc = ret_act_inc;
}

const DataAndStatistics& Solver::get_stats() const
{
  return stats;
}

bool Solver::failedLitProbeInternal() {
  print_debug(COLRED "Failed literal probing START");

  uint32_t trail_ofs = decision_stack_.top().trail_ofs();
  uint32_t num_curr_lits = 0;
  while (trail_ofs < trail.size()) {
    test_lits.clear();
    for (auto it = trail.begin() + trail_ofs;
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
    num_curr_lits = trail.size() - trail_ofs;
    trail_ofs = trail.size();
    for (auto jt = test_lits.begin(); jt != test_lits.end(); jt++) {
      viewed_lits[*jt] = false;
    }

    // Figure out which literals to probe
    vector<double> scores;
    scores.clear();
    for (auto jt = test_lits.begin(); jt != test_lits.end(); jt++) {
      scores.push_back(variables_[jt->var()].activity);
    }
    sort(scores.begin(), scores.end());
    num_curr_lits = 10 + num_curr_lits / 20;
    double threshold = 0.0;
    if (scores.size() > num_curr_lits) {
      threshold = scores[scores.size() - num_curr_lits];
    }
    stats.num_failed_lit_tests_ += test_lits.size();

    // Do the probing
    for (auto lit : test_lits) {
      if (isUnknown(lit) && threshold <= variables_[lit.var()].activity) {
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
            setLiteralIfFree(it->front(), addUIPConflictClause(*it));
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

void Solver::minimizeAndStoreUIPClause(Lit uipLit, vector<Lit> &cl, const vector<uint8_t>& seen) {
  tmp_clause_minim.clear();
  assertion_level_ = 0;
  for (const auto& lit : cl) {
    if (existsUnitClauseOf(lit.var())) continue;
    bool resolve_out = false;
    if (hasAntecedent(lit)) {
      resolve_out = true;
      if (getAntecedent(lit).isAClause()) {
        for (auto it = beginOf(getAntecedent(lit).asCl()) + 1; *it != SENTINEL_LIT; it++) {
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
  tmp_seen.resize(nVars()+1, false);
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
    increaseActivity(l);
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
           *it != SENTINEL_LIT; it++) {
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
      increaseActivity(alit);
      increaseActivity(curr_lit);
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
  tmp_seen.resize(nVars()+1, false);
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
    increaseActivity(l);
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
           *it != SENTINEL_LIT; it++) {
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
      increaseActivity(alit);
      increaseActivity(curr_lit);
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
    cout << "/" << stats.num_binary_red_clauses_ << "/"
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
