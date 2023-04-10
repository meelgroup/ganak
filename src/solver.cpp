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

void Counter::simplePreProcess()
{
  for (auto lit : unit_clauses_) {
    assert(!isUnitClause(lit.neg()) && "Formula is not UNSAT, we ran CMS before");;
    setLiteralIfFree(lit);
  }


  bool succeeded = propagate(0);
  release_assert(succeeded && "We ran CMS before, so it cannot be UNSAT");
  viewed_vars.resize(nVars() + 1, 0);
  stats.num_unit_irred_clauses_ = unit_clauses_.size();
  irred_lit_pool_size_ = lit_pool_.size();
  init_decision_stack();
}

void Counter::set_indep_support(const set<uint32_t> &indeps)
{
  indep_support_ = indeps;
  indep_support_lookup.clear();
  indep_support_lookup.resize(nVars()+1, 0);
  for(const auto& v: indep_support_) {
    indep_support_lookup[v] = 1;
  }
}

void Counter::init_activity_scores()
{

  for (auto l = Lit(1, false); l != watches_.end_lit(); l.inc())
  {
    watches_[l].activity = watches_[l].binary_links_.size() + occ_lists_[l].size();
  }
}

void Counter::end_irred_cls()
{
  tmp_seen.resize(nVars()+1, 0);
  comp_manager_ = new ComponentManager(config_,stats, lit_values_, indep_support_, this);
  if (config_.do_pcc) comp_manager_->getrandomseedforclhash();
  depth_queue.clearAndResize(config_.first_restart);
  cache_miss_rate_queue.clearAndResize(config_.first_restart);
  comp_size_queue.clearAndResize(config_.first_restart);

  release_assert(!ended_irred_cls && "ERROR *must not* call end_irred_cls() twice");
  stats.maximum_cache_size_bytes_ = config_.maximum_cache_size_bytes_;
  init_decision_stack();
  simplePreProcess();
  ended_irred_cls = true;

  if (config_.verb) stats.printShortFormulaInfo();
  comp_manager_->initialize(watches_, lit_pool_);
}

void Counter::add_red_cl(const vector<Lit>& lits) {
  release_assert(ended_irred_cls && "ERROR *must* call end_irred_cls() before add_red_cl()");
  assert(lits.size() <= 2); // TODO longer clauses -- but then have to add activity
  addUIPConflictClause(lits);
}

void Counter::get_unit_cls(vector<Lit>& units) const
{
  assert(units.empty());
  units = unit_clauses_;
}

mpz_class Counter::count(vector<Lit>& largest_cube_ret)
{
  release_assert(ended_irred_cls && "ERROR *must* call end_irred_cls() before solve()");
  largest_cube.clear();
  largest_cube_val = 0;
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

  const auto exit_state = countSAT();
  stats.num_long_red_clauses_ = num_conflict_clauses();
  comp_manager_->gatherStatistics();
  if (config_.verb) stats.printShort(this, &comp_manager_->get_cache());
  if (exit_state == RESTART) {
    largest_cube_ret = largest_cube;
    return largest_cube_val;
  } else {
    assert(exit_state == SUCCESS);
    largest_cube_ret.clear();
    return decision_stack_.top().getTotalModelCount();
  }
}

void Counter::set_target_polar(const vector<CMSat::lbool>& model) {
  assert(target_polar.size() > nVars());
  for(uint32_t i = 0; i < nVars(); i++) {
    target_polar[i+1] = model[i] == CMSat::l_True;
  }
  counted_bottom_comp = false;
}

void Counter::print_all_levels() {
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

    const auto& c = comp_manager_->at(sup_at);
    cout << COLORG "-> Variables in comp_manager_->at(" << sup_at << ")."
      << " num vars: " << c->nVars() << " vars: ";
    for(uint32_t i = 0; i < c->nVars(); i++) cout << c->varsBegin()[i] << " ";
    cout << endl;
    dec_lev++;
  }
  cout << COLORG "--- Went through all levels now --" << COLDEF << endl;
}

void Counter::print_stat_line() {
  if (next_print_stat > stats.num_cache_hits_) return;
  stats.printShort(this, &comp_manager_->get_cache());
  next_print_stat = stats.num_cache_hits_ + 9000000;
}

SOLVER_StateT Counter::countSAT() {
  retStateT state = RESOLVED;

  while (true) {
    print_debug("var top of decision stack: " << decision_stack_.top().getbranchvar());
    // NOTE: findNextRemainingComponentOf finds disjoint comps
    // we then solve them all with the decideLiteral & calling findNext.. again
    while (comp_manager_->findNextRemainingComponentOf(decision_stack_.top())) {
      // It's a component. It will ONLY fall into smaller pieces if we decide on a literal
      checkProbabilisticHashSanity();
      if (restart_if_needed()) {return RESTART;}
      decideLiteral();
      VERBOSE_DEBUG_DO(print_all_levels());
      print_stat_line();

      while (!prop_and_probe()) {
        state = resolveConflict();
        if (state == BACKTRACK) break;
      }
      if (state == BACKTRACK) break;
    }
    // we are here because there is no next component, or we had to backtrack

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

bool Counter::get_polarity(const uint32_t v)
{
  bool polarity;
  if (config_.do_restart && decision_stack_.top().on_path_to_target_) polarity = target_polar[v];
  else {
    if (var(Lit(v, false)).set_once) {
      polarity = var(Lit(v, false)).last_polarity;
    } else {
      return false;
    }
  }
  return polarity;
}

void Counter::decideLiteral(Lit lit) {
  print_debug("new decision level is about to be created, lev now: " << decision_stack_.get_decision_level() << " on path: " << decision_stack_.top().on_path_to_target_ << " branch: " << decision_stack_.top().is_right_branch());
  bool on_path = true;
  if (decision_stack_.size() != 1)
    on_path = decision_stack_.top().on_path_to_target_ && !decision_stack_.top().is_right_branch();
  decision_stack_.push_back(
    StackLevel(decision_stack_.top().currentRemainingComponent(),
               trail.size(),
               comp_manager_->comp_stack_size()));
  decision_stack_.top().on_path_to_target_ = on_path;

  // The decision literal is now ready. Deal with it.
  if (lit == NOT_A_LIT) {
    const uint32_t v = find_best_branch();
    lit = Lit(v, get_polarity(v));
  }
  print_debug(COLYEL "decideLiteral() is deciding: " << lit << " dec level: "
      << decision_stack_.get_decision_level());
  decision_stack_.top().setbranchvariable(lit.var());
  setLiteralIfFree(lit);
  stats.num_decisions_++;
  if (stats.num_decisions_ % 128 == 0) decayActivities(config_.exp == 1.0);
  assert( decision_stack_.top().remaining_comps_ofs() <= comp_manager_->comp_stack_size());
}

double Counter::alternate_score(uint32_t v, bool val)
{
  double score = 0;

  auto before = decision_stack_.top();
  decision_stack_.top().setbranchvariable(v);
  setLiteralIfFree(Lit(v, val));
  const uint32_t start_ofs = trail.size() - 1;
  bool bSucceeded = propagate(start_ofs);
  if (!bSucceeded) score = 30000;
  else {
    auto& top = decision_stack_.top();
    score = comp_manager_->get_comp_score(top);
  }
  /* uint32_t diff = trail.size() - (start_ofs); */
  reactivate_comps_and_backtrack_trail();
  decision_stack_.pop_back();
  decision_stack_.push_back(before);

  double c = watches_[Lit(v, val)].activity;
  return score*c;
}

uint32_t Counter::find_best_branch()
{
  vars_scores.clear();
  bool lookahead_try = decision_stack_.size() > depth_queue.getLongtTerm().avg()*config_.lookahead_depth;
  uint32_t best_var = 0;
  double best_var_score = -1;
  auto it = comp_manager_->getSuperComponentOf(decision_stack_.top()).varsBegin();
  while (*it != varsSENTINEL) {
    if (indep_support_lookup[*it]) {
      const double score = scoreOf(*it) ;
      if (lookahead_try) vars_scores.push_back(VS(*it, score, *it));
      if (best_var_score == -1 || score > best_var_score) {
        best_var = *it;
        best_var_score = score;
      }
    }
    it++;
  }
  assert(best_var != 0 && best_var_score != -1);

  if (vars_scores.size() > 30 && lookahead_try) {
    std::sort(vars_scores.begin(), vars_scores.end());
    best_var = vars_scores[0].v;
    stats.lookaheads++;
    stats.lookahead_computes++;
    double best_score = alternate_score(vars_scores[0].v, true) *
      alternate_score(vars_scores[0].v, false);
    for(uint32_t i = 1; i < config_.lookahead_num && i < vars_scores.size(); i ++) {
      stats.lookahead_computes++;
      double score = alternate_score(vars_scores[i].v, true) *
        alternate_score(vars_scores[i].v, false);
      if (score > best_score) {
        best_score = score;
        best_var = vars_scores[i].v;
      }
    }
  }

  /* cout << "Decided on var: " << v << " score: " << max_score << endl; */
  return best_var;
}

void Counter::shuffle_activities(MTRand &mtrand2) {
  for(auto& v: watches_) v.activity=act_inc*mtrand2.randExc();
}

void Counter::computeLargestCube()
{
  assert(config_.do_restart);
  largest_cube.clear();
  print_debug(COLWHT "-- computeLargestCube BEGIN");
  print_debug_noendl(COLWHT "Decisions in the cube: ");

  // add decisions, components, and counts
  largest_cube_val = decision_stack_.top().getTotalModelCount();
#ifdef VERBOSE_DEBUG
  bool error = false;
#endif
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
      const auto& c = comp_manager_->at(i2);
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
      assert(i2 < comp_manager_->comp_stack_size());
      const auto& c = comp_manager_->at(i2);
      cout << COLWHT "-> comp at: " << std::setw(3) << i2 << " ID: " << c->id() << " -- vars : ";
      for(auto v = c->varsBegin(); *v != varsSENTINEL; v++) cout << *v << " ";
      cout << COLDEF << endl;
    }
  }

  // All comps
  print_debug(COLWHT "== comp list START");
  for(uint32_t i2 = 0; i2 < comp_manager_->comp_stack_size(); i2++) {
    const auto& c = comp_manager_->at(i2);
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

bool Counter::restart_if_needed() {
  cache_miss_rate_queue.push(stats.cache_miss_rate());
  depth_queue.push(decision_stack_.size());
  /* if (cache_miss_rate_queue.isvalid()) { */
  /*     cout << " Lterm miss avg: " << cache_miss_rate_queue.getLongtTerm().avg() */
  /*     << " Sterm miss avg: " << cache_miss_rate_queue.avg() */
  /*     << endl; */
  /* } */
  /* if (comp_size_queue.isvalid()) { */
  /*     cout << " Lterm comp size avg: " << comp_size_queue.getLongtTerm().avg() */
  /*     << " Sterm comp size avg: " << comp_size_queue.avg() */
  /*     << endl; */
  /* } */
  /* if (depth_queue.isvalid()) { */
  /*     cout << " Lterm dec avg: " << std::setw(5) << depth_queue.getLongtTerm().avg() */
  /*     << " Sterm dec avg: " << std::setw(5) << depth_queue.avg() */
  /*     << endl; */
  /* } */

  if (!config_.do_restart || largest_cube.empty()) return false;
  bool restart = false;
  if (config_.restart_type == 0
      && comp_size_queue.isvalid() && comp_size_queue.avg() < comp_size_queue.getLongtTerm().avg())
    restart = true;
  if (config_.restart_type == 1
      && cache_miss_rate_queue.isvalid() && cache_miss_rate_queue.avg() > cache_miss_rate_queue.getLongtTerm().avg()*0.95)
    restart = true;

  if (config_.restart_type == 2
      && depth_queue.isvalid() && depth_queue.avg() > depth_queue.getLongtTerm().avg()*1.1)
    restart = true;

  if (config_.restart_type == 3 && stats.num_decisions_ > config_.next_restart)
    restart = true;

  if (config_.restart_type == 4 && stats.total_num_cached_comps_ > config_.next_restart)
    restart = true;

  if (restart) {
    cout << "c Restarting. ";
    if (comp_size_queue.isvalid()) {
      cout << " Lterm comp size avg: " << std::setw(5) << comp_size_queue.getLongtTerm().avg()
      << std::setw(5)<< " Sterm comp size avg: " << comp_size_queue.avg();
    }
    if (cache_miss_rate_queue.isvalid()) {
      cout << " Lterm miss avg: " << std::setw(5) << cache_miss_rate_queue.getLongtTerm().avg()
      << " Sterm miss avg: " << std::setw(5) << cache_miss_rate_queue.avg();
    }
    if (depth_queue.isvalid()) {
      cout << " Lterm dec avg: " << std::setw(5) << depth_queue.getLongtTerm().avg()
      << " Sterm dec avg: " << std::setw(5) << depth_queue.avg();
    }
    cout << " Num units: " << unit_clauses_.size();
    cout << " CC: " << conflict_clauses_.size();
    cout << " Num decisions since last restart: "
      << stats.num_decisions_-stats.last_restart_num_decisions
      << endl;

    stats.last_restart_num_decisions = stats.num_decisions_;
    depth_queue.clear();
    cache_miss_rate_queue.clear();
    comp_size_queue.clear();

    while (decision_stack_.size() > 1) {
      if (decision_stack_.top().branch_found_unsat()
          || !decision_stack_.top().on_path_to_target_) {
        comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
      }
      reactivate_comps_and_backtrack_trail();
      decision_stack_.pop_back();
    }
    stats.total_num_cached_comps_ = 0;
    stats.num_decisions_ = 0;

    // experimental for deleting polluted cubes and re-using GANAK
    /* set<uint32_t> vars; */
    /* for(const auto& v: largest_cube) vars.insert(v.var()); */
    /* comp_manager_->delete_comps_with_vars(vars); */
    return true;
  }
  return false;
}

retStateT Counter::backtrack() {
  assert(decision_stack_.top().remaining_comps_ofs() <= comp_manager_->comp_stack_size());

  do {
    if (decision_stack_.top().branch_found_unsat()) {
      comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
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
    comp_manager_->cacheModelCountOf(decision_stack_.top().super_comp(),
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

    assert(dst.remaining_comps_ofs() < comp_manager_->comp_stack_size() + 1);
  } while (true);
  return EXIT;
}

retStateT Counter::resolveConflict() {
  recordLastUIPCauses();
  act_inc *= 1.0/config_.exp;

  if (stats.num_clauses_learned_ - last_ccl_deletion_decs_ > stats.clause_deletion_interval()) {
    deleteConflictClauses();
    last_ccl_deletion_decs_ = stats.num_clauses_learned_;
  }

  if (stats.num_clauses_learned_ - last_ccl_cleanup_decs_ > 100000) {
    compactConflictLiteralPool();
    last_ccl_cleanup_decs_ = stats.num_clauses_learned_;
  }

  stats.num_conflicts_++;
  assert(decision_stack_.top().remaining_comps_ofs() <= comp_manager_->comp_stack_size());
  assert(uip_clauses_.size() == 1);
  if (uip_clauses_.back().empty()) { cout << "c EMPTY CLAUSE FOUND" << endl; }
  decision_stack_.top().mark_branch_unsat();

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
  assert( decision_stack_.top().remaining_comps_ofs() == comp_manager_->comp_stack_size());

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

bool Counter::prop_and_probe() {
  // the asserted literal has been set, so we start
  // bcp on that literal
  assert(trail.size() > 0 && "Mate added this, but it seems OK");

  const uint32_t start_ofs = trail.size() - 1;
  print_debug("--> Setting units of this comp...");
  for (const auto& lit : unit_clauses_) setLiteralIfFree(lit);
  print_debug("--> Units of this comp set, propagating");

  bool bSucceeded = propagate(start_ofs);
  if (config_.failed_lit_probe_type == 2 && bSucceeded &&
      (double)decision_stack_.size() >
        depth_queue.getLongtTerm().avg()*config_.ratio_flitprobe) {
    bSucceeded = failedLitProbeInternal();
  }
  else if (config_.failed_lit_probe_type == 1 && bSucceeded) {
    bSucceeded = failedLitProbeInternal();
  }
  return bSucceeded;
}

bool Counter::propagate(const uint32_t start_at_trail_ofs) {
  for (auto i = start_at_trail_ofs; i < trail.size(); i++) {
    const Lit unLit = trail[i].neg();

    //Propagate bin clauses
    for (const auto& l : litWatchList(unLit).binary_links_) {
      if (isFalse(l)) {
        setConflictState(unLit, l);
        return false;
      }
      setLiteralIfFree(l, Antecedent(unLit));
    }

    //Propagate long clauses
    auto& ws = litWatchList(unLit).watch_list_;
    auto it2 = ws.begin();
    for (auto it = ws.begin(); it != ws.end(); it++) {
      if (isTrue(it->blckLit)) { *it2++ = *it; continue; }

      const auto ofs = it->ofs;
      bool isLitA = (*beginOf(ofs) == unLit);
      auto p_watchLit = beginOf(ofs) + 1 - isLitA;
      auto p_otherLit = beginOf(ofs) + isLitA;
      if (isTrue(*p_otherLit)) { *it2++ = *it; continue; }

      auto itL = beginOf(ofs) + 2;
      while (isFalse(*itL)) itL++;
      // either we found a free or satisfied lit
      if (*itL != SENTINEL_LIT) {
        litWatchList(*itL).addWatchLinkTo(ofs, *p_watchLit);
        std::swap(*itL, *p_watchLit);
      } else {
        // or p_unLit stays resolved
        // and we have hence no free literal left
        // for p_otherLit remain poss: Active or Resolved
        if (setLiteralIfFree(*p_otherLit, Antecedent(ofs))) { // implication
          if (isLitA) std::swap(*p_otherLit, *p_watchLit);
          *it2++ = *it;
        } else {
          setConflictState(ofs);
          while(it != ws.end()) *it2++ = *it++;
          ws.resize(it2-ws.begin());
          return false;
        }
      }
    }
    ws.resize(it2-ws.begin());
  }
  return true;
}

void Counter::get_activities(vector<double>& acts, vector<uint8_t>& polars,
    double& ret_act_inc, vector<uint32_t>& comp_acts) const
{

  acts.resize((nVars()+1)*2);
  for (auto l = Lit(1, false); l != watches_.end_lit(); l.inc())
    acts[l.raw()] = watches_[l].activity;
  polars.clear();
  for(const auto& v: variables_) polars.push_back(v.last_polarity);
  comp_acts.clear();
  for(uint32_t i = 0; i < nVars()+1; i++) comp_acts.push_back(comp_manager_->scoreOf(i));
  ret_act_inc = act_inc;
}

void Counter::set_activities(const vector<double>& acts, const vector<uint8_t>& polars,
    double ret_act_inc, vector<uint32_t>& comp_acts)
{
  for (auto l = Lit(1, false); l != watches_.end_lit(); l.inc())
    watches_[l].activity = acts[l.raw()];

  for(uint32_t i = 0; i < nVars()+1; i++) comp_manager_->scoreOf(i) = comp_acts[i];

  uint32_t i = 0;
  for(auto& v: variables_) {
    v.set_once = true;
    v.last_polarity = polars[i];
    i++;
  }

  act_inc = ret_act_inc;
}

const DataAndStatistics& Counter::get_stats() const
{
  return stats;
}

bool Counter::failedLitProbeInternal() {
  print_debug(COLRED "Failed literal probing START");

  uint32_t trail_ofs = decision_stack_.top().trail_ofs();
  uint32_t num_curr_lits = 0;
  while (trail_ofs < trail.size()) {
    test_lits.clear();
    for (auto it = trail.begin() + trail_ofs; it != trail.end(); it++) {
      // Only going through the long, the binary clauses have set the variables already
      for (auto cl_ofs : occ_lists_[it->neg()]) {
        if (!isSatisfied(cl_ofs)) {
          for (auto lt = beginOf(cl_ofs); *lt != SENTINEL_LIT; lt++) {
            if (isUnknown(*lt) && !viewed_vars[lt->var()]) {
              test_lits.push_back(*lt);
              print_debug("-> potential lit to test: " << lt->neg());
              viewed_vars[lt->var()] = true;
            }
          }
        }
      }
    }
    num_curr_lits = trail.size() - trail_ofs;
    trail_ofs = trail.size();
    for (const auto& l: test_lits) viewed_vars[l.var()] = false;

    // Figure out which literals to probe
    scores.clear();
    for (const auto& l: test_lits) {
      scores.push_back(std::max(watches_[l].activity, watches_[l.neg()].activity));
    }
    sort(scores.begin(), scores.end());
    num_curr_lits = 5 + num_curr_lits / 40;
    if (!config_.bprop) num_curr_lits*=2;
    double threshold = 0.0;
    if (scores.size() > num_curr_lits) {
      threshold = scores[scores.size() - num_curr_lits];
    }

    // Do the probing
    // TODO: bprop is buggy!! We undercount sometimes.
    toSet.clear();
    for (auto& l : test_lits) {
      if (isUnknown(l) && threshold <=
          std::max(watches_[l].activity, watches_[l.neg()].activity)) {
        if (watches_[l].activity < watches_[l.neg()].activity) l = l.neg();
        if (!one_lit_probe(l, true)) return false;
        if (config_.bprop && isUnknown(l) && !one_lit_probe(l.neg(), false)) return false;
      }
    }
    if (config_.bprop) {
      for(const auto& l: toSet) {
        if (isUnknown(l)) {
          auto sz = trail.size();
          setLiteralIfFree(l, Antecedent(NOT_A_CLAUSE), true);
          bool bSucceeded = propagate(sz);
          if (!bSucceeded) return false;
          stats.num_failed_bprop_literals_failed++;
        }
      }
    }
  }
  print_debug(COLRED "Failed literal probing END -- no UNSAT, gotta check this branch");
  return true;
}

bool Counter::one_lit_probe(Lit lit, bool set)
{
  stats.num_failed_lit_tests_++;
  uint32_t sz = trail.size();
  // we increase the decLev artificially
  // s.t. after the tentative BCP call, we can learn a conflict clause
  // relative to the assignment of *jt
  decision_stack_.startFailedLitTest();
  setLiteralIfFree(lit);
  assert(!hasAntecedent(lit));
  if (set == true) assert(toClear.empty());

  bool bSucceeded = propagate(sz);
  if (!bSucceeded) {
    for(const auto& v: toClear) tmp_seen[v] = 0;
    toClear.clear();
    recordAllUIPCauses();
  }
  decision_stack_.stopFailedLitTest();

  // backtracking
  while (trail.size() > sz) {
    Lit l = trail.back();
    unSet(l);
    if (config_.bprop && bSucceeded) {
      if (set) {
        tmp_seen[l.var()] = 1U | ((uint8_t)l.sign() << 1);
        toClear.push_back(l.var());
      } else {
        if (tmp_seen[l.var()] && (tmp_seen[l.var()] >> 1) == (uint8_t)l.sign()) {
          toSet.insert(l);
        }
      }
    }
    trail.pop_back();
  }
  if (!set) {
    for(const auto& v: toClear) tmp_seen[v] = 0;
    toClear.clear();
  }

  if (!bSucceeded) {
    stats.num_failed_literals_detected_++;
    print_debug("-> failed literal detected");
    sz = trail.size();
    for (auto it = uip_clauses_.rbegin(); it != uip_clauses_.rend(); it++) {
      if (it->size() == 0) cout << "c EMPTY CLAUSE FOUND" << endl;
      setLiteralIfFree(it->front(), addUIPConflictClause(*it));
    }
    if (!propagate(sz)) {
      print_debug("Failed literal probing END -- this comp/branch is UNSAT");
      return false;
    }
  }

  return true;
}

void Counter::minimizeAndStoreUIPClause(Lit uipLit, vector<Lit> &cl) {
  tmp_clause_minim.clear();
  assertion_level_ = 0;
  for (const auto& lit : cl) {
    if (existsUnitClauseOf(lit.var())) continue;
    bool resolve_out = false;
    if (hasAntecedent(lit)) {
      resolve_out = true;
      if (antedecentBProp(lit)) {
        // BProp is the reason
        for (int32_t i = 1; i < variables_[lit.var()].decision_level+1; i++) {
          const Lit l = trail[decision_stack_[i].trail_ofs()];
          if (!tmp_seen[l.var()]) {
            resolve_out = false;
            break;
          }
        }
      } else if (getAntecedent(lit).isAClause()) {
        // Long clause reason
        for (auto it = beginOf(getAntecedent(lit).asCl()) + 1; *it != SENTINEL_LIT; it++) {
          if (!tmp_seen[it->var()]) {
            resolve_out = false;
            break;
          }
        }
      } else if (!tmp_seen[getAntecedent(lit).asLit().var()]) {
        // Binary clause reason
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

void Counter::recordLastUIPCauses() {
  // note:
  // variables of lower dl: if seen we dont work with them anymore
  // variables of this dl: if seen we incorporate their
  // antecedent and set to unseen
  tmp_clause.clear();
  assert(toClear.empty());

  assertion_level_ = 0;
  uip_clauses_.clear();

  uint32_t trail_ofs = trail.size();
  const uint32_t DL = decision_stack_.get_decision_level();
  uint32_t lits_at_current_dl = 0;

  for (const auto& l: violated_clause) {
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
      // perform UIP stuff
      if (!hasAntecedent(curr_lit)) {
        // this should be the decision literal when in first branch
        // or it is a literal decided to explore in failed literal testing
        break;
      }
    }

    assert(hasAntecedent(curr_lit));
    if (antedecentBProp(curr_lit)) {
      // BProp is the reason
      for (int32_t i = 1; i < variables_[curr_lit.var()].decision_level+1; i++) {
        const Lit l = trail[decision_stack_[i].trail_ofs()];
        if (tmp_seen[l.var()] || (var(l).decision_level == 0) ||
            existsUnitClauseOf(l.var())) {
          continue;
        }
        if (var(l).decision_level < (int)DL) {
          tmp_clause.push_back(l);
        } else {
          lits_at_current_dl++;
        }
        tmp_seen[l.var()] = true;
        toClear.push_back(l.var());
      }
    } else if (getAntecedent(curr_lit).isAClause()) {
      // Long clause is the reason
      updateActivities(getAntecedent(curr_lit).asCl());
      assert(curr_lit == *beginOf(getAntecedent(curr_lit).asCl()));

      for (auto it = beginOf(getAntecedent(curr_lit).asCl()) + 1; *it != SENTINEL_LIT; it++) {
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
      // Binary clause is reason
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
  minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause);
  for(const auto& v: toClear) tmp_seen[v] = false;
  toClear.clear();
}

void Counter::recordAllUIPCauses() {
  // note:
  // variables of lower dl: if seen we dont work with them anymore
  // variables of this dl: if seen we incorporate their
  // antecedent and set to unseen
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
    if (config_.alluip_inc_act) increaseActivity(l);
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
      minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause);
    }

    assert(hasAntecedent(curr_lit));
    if (antedecentBProp(curr_lit)) {
      // BProp is the reason
      for (int32_t i = 1; i < variables_[curr_lit.var()].decision_level+1; i++) {
        const Lit l = trail[decision_stack_[i].trail_ofs()];
        if (tmp_seen[l.var()] || (var(l).decision_level == 0) ||
            existsUnitClauseOf(l.var())) {
          continue;
        }
        if (var(l).decision_level < (int)DL) {
          tmp_clause.push_back(l);
        } else {
          lits_at_current_dl++;
        }
        tmp_seen[l.var()] = true;
        toClear.push_back(l.var());
      }
    } else if (getAntecedent(curr_lit).isAClause()) {
      // Long clause is the reason
      if (config_.alluip_inc_act) updateActivities(getAntecedent(curr_lit).asCl());
      assert(curr_lit == *beginOf(getAntecedent(curr_lit).asCl()));

      for (auto it = beginOf(getAntecedent(curr_lit).asCl()) + 1; *it != SENTINEL_LIT; it++) {
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
      // Binary clause is the reason
      Lit alit = getAntecedent(curr_lit).asLit();
      if (config_.alluip_inc_act) {
        increaseActivity(alit);
        increaseActivity(curr_lit);
      }
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
    minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause);
  }
  for(const auto& v: toClear) tmp_seen[v] = false;
  toClear.clear();
}

void Counter::printOnlineStats() {
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

    comp_manager_->gatherStatistics();

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

Counter::Counter(const CounterConfiguration& conf)
{
  config_ = conf;
  mtrand.seed(conf.seed);
}

Counter::~Counter()
{
  delete comp_manager_;
}
