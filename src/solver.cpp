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

#include "solver.h"

#include <algorithm>
#include <complex>
#include <ios>
#include <iomanip>
#include <numeric>
#include "common.h"
#include "comp_types/comp.h"
#include "cryptominisat5/cryptominisat.h"
#include "cryptominisat5/solvertypesmini.h"
#include "primitive_types.h"
#include "stack.h"
#include "structures.h"
#include "time_mem.h"
#include "IFlowCutter.h"
#include "graph.hpp"

void Counter::simplePreProcess()
{
  for (auto lit : unit_clauses_) {
    assert(!isUnitClause(lit.neg()) && "Formula is not UNSAT, we ran CMS before");;
    setLiteralIfFree(lit);
  }

  bool succeeded = propagate(0);
  release_assert(succeeded && "We ran CMS before, so it cannot be UNSAT");
  viewed_lits.resize(2*(nVars() + 1), 0);
  stats.num_unit_irred_clauses_ = unit_clauses_.size();
  irred_lit_pool_size_ = lit_pool_.size();
  init_decision_stack();
}

void Counter::set_indep_support(const set<uint32_t> &indeps)
{
  if (indeps.count(0)) {
    cout << "ERROR: variable 0 does NOT exist!!" << endl;
    exit(-1);
  }
  vector<uint32_t> tmp(indeps.begin(), indeps.end());
  std::sort(tmp.begin(), tmp.end());
  for(uint32_t i = 0; i < tmp.size(); i++) {
    if (tmp[i] > nVars()) {
      cout << "ERROR: sampling set contains a variable larger than nVars()" << endl;
      exit(-1);
    }
    if (tmp[i] != i+1) {
      cout << "ERROR: independent support MUST start from variable 1 and be consecutive, e.g. 1,2,3,4,5. It cannot skip any variables. You skipped variable: " << i << endl;
      exit(-1);
    }
  }
  if (tmp.size() == 0) indep_support_end = 0;
  else indep_support_end = tmp.back()+1;
  if (indep_support_end == nVars()+1) perform_projected_counting = false;
  else perform_projected_counting = true;
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
  tmp_seen.resize(2*(nVars()+2), 0);
  comp_manager_ = new ComponentManager(config_,stats, lit_values_, indep_support_end, this);
  comp_manager_->getrandomseedforclhash();
  depth_q.clearAndResize(config_.first_restart);
  cache_miss_rate_q.clearAndResize(config_.first_restart);
  comp_size_q.clearAndResize(config_.first_restart);

  release_assert(!ended_irred_cls && "ERROR *must not* call end_irred_cls() twice");
  stats.maximum_cache_size_bytes_ = config_.maximum_cache_size_MB*1024*1024;
  init_decision_stack();
  simplePreProcess();
  ended_irred_cls = true;

  if (config_.verb) stats.printShortFormulaInfo();
  comp_manager_->initialize(watches_, lit_pool_, nVars());
}

void Counter::add_red_cl(const vector<Lit>& lits, int lbd) {
  assert(ended_irred_cls);
  for(const auto& l: lits) release_assert(l.var() <= nVars());
  for(const auto& l: lits) release_assert(isUnknown(l));
  assert(lits.size() >= 2 && "No unit or empty clauses please");

  // NOTE: since we called end_irred_cls, this binary will NOT end up
  //       through ComponentAnalyzer::initialize in analyzer's
  //       unified_variable_links_lists_pool_ which means it will NOT
  //       connect components -- which is what we want
  ClauseOfs cl_ofs = addClause(lits, true);
  if (cl_ofs != 0) {
    red_cls.push_back(cl_ofs);
    auto& header = getHeaderOf(cl_ofs);
    if (lbd == -1) lbd = lits.size();
    header = ClHeader(lbd, true);
  }
}

void Counter::compute_score(TreeDecomposition& tdec) {
  int weight = -1;
  const int n = nVars();
  assert(tdscore.empty());
  tdscore.resize(nVars()+1);
  if (n <= 2) return;
  const auto& bags = tdec.Bags();
  const auto& adj = tdec.get_adj_list();
#if 0
  for(uint32_t i = 0; i < bags.size(); i++) {
    const auto& b = bags[i];
    cout << "bag id:" << i << endl;
    for(const auto& bb: b) {
      cout << bb << " ";
    }
    cout << endl;
  }
  for(uint32_t i = 0; i < adj.size(); i++) {
    const auto& a = adj[i];
    for(const auto& nn: a) cout << i << " " << nn << endl;
  }
#endif
  sspp::TreeDecomposition dec(bags.size(),nVars()+1);
  for(uint32_t i = 0; i < bags.size();i++) {
    dec.SetBag(i+1, bags[i]);
  }
  for(uint32_t i = 0; i < adj.size(); i++) {
    const auto& a = adj[i];
    for(const auto& nn: a) {
      dec.AddEdge(i+1, nn+1);
    }
  }

  auto width = dec.Width();
  auto ord = dec.GetOrd();
  // We use 1-indexing, ignore index 0

  assert(ord.size() == tdscore.size());
  int max_ord = 0;
  for (int i = 1; i <= n; i++) {
    /* assert(ord[i] >= 1); */
    max_ord = std::max(max_ord, ord[i]);
  }
  assert(max_ord >= 1);
  // Normalize
  for (int i = 1; i <= n; i++) {
    tdscore[i] = max_ord - ord[i];
    tdscore[i] /= (double)max_ord;
    assert(tdscore[i] > -0.01 && tdscore[i] < 1.01);
  }
  // Now scores are between 0..1
  double coef = 1;
  if (weight > 0) {
    double rt = (double)n/(double)width;
    /* cout << "RT" << rt << endl; */
    if (rt > 40) {
      coef = 1e7;
    } else {
      coef = weight*exp(rt)/(double)n;
    }
  } else {
    coef = 1e7;
  }
  coef = std::min(coef, 1e7);
  /* cout << "c o COEF: " << coef << " Width: " << width << endl; */
  for (int i = 1; i <= n; i++) {
    tdscore[i] *= coef;
  }
#ifdef VERBOSE_DEBUG
  for(int i = 1; i <= n; i++) {
    cout << "TD var: " << i << " tdscore: " << tdscore[i] << endl;
  }
#endif
}

void Counter::get_unit_cls(vector<Lit>& units) const
{
  assert(units.empty());
  units = unit_clauses_;
}

void Counter::td_decompose()
{
  bool conditionOnCNF = indep_support_end > 3 && nVars() > 20 && nVars() <= config_.td_varlim;
  if (!conditionOnCNF) {
    verb_print(1, "skipping TD, too many/few vars. Setting branch to fallback");
    config_.branch_type = config_.branch_fallback_type;
    return;
  }

  Graph primal(nVars()+1);
  for(uint32_t i = 2; i < (nVars()+1)*2; i++) {
    Lit l(i/2, i%2);
    for(const auto& l2: watches_[l].binary_links_) {
      if (l < l2) {
        print_debug("v1: " << l.var());
        print_debug("v2: " << l2.var());
        primal.addEdge(l.var(), l2.var());
      }
    }
  }

  for(uint32_t i = ClHeader::overheadInLits()+1; i < irred_lit_pool_size_;
      i+=ClHeader::overheadInLits()) {
    for(; lit_pool_[i] != SENTINEL_LIT; i++) {
      for(uint32_t i2 = i+1; lit_pool_[i2] != SENTINEL_LIT; i2++) {
        print_debug("v1: " << lit_pool_[i].var());
        print_debug("v2: " << lit_pool_[i2].var());
        primal.addEdge(lit_pool_[i].var(), lit_pool_[i2].var());
      }
    }
    i++;
  }
  verb_print(1, "Primal graph: nodes: " << nVars()+1 << ", edges " <<  primal.numEdges());

  double density = (double)primal.numEdges()/(double)(nVars() * nVars());
  double edge_var_ratio = (double)primal.numEdges()/(double)nVars();
  verb_print(1, "Primal graph density: "
    << std::fixed << std::setw(9) << std::setprecision(3) << density
    << " edge/var: "
    << std::fixed << std::setw(9) << std::setprecision(3) << edge_var_ratio);
  /* bool conditionOnPrimalGraph = */
  /*     density <= config_.td_denselim && */
  /*     edge_var_ratio <= config_.td_ratiolim; */

  /* if (!conditionOnPrimalGraph) { */
  /*   verb_print(1, "skipping td, primal graph is too large or dense." */
  /*       " Setting branch to fallback"); */
  /*   config_.branch_type = config_.branch_fallback_type; */
  /*   return; */
  /* } */

  // run FlowCutter
  verb_print(1, "FlowCutter is running...");
  IFlowCutter FC(nVars()+1, primal.numEdges(), 0); //TODO: fix time limit
  FC.importGraph(primal);
  TreeDecomposition td = FC.constructTD();

  td.centroid(nVars()+1);
  compute_score(td);
}

mpz_class Counter::count(vector<Lit>& largest_cube_ret, CMSat::SATSolver* _sat_solver)
{
  sat_solver = _sat_solver;
  release_assert(ended_irred_cls && "ERROR *must* call end_irred_cls() before solve()");
  if (indep_support_end == std::numeric_limits<uint32_t>::max()) indep_support_end = nVars()+2;
  largest_cube.clear();
  largest_cube_val = 0;
  start_time = cpuTime();
  if (config_.verb) { cout << "c Sampling set size: " << indep_support_end-1 << endl; }

  if (config_.branch_type == branch_t::sharptd ||
      config_.branch_type == branch_t::gpmc) td_decompose();

  verb_print(1, "branch type: " << config_.get_branch_type_str());

  const auto exit_state = countSAT();
  stats.num_long_red_clauses_ = red_cls.size();
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
  if (next_print_stat_cache > stats.num_cache_look_ups_ &&
      next_print_stat_confl > stats.conflicts) return;
  if (config_.verb) {
    cout << "c total time: " << (cpuTime() - start_time) << endl;
    stats.printShort(this, &comp_manager_->get_cache());
  }
  next_print_stat_cache = stats.num_cache_look_ups_ + (2ULL*1000LL*1000LL);
  next_print_stat_confl = stats.conflicts + (15LL*1000LL);
}

SOLVER_StateT Counter::countSAT() {
  retStateT state = RESOLVED;

  while (true) {
    print_debug("var top of decision stack: " << decision_stack_.top().getbranchvar());
    // NOTE: findNextRemainingComponentOf finds disjoint comps
    // we then solve them all with the decideLiteral & calling findNext.. again
    while (comp_manager_->findNextRemainingComponentOf(decision_stack_.top())) {
      // It's a component. It will ONLY fall into smaller pieces if we decide on a literal
      /* checkProbabilisticHashSanity(); -- no need, there is no way we get to 2**45 lookups*/
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

bool Counter::standard_polarity(const uint32_t v) const {
    return watches_[Lit(v, true)].activity >
      watches_[Lit(v, false)].activity;
}

bool Counter::get_polarity(const uint32_t v) const
{
  bool polarity;
  if (config_.do_restart && decision_stack_.top().on_path_to_target_) polarity = target_polar[v];
  else {
    if (config_.polar_type == 0) {
      if (var(Lit(v, false)).set_once) {
        polarity = var(Lit(v, false)).last_polarity;
        // TODO ** ONLY ** do it in case it's non-exact, right??
        /* if (config_.do_restart) polarity = !polarity; */
      } else polarity = standard_polarity(v);
    } else if (config_.polar_type == 1) polarity = standard_polarity(v);
    else if (config_.polar_type == 2) polarity = false;
    else if (config_.polar_type == 3) polarity = true;
    else assert(false);
  }
  return polarity;
}

void Counter::decideLiteral() {
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
  uint32_t v;
  isindependent = true;
  if (config_.branch_type == branch_t::gpmc) v = find_best_branch_gpmc(true);
  else v = find_best_branch(true);
  if (v == 0 && perform_projected_counting) {
    isindependent = false;
    if (config_.branch_type == branch_t::gpmc) v = find_best_branch_gpmc(false);
    else v = find_best_branch(false);
  }
  assert(v != 0);
  Lit lit = Lit(v, get_polarity(v));
  print_debug(COLYEL "decideLiteral() is deciding: " << lit << " dec level: "
      << decision_stack_.get_decision_level());
  decision_stack_.top().setbranchvariable(lit.var());
  setLiteralIfFree(lit);
  stats.decisions++;
  if (stats.decisions % 128 == 0) {
    comp_manager_->rescale_cache_scores();
  }
  assert( decision_stack_.top().remaining_comps_ofs() <= comp_manager_->comp_stack_size());
}

double Counter::alternate_score(uint32_t v, bool val)
{
  assert(!perform_projected_counting && "Cannot work with projected model counting");
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

uint32_t Counter::find_best_branch_gpmc(bool do_indep)
{
  uint32_t maxv = 0;
  double max_score_a = -1;
  double max_score_f = -1;
  double max_score_td = -1;

  for (auto it = comp_manager_->getSuperComponentOf(decision_stack_.top()).varsBegin();
      *it != varsSENTINEL; it++) if (!do_indep || *it < indep_support_end) {
    uint32_t v = *it;
    double score_td = tdscore[v];
    double score_f = comp_manager_->scoreOf(v);
    double score_a = watches_[Lit(v, false)].activity + watches_[Lit(v, true)].activity;

    if(score_td > max_score_td) {
      max_score_td = score_td;
      max_score_f = score_f;
      max_score_a = score_a;
      maxv = v;
    }
    else if( score_td == max_score_td) {
      if(score_f > max_score_f) {
        max_score_f = score_f;
        max_score_a = score_a;
        maxv = v;
      } else if (score_f == max_score_f && score_a > max_score_a) {
        max_score_a = score_a;
        maxv = v;
      }
    }
  }
  return maxv;
}

uint32_t Counter::find_best_branch(bool do_indep)
{
  assert(!(config_.do_cache_score && config_.do_lookahead) && "can't have both active");

  vars_scores.clear();
  bool lookahead_try = !config_.do_cache_score && config_.do_lookahead &&
    decision_stack_.size() > depth_q.getLongtTerm().avg()*config_.lookahead_depth;
  uint32_t best_var = 0;
  double best_var_score = -1;
  for (auto it = comp_manager_->getSuperComponentOf(decision_stack_.top()).varsBegin();
      *it != varsSENTINEL; it++) {
    if (!do_indep || *it < indep_support_end) {
      const double score = scoreOf(*it) ;
      if (lookahead_try) vars_scores.push_back(VS(*it, score, *it));
      if (best_var_score == -1 || score > best_var_score) {
        best_var = *it;
        best_var_score = score;
      }
    }
  }

  if (!config_.do_lookahead && config_.do_cache_score && best_var != 0) {
    double cachescore = comp_manager_->cacheScoreOf(best_var);
    for (auto it = comp_manager_->getSuperComponentOf(decision_stack_.top()).varsBegin();
         *it != varsSENTINEL; it++) {
      if (!do_indep || *it < indep_support_end) {
        const double score = scoreOf(*it);
        if (score > best_var_score * 0.9) {
          if (comp_manager_->cacheScoreOf(*it) > cachescore) {
            best_var = *it;
            cachescore = comp_manager_->cacheScoreOf(*it);
          }
        }
      }
    }
  }

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

void Counter::print_restart_data() const
{
  if (comp_size_q.isvalid()) {
    verb_print(1, std::setw(30) << std::left
       << "c Lterm comp size avg: " << std::setw(9) << comp_size_q.getLongtTerm().avg()
       << std::right  << std::setw(30) << std::left
       << std::left   << " Sterm comp size avg: " << comp_size_q.avg());
  }
  if (cache_miss_rate_q.isvalid()) {
    verb_print(1, std::setw(30) << std::left
      << "c Lterm miss avg: " << std::setw(9) << cache_miss_rate_q.getLongtTerm().avg()
      << std::right  << std::setw(30) << std::left
      << std::left   << " Sterm miss avg: " << std::setw(9) << cache_miss_rate_q.avg());
  }
  if (depth_q.isvalid()) {
    verb_print(1, std::setw(30) << std::left
      << "c Lterm dec avg: " << std::setw(9) << depth_q.getLongtTerm().avg()
      << std::right << std::setw(30) << std::left
      << std::left  << " Sterm dec avg: " << std::setw(9) << depth_q.avg());
  }
  if (stats.cache_hits_misses_q.isvalid()) {
    verb_print(1, std::setw(30) << std::left
      << "c Lterm hit avg: " << std::setw(9) << stats.cache_hits_misses_q.getLongtTerm().avg()
      << std::right  << std::setw(30) << std::left
      << std::left   << " Sterm hit avg: " << std::setw(5) << stats.cache_hits_misses_q.avg());
  }
  if (stats.comp_size_times_depth_q.isvalid()) {
    verb_print(1, std::setw(30) << std::left
      << "c Lterm compsz/depth avg: " << std::setw(9)
      << stats.comp_size_times_depth_q.getLongtTerm().avg()
      << std::right  << std::setw(30) << std::left
      << std::left << " Sterm compsz/depth avg: " << std::setw(9) << stats.comp_size_times_depth_q.avg()
      << " depth: " << decision_stack_.size()-1);
  }
  cout << std::right;
}

bool Counter::restart_if_needed() {
  cache_miss_rate_q.push(stats.cache_miss_rate());
  depth_q.push(decision_stack_.size());
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
  /* if (stats.comp_size_per_depth.isvalid()) { */
  /*     cout */
  /*       << " Lterm compsz/depth avg: " */
  /*       << std::setw(9) << stats.comp_size_per_depth.getLongtTerm().avg() */
  /*     << " Sterm compsz/depth avg: " */
  /*     << std::setw(9) << stats.comp_size_per_depth.avg() */
  /*     << " depth: " << decision_stack_.size()-1 */
  /*     << endl; */
  /* } */

  if (!config_.do_restart || largest_cube.empty()) return false;
  bool restart = false;
  if (config_.restart_type == 0
      && comp_size_q.isvalid() &&
      comp_size_q.avg() < comp_size_q.getLongtTerm().avg()*config_.restart_cutoff_mult)
    restart = true;
  if (config_.restart_type == 1
      && cache_miss_rate_q.isvalid() &&
      cache_miss_rate_q.avg() > cache_miss_rate_q.getLongtTerm().avg()*0.95)
    restart = true;

  if (config_.restart_type == 2
      && depth_q.isvalid() &&
      depth_q.avg() > depth_q.getLongtTerm().avg()*(1.0/config_.restart_cutoff_mult))
    restart = true;

  if (config_.restart_type == 3 &&
      (stats.decisions-stats.last_restart_num_decisions) > config_.next_restart)
    restart = true;

  if (config_.restart_type == 4 && stats.cache_hits_misses_q.isvalid()
      && stats.cache_hits_misses_q.avg() <
      stats.cache_hits_misses_q.getLongtTerm().avg()*config_.restart_cutoff_mult)
      restart = true;

  if (config_.restart_type == 5 && stats.comp_size_times_depth_q.isvalid() &&
        stats.comp_size_times_depth_q.avg() >
          stats.comp_size_times_depth_q.getLongtTerm().avg()*(1.0/config_.restart_cutoff_mult))
      restart = true;

  // don't restart if we didn't change the scores
  if (stats.last_restart_num_conflicts == stats.conflicts)
    restart = false;

  if (restart) {
    cout << "c  ************* Restarting.  **************" << endl;
    print_restart_data();
    /* cout << "c Num units: " << unit_clauses_.size(); */
    /* cout << " CC: " << conflict_clauses_.size(); */
    cout << "c Num decisions since last restart: "
      << stats.decisions-stats.last_restart_num_decisions
      << endl;
    cout << "c Num cache lookups since last restart: "
      << stats.num_cache_look_ups_-stats.last_restart_num_cache_look_ups
      << endl;

    depth_q.clear();
    cache_miss_rate_q.clear();
    comp_size_q.clear();
    stats.cache_hits_misses_q.clear();
    stats.comp_size_times_depth_q.clear();

    while (decision_stack_.size() > 1) {
      bool on_path = true;
      if (decision_stack_.top().branch_found_unsat()
          || !decision_stack_.top().on_path_to_target_) {
        comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
        on_path = false;
      }
      if (config_.do_on_path_print) cout << "ON PATH: " << on_path << " -- ";
      reactivate_comps_and_backtrack_trail(config_.do_on_path_print);
      decision_stack_.pop_back();
    }
    stats.last_restart_num_conflicts = stats.conflicts;
    stats.last_restart_num_decisions = stats.decisions;
    stats.last_restart_num_cache_look_ups = stats.num_cache_look_ups_;

    // experimental for deleting polluted cubes and re-using GANAK
    /* set<uint32_t> vars; */
    /* for(const auto& v: largest_cube) vars.insert(v.var()); */
    /* comp_manager_->delete_comps_with_vars(vars); */
    return true;
  }
  return false;
}

retStateT Counter::backtrack_nonindep() {
  assert(!isindependent && perform_projected_counting);
  print_debug("[nonindep] Backtrack nonindep");

  do {
    if (decision_stack_.top().branch_found_unsat()) {
      comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
    } else if (decision_stack_.top().anotherCompProcessible()) {
      print_debug("[nonindep] Processing another comp at dec lev "
          << decision_stack_.get_decision_level()
          << " instead of backtracking." << " Num unprocessed comps: "
          << decision_stack_.top().numUnprocessedComponents()
          << " so far the count: " << decision_stack_.top().getTotalModelCount());
      return PROCESS_COMPONENT;
    }

    // Processed all components, and there is at least 1 solution.
    // Now go back until indep var
    // NOTE: none of the things above can be UNSAT (otherwise we would have exited already)
    //       but they may have unprocessed components (and hence may become UNSAT)
    if (decision_stack_.top().getBranchSols() != 0 && !isindependent) {
      while (decision_stack_.top().getbranchvar() >= indep_support_end) {
        print_debug("[nonindep] Going BACK because it's not independent and there is at least 1 solution");
        assert(!isindependent);
        if (decision_stack_.get_decision_level() <= 0) { break; }
        reactivate_comps_and_backtrack_trail();
        assert(decision_stack_.size() >= 2);
        assert(decision_stack_.top().getTotalModelCount() > 0);
        (decision_stack_.end() - 2)->includeSolution(
            decision_stack_.top().getTotalModelCount() > 0);
        decision_stack_.pop_back();
        isindependent = (decision_stack_.top().getbranchvar() < indep_support_end);
        // step to the next component not yet processed
        decision_stack_.top().nextUnprocessedComponent();
        assert( decision_stack_.top().remaining_comps_ofs() <
            comp_manager_->comp_stack_size() + 1);
        if (decision_stack_.top().anotherCompProcessible()) {
          print_debug("[nonindep] Processing another comp at dec lev "
              << decision_stack_.get_decision_level()
              << " instead of backtracking." << " Num unprocessed comps: "
              << decision_stack_.top().numUnprocessedComponents()
              << " so far the count: " << decision_stack_.top().getTotalModelCount());
          return PROCESS_COMPONENT;
        }
      }
    }

    // Maybe it's zero on this side, let's try the other side.
    if (!decision_stack_.top().is_right_branch() &&
        (isindependent || decision_stack_.top().getTotalModelCount() == 0)) {
      print_debug("[nonindep] We have NOT explored the right branch (isSecondBranch==false). Let's do it!"
          << " -- dec lev: " << decision_stack_.get_decision_level());
      Lit aLit = top_dec_lit();
      assert(decision_stack_.get_decision_level() > 0);
      decision_stack_.top().change_to_right_branch();
      reactivate_comps_and_backtrack_trail();
      print_debug("[nonindep] Flipping lit to: " << aLit.neg());
      setLiteralIfFree(aLit.neg());
      return RESOLVED;
    }
    isindependent = (decision_stack_.top().getbranchvar() < indep_support_end);
    comp_manager_->cacheModelCountOf(decision_stack_.top().super_comp(),
                                    decision_stack_.top().getTotalModelCount());
    //update cache scores
    stats.numcachedec_++;
    if (stats.numcachedec_ % 128 == 0) comp_manager_->rescale_cache_scores();
    comp_manager_->decreasecachescore(comp_manager_->getSuperComponentOf(decision_stack_.top()));

    if (decision_stack_.get_decision_level() <= 0) break;
    reactivate_comps_and_backtrack_trail();
    assert(decision_stack_.size() >= 2);
    if ((decision_stack_.top().getbranchvar() < indep_support_end))
      (decision_stack_.end() - 2)->includeSolution(decision_stack_.top().getTotalModelCount());
    else
      (decision_stack_.end() - 2)->includeSolution(decision_stack_.top().getTotalModelCount() > 0);
    print_debug("[nonindep] Backtracking from level " << decision_stack_.get_decision_level()
        << " count here is: " << decision_stack_.top().getTotalModelCount());
    decision_stack_.pop_back();
    isindependent = (decision_stack_.top().getbranchvar() < indep_support_end);
    // step to the next component not yet processed
    decision_stack_.top().nextUnprocessedComponent();

    assert( decision_stack_.top().remaining_comps_ofs() <
        comp_manager_->comp_stack_size() + 1);
  } while (true);
  return EXIT;
}

retStateT Counter::backtrack() {
  assert(decision_stack_.top().remaining_comps_ofs() <= comp_manager_->comp_stack_size());

  if (!isindependent && perform_projected_counting) return backtrack_nonindep();
  print_debug("[indep] Backtrack");
  do {
    print_debug("[indep] top count here: " << decision_stack_.top().getTotalModelCount());
    if (decision_stack_.top().branch_found_unsat()) {
      comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
    } else if (decision_stack_.top().anotherCompProcessible()) {
      print_debug("[indep] Processing another comp at dec lev "
          << decision_stack_.get_decision_level()
          << " instead of backtracking." << " Num unprocessed comps: "
          << decision_stack_.top().numUnprocessedComponents()
          << " so far the count: " << decision_stack_.top().getTotalModelCount());
      return PROCESS_COMPONENT;
    }

    // We have NOT explored the other side! Let's do it now!
    if (!decision_stack_.top().is_right_branch()) {
      print_debug("[indep] We have NOT explored the right branch (isSecondBranch==false). Let's do it!"
          << " -- dec lev: " << decision_stack_.get_decision_level());
      const Lit aLit = top_dec_lit();
      assert(decision_stack_.get_decision_level() > 0);
      decision_stack_.top().change_to_right_branch();
      reactivate_comps_and_backtrack_trail();
      print_debug("[indep] Flipping lit to: " << aLit.neg());
      setLiteralIfFree(aLit.neg());
      print_debug(COLORGBG "[indep] Backtrack finished -- we flipped the branch");
      return RESOLVED;
    }
    print_debug(COLORGBG "[indep] We have explored BOTH branches, actually BACKTRACKING."
        << " -- dec lev: " << decision_stack_.get_decision_level());
    comp_manager_->cacheModelCountOf(decision_stack_.top().super_comp(),
                                    decision_stack_.top().getTotalModelCount());


    //Cache score should be decreased since the component is getting added to cache
    if (config_.do_cache_score) {
      stats.numcachedec_++;
      if (stats.numcachedec_ % 128 == 0) comp_manager_->rescale_cache_scores();
      comp_manager_->decreasecachescore(
          comp_manager_->getSuperComponentOf(decision_stack_.top()));
    }

    // Backtrack from end, i.e. finished.
    if (decision_stack_.get_decision_level() == 0) {
      print_debug("[indep] Backtracking from lev 0, i.e. ending");
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
    print_debug("[indep] Backtracking from level " << decision_stack_.get_decision_level()
        << " count here is: " << decision_stack_.top().getTotalModelCount());
    decision_stack_.pop_back();
    isindependent = (decision_stack_.top().getbranchvar() < indep_support_end);
    auto& dst = decision_stack_.top();
    print_debug("[indep] -> Backtracked to level " << decision_stack_.get_decision_level()
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

void Counter::print_dec_info() const
{
  cout << "dec lits: " << endl;
  for(uint32_t i = 1; i < decision_stack_.size(); i ++) {
    Lit l = *(trail.begin()+ decision_stack_[i].trail_ofs());
    cout << "dec lev: " << std::setw(3) << i <<
      " lit: " << std::setw(6)
      << l
      << " is right: "
      << (int)decision_stack_[i].is_right_branch()
      << " ante: " << std::setw(10) << var(l).ante
      << " lev: " << var(l).decision_level
      << endl;
  }
}

void Counter::print_conflict_info() const
{
  print_dec_info();
  cout << "UIP cl lits: " << endl;
  for(uint32_t i = 0; i < uip_clause.size(); i ++) {
    const auto l = uip_clause[i];
    cout << "lit " << std::setw(6) << l
      << " lev: " << std::setw(4) << var(l).decision_level
      << " ante: " << std::setw(5) << std::left << var(l).ante
    << " val: " << lit_val_str(l) << endl;
  }
  cout << " ---- " << endl;
  cout << "top_dec_lit().neg(): " << top_dec_lit().neg() << endl;
  cout << "uip_clause[0]: " << uip_clause[0] << endl;
  cout << "uip_clause.front(): " << uip_clause.front() << endl;
}

void Counter::print_comp_stack_info() const {
    cout << "decision_stack_.top().remaining_comps_ofs(): "
      << decision_stack_.top().remaining_comps_ofs() << endl;
    cout << "comp_manager_->comp_stack_size(): " <<
      comp_manager_->comp_stack_size() << endl;
}

struct UIPFixer {
  UIPFixer(vector<Variable>& _vars) : vars(_vars){
  }
  bool operator()(const Lit& a, const Lit& b) const {
    auto a_dec = vars[a.var()].decision_level;
    auto b_dec = vars[b.var()].decision_level;
    if (a_dec != b_dec) return a_dec > b_dec;
    auto a_ante = vars[a.var()].ante;
    /* auto b_ante = vars[b.var()].ante; */
    if (!a_ante.isAnt()) return true;
    return false;
  }
  vector<Variable>& vars;
};

bool Counter::uip_clause_is_implied() {
    assert(sat_solver);
    vector<CMSat::Lit> lits;
    for(const auto& l: uip_clause) {
      lits.push_back(CMSat::Lit(l.var()-1, l.sign()));
    }
    cout << "to check lits: " << lits << endl;
    auto ret = sat_solver->solve(&lits);
    cout << "Ret: " << ret << endl;
    return ret == CMSat::l_False;
}

retStateT Counter::resolveConflict() {
  cout << "****** RECORD START" << endl;
  recordLastUIPCauses();
  cout << "*RECORD FINISHED*" << endl;
  act_inc *= 1.0/config_.act_exp;

  if (stats.conflicts > last_reduceDB_conflicts+10000) {
    reduceDB();
    if (stats.cls_deleted_since_compaction > 50000) compactConflictLiteralPool();
    last_reduceDB_conflicts = stats.conflicts;
  }
  print_conflict_info();
  VERBOSE_DEBUG_DO(cout << "NOW SORTING...." << endl);
  std::stable_sort(uip_clause.begin(), uip_clause.end(), UIPFixer(variables_));
  VERBOSE_DEBUG_DO(print_conflict_info());

  stats.conflicts++;
  assert(decision_stack_.top().remaining_comps_ofs() <= comp_manager_->comp_stack_size());
  // TODO deal with empty!
  if (uip_clause.empty()) { cout << "c EMPTY CLAUSE FOUND" << endl; }
  decision_stack_.top().mark_branch_unsat();
  assert(uip_clause.front() != NOT_A_LIT);

  cout << "backwards cleaning" << endl;
  print_comp_stack_info();
  uint32_t old_level = dec_level();
  int32_t backj = var(uip_clause.front()).decision_level;
  cout << "going back to lev: " << backj << " dec level now: " << dec_level()-1 << endl;
  while(dec_level()-1 > backj) {
    VERBOSE_DEBUG_DO(cout << "at dec lit: " << top_dec_lit() << endl);
    VERBOSE_DEBUG_DO(print_comp_stack_info());
    decision_stack_.top().mark_branch_unsat();
    reactivate_comps_and_backtrack_trail();
    decision_stack_.pop_back();
    decision_stack_.top().zero_out_branch_sol();
    comp_manager_->cleanRemainingComponentsOf(decision_stack_.top());
    comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
  }
  cout << "last dec lit: " << top_dec_lit() << endl;
  decision_stack_.top().mark_branch_unsat();
  decision_stack_.top().resetRemainingComps();
  print_comp_stack_info();
  VERBOSE_DEBUG_DO(cout << "DONE backw cleaning" << endl);
  VERBOSE_DEBUG_DO(print_conflict_info());

  Antecedent ant(NOT_A_CLAUSE);
  if (!uip_clause.empty()) {
    if (top_dec_lit().neg() == uip_clause[0]) {
      bool implied = uip_clause_is_implied();
      if (!implied) {
        cout << "old level: " << old_level << " went back to lev: " << backj << " dec level now: " << dec_level()-1 << endl;
        cout << "last dec lit: " << top_dec_lit() << endl;
        print_comp_stack_info();
        print_conflict_info();
        assert(false);
      }
      VERBOSE_DEBUG_DO(cout << "Setting reason the conflict cl" << endl);
      assert(var(uip_clause[0]).decision_level != -1);
      var(top_dec_lit().neg()).ante = addUIPConflictClause(uip_clause);
      ant = var(top_dec_lit()).ante;
      // oh wow, we set this decision variable to a propagated one
      // but we don't change its level! even though it's set due to previous level
    } else {
      addUIPConflictClause(uip_clause);
      stats.uip_not_added++;
    }
  }
  VERBOSE_DEBUG_DO(cout << "AFTER conflict, setup: ");
  VERBOSE_DEBUG_DO(print_conflict_info());
  VERBOSE_DEBUG_DO(cout << "is right here? " << decision_stack_.top().is_right_branch() << endl);

  if (decision_stack_.top().is_right_branch()) {
    // Backtracking since finished with this AND the other branch.
    return BACKTRACK;
  }

  assert(decision_stack_.get_decision_level() > 0);
  assert(decision_stack_.top().branch_found_unsat());

  // we do not have to remove pollutions here,
  // since conflicts only arise directly before
  // remaining comps are stored hence
  assert(decision_stack_.top().remaining_comps_ofs() == comp_manager_->comp_stack_size());

  decision_stack_.top().change_to_right_branch();
  const Lit lit = top_dec_lit();
  reactivate_comps_and_backtrack_trail();
  if (ant == Antecedent(NOT_A_CLAUSE)) {
    print_debug("Conflict pushes us to: " << lit<< " and due to failed literal probling, we can't guarantee it's due to the 1UIP, so setting it as a decision instead");
  } else {
    print_debug("Conflict pushes us to: " << lit);
  }
  setLiteralIfFree(lit.neg(), ant, false);
  if (ant.isAnt()) {
    cout << "lev: " << var(before_top_dec_lit()).decision_level << endl;
    cout << "other : " << var(lit).decision_level-1 << endl;
    int32_t dec_to_set = var(before_top_dec_lit()).decision_level;
    //Rewriting levels now.
    for(uint32_t i = decision_stack_[decision_stack_.size()-2].trail_ofs();
        i < trail.size(); i++) {
      cout << "dec level rewritten of lit: " << trail[i] << endl;
      var(trail[i]).decision_level = dec_to_set;
    }

    /* cout << "decision_stack_[decision_stack_.size()-2].trail_ofs(): " << decision_stack_[decision_stack_.size()-2].trail_ofs() << endl; */
    /* var(lit).decision_level = var(before_top_dec_lit()).decision_level; */
    assert(var(lit).decision_level == dec_to_set);
  }
  cout << "Returning from resolveConflict() with:";
  print_conflict_info();
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
  if (bSucceeded && config_.num_probe_multi > 0 && config_.failed_lit_probe_type != 0) {
    if (config_.failed_lit_probe_type == 2  &&
      (double)decision_stack_.size() >
        depth_q.getLongtTerm().avg()*config_.probe_only_after_ratio) {
      bSucceeded = failed_lit_probe();
    } else if (config_.failed_lit_probe_type == 1) {
      bSucceeded = failed_lit_probe();
    }
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
      if (isTrue(*p_otherLit)) {
        *it2++ = ClOffsBlckL(ofs, *p_otherLit);
        continue;
      }

      auto itL = beginOf(ofs) + 2;
      while (isFalse(*itL)) itL++;
      // either we found a free or satisfied lit
      if (*itL != SENTINEL_LIT) {
        litWatchList(*itL).addWatchLinkTo(ofs, *p_otherLit);
        std::swap(*itL, *p_watchLit);
      } else {
        // or unLit stays resolved
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

  // TODO get learnt clauses too
    /* for(auto cl_ofs: conflict_clauses_) { */
    /*     const ClHeader* ch = (const ClHeader *) */
    /*       ( &lit_pool_[cl_ofs - ClHeader::overheadInLits()]); */
    /*     auto sz = ch->length(); */
    /* } */
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

bool Counter::failed_lit_probe() {
  if (config_.bprop) return failed_lit_probe_with_bprop();
  else return failed_lit_probe_no_bprop();
}

bool Counter::failed_lit_probe_no_bprop()
{
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
            if (isUnknown(*lt) && !viewed_lits[lt->raw()]) {
              test_lits.push_back(*lt);
              print_debug("-> potential lit to test: " << lt->neg());
              viewed_lits[lt->raw()] = true;
            }
          }
        }
      }
    }
    num_curr_lits = trail.size() - trail_ofs;
    trail_ofs = trail.size();
    for (const auto& l: test_lits) viewed_lits[l.raw()] = false;

    // Figure out which literals to probe
    scores.clear();
    for (const auto& l: test_lits) scores.push_back(watches_[l].activity);
    sort(scores.begin(), scores.end());
    num_curr_lits = 10 + num_curr_lits / 20;
    num_curr_lits *= config_.num_probe_multi;
    double threshold = 0.0;
    if (scores.size() > num_curr_lits) {
      threshold = scores[scores.size() - num_curr_lits];
    }

    // Do the probing
    toSet.clear();
    for (auto& l : test_lits) if (isUnknown(l) && threshold <= watches_[l].activity) {
        if (!one_lit_probe(l, false)) return false;
        SLOW_DEBUG_DO(for(const auto& s: tmp_seen) assert(s == 0););
      }
  }
  print_debug(COLRED "Failed literal probing END -- no UNSAT, gotta check this branch");
  return true;
}

bool Counter::failed_lit_probe_with_bprop() {
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
            if (isUnknown(*lt) && !viewed_lits[lt->var()]) {
              test_lits.push_back(*lt);
              print_debug("-> potential lit to test: " << lt->neg());
              viewed_lits[lt->var()] = true;
            }
          }
        }
      }
    }
    num_curr_lits = trail.size() - trail_ofs;
    trail_ofs = trail.size();
    for (const auto& l: test_lits) viewed_lits[l.var()] = false;

    // Figure out which literals to probe
    scores.clear();
    for (const auto& l: test_lits) {
      scores.push_back(watches_[l].activity + watches_[l.neg()].activity);
    }
    sort(scores.begin(), scores.end());
    num_curr_lits = 5 + num_curr_lits / 40;
    num_curr_lits *= config_.num_probe_multi;
    double threshold = 0.0;
    if (scores.size() > num_curr_lits) {
      threshold = scores[scores.size() - num_curr_lits];
    }

    // Do the probing
    toSet.clear();
    for (auto& l : test_lits) {
      if (isUnknown(l) && threshold <=
          watches_[l].activity + watches_[l.neg()].activity) {
        if (watches_[l].activity < watches_[l.neg()].activity) l = l.neg();
        if (!one_lit_probe(l, true)) return false;
        if (isUnknown(l) && !one_lit_probe(l.neg(), false)) return false;
        SLOW_DEBUG_DO(for(const auto& s: tmp_seen) assert(s == 0););
      }
    }

    // Finally set what we came to set
    for(const auto& l: toSet) {
      if (isUnknown(l)) {
        auto sz = trail.size();
        setLiteralIfFree(l, Antecedent(Lit(), true));
        bool bSucceeded = propagate(sz);
        if (!bSucceeded) return false;
        stats.num_failed_bprop_literals_failed++;
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
    setLiteralIfFree(lit.neg(), Antecedent(Lit(), true));
    for(const auto& v: toClear) tmp_seen[v] = 0;
    toClear.clear();
    if (!propagate(sz)) {
      print_debug("Failed literal probing END -- this comp/branch is UNSAT");
      return false;
    }
  }

  return true;
}

void Counter::minimizeAndStoreUIPClause(Lit uipLit, vector<Lit> &cl) {
  tmp_clause_minim.clear();
  /*assertion_level_ = 0;
  for (const auto& lit : cl) {
    if (existsUnitClauseOf(lit.var())) continue;
    bool resolve_out = false;
    if (hasAntecedent(lit)) {
      resolve_out = true;
      if (getAntecedent(lit).isFake()) {
        // Probe/BProp is the reason
        for (int32_t i = 1; i < variables_[lit.var()].decision_level+1; i++) {
          const Lit l = trail[decision_stack_[i].trail_ofs()].neg();
          assert(decision_stack_[i].getbranchvar() == l.var());
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
  }*/

  /* if (uipLit.var()) { */
  /*   assert(var(uipLit).decision_level >= 0); */
  /*   /1* assert(var(uipLit).decision_level == decision_stack_.get_decision_level()); *1/ */
  /* } */


  tmp_clause_minim.clear();
  for(const auto& l:uip_clause) tmp_clause_minim.push_back(l);
  // Clearing
  for(const auto& v: toClear) tmp_seen[v] = false;
  toClear.clear();

  //assert(uipLit.var() != 0);
  stats.uip_lits_learned+=tmp_clause_minim.size();
  /* if (uipLit.var() != 0) { */
    stats.uip_lits_learned++;
    /* tmp_clause_minim.push_front(uipLit); */

    /* uint32_t lbd = calc_lbd(tmp_clause_minim); */
    /* if (lbd < 6) */
    if (stats.rem_lits_tried <= (200ULL*1000ULL) ||
        (stats.rem_lits_tried > (200ULL*1000ULL) &&
        ((double)stats.rem_lits_with_bins/(double)stats.rem_lits_tried > 3)))
      minimize_uip_cl_with_bins(tmp_clause_minim);
  /* } */
  stats.uip_cls++;
  stats.final_cl_sz+=tmp_clause_minim.size();
  uip_clause.clear();
  for(const auto& l: tmp_clause_minim) uip_clause.push_back(l);
}

void Counter::recordLastUIPCauses() {
  // note:
  // variables of lower dl: if seen we dont work with them anymore
  // variables of this dl: if seen we incorporate their
  // antecedent and set to unseen
  tmp_clause.clear();
  assert(toClear.empty());

  assertion_level_ = 0;
  uip_clause.clear();
  uip_clause.push_back(Lit(0, false));;
  Lit p = NOT_A_LIT;

  const int32_t DL = var(top_dec_lit()).decision_level;
  cout << "orig DL: " << decision_stack_.get_decision_level() << endl;
  cout << "new DL : " << DL << endl;
  print_dec_info();

  cout << "Doing loop:" << endl;
  uint32_t index = trail.size()-1;
  uint32_t pathC = 0;
  vector<Lit> c;
  do {
    if (confl.isAClause()) {
      assert(confl.asCl() != NOT_A_CLAUSE);
      c.clear();
      for(auto l = beginOf(confl.asCl()); *l != NOT_A_LIT; l++) {
        c.push_back(*l);
      }
    } else if (confl.isFake()) {
      assert(false);
    } else {
      assert(!confl.isAClause());
      c.clear();
      c.push_back(conflLit); //ONLY valid for 1st
      c.push_back(confl.asLit());
    }
    cout << "next cl: " << endl;
    for(const auto& l: c) {
        cout << std::setw(5) << l<< " lev: " << std::setw(3) << var(l).decision_level
          << " ante: " << std::setw(8) << var(l).ante
          << " val : " << std::setw(7) << lit_val_str(l)
          << endl;
    }

    cout << "For loop." << endl;
    for(uint32_t j = ((p == NOT_A_LIT) ? 0 : 1); j < c.size() ;j++) {
      Lit q = c[j];
      if (!tmp_seen[q.var()] && var(q).decision_level > 0){
        tmp_seen[q.var()] = 1;
        toClear.push_back(q.var());

        cout << std::setw(5) << q << " lev: " << std::setw(3) << var(q).decision_level
          << " ante: " << std::setw(8) << var(q).ante
          << " val : " << std::setw(7) << lit_val_str(q)
          << endl;
        if (var(q).decision_level >= DL) {
          pathC++;
          cout << "pathc inc." << endl;
        } else {
          uip_clause.push_back(q);
          cout << "added to cl." << endl;
        }
      }
    }
    cout << "PathC: " << pathC << endl;

    while (!tmp_seen[trail[index--].var()]) {}
    p     = trail[index+1];
    cout << "Next p: " << p << endl;
    confl = var(p).ante;
    tmp_seen[p.var()] = 0;
    pathC--;
  } while (pathC > 0);
  uip_clause[0] = p.neg();
  cout << "UIP cl: " << endl;
  for(const auto& l: uip_clause) {
      cout << std::setw(5) << l<< " lev: " << std::setw(3) << var(l).decision_level
        << " ante: " << std::setw(14) << var(l).ante
        << " val : " << std::setw(7) << lit_val_str(l)
        << endl;
  }

  //minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause);
  for(const auto& v: toClear) tmp_seen[v] = 0;
  toClear.clear();

  SLOW_DEBUG_DO(for(const auto& s: tmp_seen) assert(s == 0));
}

Counter::Counter(const CounterConfiguration& conf) : Instance(conf)
{
  mtrand.seed(conf.seed);
}

Counter::~Counter()
{
  delete comp_manager_;
}
