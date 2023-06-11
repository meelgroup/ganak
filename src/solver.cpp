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
#include <solver.h>
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
    setLiteral(lit, 0);
  }

  bool succeeded = propagate();
  release_assert(succeeded && "We ran CMS before, so it cannot be UNSAT");
  viewed_lits.resize(2*(nVars() + 1), 0);
  stats.num_unit_irred_clauses_ = unit_clauses_.size();
  init_decision_stack();
  qhead = 0;
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
  comp_manager_->initialize(watches_, alloc, longIrredCls, nVars());
}

void Counter::add_irred_cl(const vector<Lit>& lits) {
  if (lits.empty()) {
    cout << "ERROR: UNSAT should have been caught by external SAT solver" << endl;
    exit(-1);
  }
  for(const auto& l: lits) assert(l.var() <= nVars());
  stats.incorporateIrredClauseData(lits);
  Clause* cl = addClause(lits, false);
  auto off = alloc->get_offset(cl);
  if (cl) {
    if (off) longIrredCls.push_back(off);
    if (lits.size() >= 3) {
      for (const auto& l : lits)
        occ_lists_[l].push_back(off);
    }
  }


#ifdef SLOW_DEBUG
  debug_irred_cls.push_back(lits);
#endif
}

void Counter::add_red_cl(const vector<Lit>& lits, int lbd) {
  assert(ended_irred_cls);
  for(const auto& l: lits) release_assert(l.var() <= nVars());
  for(const auto& l: lits) release_assert(isUnknown(l));
  assert(lits.size() >= 2 && "No unit or empty clauses please");

  Clause* cl = addClause(lits, true);
  if (cl) {
    auto off = alloc->get_offset(cl);
    longRedCls.push_back(off);
    if (lbd == -1) lbd = lits.size();
    cl->lbd = lbd;
    cl->red = true;
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

  for(const auto& off: longIrredCls) {
    Clause& cl = *alloc->ptr(off);
    for(uint32_t i = 0; i < cl.sz; i++) {
      for(uint32_t i2 = i+1; i2 < cl.sz; i2++) {
        print_debug("v1: " << cl[i].var());
        print_debug("v2: " << cl[i2].var());
        primal.addEdge(cl[i].var(), cl[i2].var());
      }
    }
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
      << " branch var here: " << decision_stack_.at(dec_lev).var
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
    print_debug("var top of decision stack: " << decision_stack_.top().var);
    // NOTE: findNextRemainingComponentOf finds disjoint comps
    // we then solve them all with the decideLiteral & calling findNext.. again
    while (comp_manager_->findNextRemainingComponentOf(decision_stack_.top())) {
      // It's a component. It will ONLY fall into smaller pieces if we decide on a literal
      /* checkProbabilisticHashSanity(); -- no need, there is no way we get to 2**45 lookups*/
      if (restart_if_needed()) {return RESTART;}
      if (!decideLiteral()) {
        decision_stack_.top().nextUnprocessedComponent();
        continue;
      }
      VERBOSE_DEBUG_DO(print_all_levels());
      print_stat_line();

      while (!prop_and_probe()) {
        state = resolveConflict();
        while(state == GO_AGAIN) state = resolveConflict();
        if (state == BACKTRACK) break;
      }
      if (state == BACKTRACK) break;
    }
    // we are here because there is no next component, or we had to backtrack

    state = backtrack();
    if (state == EXIT) return SUCCESS;
    while (state != PROCESS_COMPONENT && !prop_and_probe()) {
      state = resolveConflict();
      while(state == GO_AGAIN) state = resolveConflict();
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
      } else polarity = standard_polarity(v);
    } else if (config_.polar_type == 1) polarity = standard_polarity(v);
    else if (config_.polar_type == 2) polarity = false;
    else if (config_.polar_type == 3) polarity = true;
    else assert(false);
  }
  return polarity;
}

bool Counter::decideLiteral() {
  print_debug("new decision level is about to be created, lev now: " << decision_stack_.get_decision_level() << " on path: " << decision_stack_.top().on_path_to_target_ << " branch: " << decision_stack_.top().is_right_branch());
  bool on_path = true;
  if (decision_stack_.size() != 1)
    on_path = decision_stack_.top().on_path_to_target_ && !decision_stack_.top().is_right_branch();

  decision_stack_.push_back(
    StackLevel(decision_stack_.top().currentRemainingComponent(),
               comp_manager_->comp_stack_size()));

  // The decision literal is now ready. Deal with it.
  uint32_t v = 0;
  isindependent = true;
  if (config_.branch_type == branch_t::gpmc) v = find_best_branch_gpmc(true);
  else v = find_best_branch(true);
  if (v == 0 && perform_projected_counting) {
    isindependent = false;
    if (config_.branch_type == branch_t::gpmc) v = find_best_branch_gpmc(false);
    else v = find_best_branch(false);
  }
  if (v == 0) {
    // we have set all remaining var(s) from a lower decision level.
    // so there is nothing to decide. Component has a single solution.
    VERBOSE_PRINT("We have set ALL REMAINING VARS FROM LOWER LEVELS!!");
    decision_stack_.pop_back();
    return false;
  }
  assert(val(v) == X_TRI);

  decision_stack_.top().var = v;
  decision_stack_.top().on_path_to_target_ = on_path;

  Lit lit = Lit(v, get_polarity(v));
  print_debug(COLYEL "decideLiteral() is deciding: " << lit << " dec level: "
      << decision_stack_.get_decision_level());
  setLiteral(lit, decision_stack_.get_decision_level());
  stats.decisions++;
  if (stats.decisions % 128 == 0) {
    comp_manager_->rescale_cache_scores();
  }
  assert( decision_stack_.top().remaining_comps_ofs() <= comp_manager_->comp_stack_size());
  return true;
}

uint32_t Counter::find_best_branch_gpmc(bool do_indep)
{
  uint32_t maxv = 0;
  double max_score_a = -1;
  double max_score_f = -1;
  double max_score_td = -1;

  for (auto it = comp_manager_->getSuperComponentOf(decision_stack_.top()).varsBegin();
      *it != varsSENTINEL; it++) if (!do_indep || *it < indep_support_end) {
    if (val(*it) != X_TRI) continue;
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
  vars_scores.clear();
  uint32_t best_var = 0;
  double best_var_score = -1;
  for (auto it = comp_manager_->getSuperComponentOf(decision_stack_.top()).varsBegin();
      *it != varsSENTINEL; it++) {
    if (val(*it) != X_TRI) continue;
    if (!do_indep || *it < indep_support_end) {
      const double score = scoreOf(*it) ;
      if (best_var_score == -1 || score > best_var_score) {
        best_var = *it;
        best_var_score = score;
      }
    }
  }

  if (config_.do_cache_score && best_var != 0) {
    double cachescore = comp_manager_->cacheScoreOf(best_var);
    for (auto it = comp_manager_->getSuperComponentOf(decision_stack_.top()).varsBegin();
         *it != varsSENTINEL; it++) {
      if (val(*it) != X_TRI) continue;
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
    const Lit dec_lit = Lit(dec.var, val(dec.var) == T_TRI);
    // Add decision
    if (i > 0) {
      const auto dec_lit2 = (target_polar[dec.var] ? 1 : -1)*(int)dec.var;
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
    const auto dec_lit = (target_polar[dst.var] ? 1 : -1)*(int)dst.var;
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
      reactivate_comps_and_backtrack_trail();
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
      while (decision_stack_.top().var >= indep_support_end) {
        print_debug("[nonindep] Going BACK because it's not independent and there is at least 1 solution");
        assert(!isindependent);
        if (decision_stack_.get_decision_level() <= 0) { break; }
        reactivate_comps_and_backtrack_trail();
        assert(decision_stack_.size() >= 2);
        assert(decision_stack_.top().getTotalModelCount() > 0);
        (decision_stack_.end() - 2)->includeSolution(
            decision_stack_.top().getTotalModelCount() > 0);
        decision_stack_.pop_back();
        isindependent = (decision_stack_.top().var < indep_support_end);
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
      setLiteral(aLit.neg(), decision_stack_.get_decision_level());
      return RESOLVED;
    }
    isindependent = (decision_stack_.top().var < indep_support_end);
    comp_manager_->cacheModelCountOf(decision_stack_.top().super_comp(),
                                    decision_stack_.top().getTotalModelCount());
    //update cache scores
    stats.numcachedec_++;
    if (stats.numcachedec_ % 128 == 0) comp_manager_->rescale_cache_scores();
    comp_manager_->decreasecachescore(comp_manager_->getSuperComponentOf(decision_stack_.top()));

    if (decision_stack_.get_decision_level() <= 0) break;
    reactivate_comps_and_backtrack_trail();
    assert(decision_stack_.size() >= 2);
    if ((decision_stack_.top().var < indep_support_end))
      (decision_stack_.end() - 2)->includeSolution(decision_stack_.top().getTotalModelCount());
    else
      (decision_stack_.end() - 2)->includeSolution(decision_stack_.top().getTotalModelCount() > 0);
    print_debug("[nonindep] Backtracking from level " << decision_stack_.get_decision_level()
        << " count here is: " << decision_stack_.top().getTotalModelCount());
    decision_stack_.pop_back();
    isindependent = (decision_stack_.top().var < indep_support_end);
    // step to the next component not yet processed
    decision_stack_.top().nextUnprocessedComponent();

    assert( decision_stack_.top().remaining_comps_ofs() <
        comp_manager_->comp_stack_size() + 1);
  } while (true);
  return EXIT;
}

uint64_t Counter::check_count(bool include_all_dec, int32_t single_var) {
    //let's get vars active
    set<uint32_t> active;

    const auto& s = decision_stack_.top();
    auto const& sup_at = s.super_comp();
    const auto& c = comp_manager_->at(sup_at);
#ifdef VERBOSE_DEBUG
    cout << "-> Checking count. Incl all dec: " << include_all_dec << " dec lev: " << decision_stack_.get_decision_level() << " var: " << single_var << endl;
    cout << "-> Variables in comp_manager_->at(" << sup_at << ")."
      << " num vars: " << c->nVars() << " vars: ";
    for(uint32_t i = 0; i < c->nVars(); i++) cout << c->varsBegin()[i] << " ";
    cout << endl;
#endif

    if (single_var == -1) {
      for(uint32_t i = 0; i < c->nVars(); i++) active.insert(c->varsBegin()[i]);
    } else {
      active.insert(single_var);
    }

#ifdef VERBOSE_DEBUG
    cout << "active: ";
    for(const auto&a: active) cout << a << " ";
    cout << endl;
#endif
    // Check dec level 0's
    vector<CMSat::Lit> cl;
    for(const auto& t: trail) {
      if (var(t).decision_level == 0) {
        cl.clear();
        cl.push_back(CMSat::Lit(t.var()-1, t.sign()));
        auto ret = sat_solver->solve(&cl);
        if (ret != CMSat::l_False) {
          cout << "ERROR: unit " << t << " is not correct!!" << endl;;
          assert(false);
        }
      }
    }

    // Checking
    VERBOSE_DEBUG_DO(print_trail());
    VERBOSE_DEBUG_DO(cout << "dec lev: " << decision_stack_.get_decision_level() << endl);
    VERBOSE_DEBUG_DO(cout << "top dec lit: " << top_dec_lit() << endl);
    CMSat::SATSolver s2;
    CMSat::copy_solver_to_solver(sat_solver, &s2);
    for(const auto& t: trail) {
      if (!include_all_dec) {
        if (var(t).decision_level >= decision_stack_.get_decision_level()) continue;
      }
      // don't include propagations or lev0 stuff
      if (!var(t).ante.isDecision() || var(t).decision_level == 0) continue;
      cl.clear();
      cl.push_back(CMSat::Lit(t.var()-1, !t.sign()));
      s2.add_clause(cl);
    }
    uint64_t num = 0;
    while(true) {
      auto ret = s2.solve();
      if (ret == CMSat::l_True) {
        num++;
        cl.clear();
        /* cout << "Blocking : "; */
        for(uint32_t i = 0; i < s2.nVars(); i++) {
          if (active.count(i+1)) {
            CMSat::Lit l = CMSat::Lit(i, s2.get_model()[i] == CMSat::l_True);
            cl.push_back(l);
            /* cout << l << " "; */
          }
        }
        /* cout << endl; */
        s2.add_clause(cl);
      } else if (ret == CMSat::l_False) break;
      else assert(false);
    }
    VERBOSE_DEBUG_DO(cout << "num                          : " << num << endl);
    if (single_var == -1) {
      VERBOSE_DEBUG_DO(cout << "ds.top().getTotalModelCount(): " << decision_stack_.top().getTotalModelCount() << endl);
      if (num != 0) assert(decision_stack_.top().getTotalModelCount() == num);
    }
    return num;
}

retStateT Counter::backtrack() {
  VERBOSE_DEBUG_DO(cout << "in " << __FUNCTION__ << " now " << endl);
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

    // We have NOT explored the other side and it hasn't been re-written to be
    // propagation.
    //    TODO: not sure we need the antecedent check here...
    //          probably doesn't hurt, but not sure.
    // Let's do it now!
    if (!decision_stack_.top().is_right_branch() &&
        var(top_dec_lit()).ante == Antecedent(NOT_A_CLAUSE)) {
      print_debug("[indep] We have NOT explored the right branch (isSecondBranch==false). Let's do it!"
          << " -- dec lev: " << decision_stack_.get_decision_level());
      const Lit aLit = top_dec_lit();
      assert(decision_stack_.get_decision_level() > 0);
      CHECK_COUNT_DO(check_count(true));
      SLOW_DEBUG_DO(assert(decision_stack_.top().get_rigth_model_count() == 0));
      decision_stack_.top().change_to_right_branch();
      // could be the flipped that's FALSEified so that would
      // mean the watchlist is not "sane". We need to propagate the flipped var and
      // then it'll be fine
      reactivate_comps_and_backtrack_trail(false);
      print_debug("[indep] Flipping lit to: " << aLit.neg());
      setLiteral(aLit.neg(), decision_stack_.get_decision_level());
      VERBOSE_DEBUG_DO(print_trail());
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

    CHECK_COUNT_DO(check_count());
    reactivate_comps_and_backtrack_trail(false);
    assert(decision_stack_.size() >= 2);
#ifdef VERBOSE_DEBUG
    const auto parent_count_before = (decision_stack_.end() - 2)->getTotalModelCount();
    const auto parent_count_before_left = (decision_stack_.end() - 2)->get_left_model_count();
    const auto parent_count_before_right = (decision_stack_.end() - 2)->get_rigth_model_count();
#endif
    (decision_stack_.end() - 2)->includeSolution(decision_stack_.top().getTotalModelCount());
    print_debug("[indep] Backtracking from level " << decision_stack_.get_decision_level()
        << " count here is: " << decision_stack_.top().getTotalModelCount());
    decision_stack_.pop_back();
    isindependent = (decision_stack_.top().var < indep_support_end);
    auto& dst = decision_stack_.top();
    print_debug("[indep] -> Backtracked to level " << decision_stack_.get_decision_level()
        // NOTE: -1 here because we have JUST processed the child
        //     ->> (see below nextUnprocessedComponent() call)
        << " num unprocessed comps here: " << dst.numUnprocessedComponents()-1
        << " current count here: " << dst.getTotalModelCount()
        << " branch: " << dst.is_right_branch()
        << " before including child it was: " <<  parent_count_before
        << " (left: " << parent_count_before_left
        << " right: " << parent_count_before_right << ")"
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
    uint32_t dvar = decision_stack_[i].var;
    Lit l = Lit(dvar, val(dvar) == T_TRI);
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

void Counter::print_cl(const vector<Lit>& cl) const {

  for(uint32_t i = 0; i < cl.size(); i ++) {
    const auto l = cl[i];
    cout << "lit " << std::setw(6) << l
      << " lev: " << std::setw(4) << var(l).decision_level
      << " ante: " << std::setw(5) << std::left << var(l).ante
    << " val: " << lit_val_str(l) << endl;
  }
}

void Counter::print_conflict_info() const
{
  print_dec_info();
  cout << "UIP cl lits: " << endl;
  print_cl(uip_clause);
  cout << " ---- " << endl;
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

size_t Counter::find_backtrack_level_of_learnt()
{
  if (uip_clause.size() == 0) return 0;
  else {
      uint32_t max_i = 0;
      for (uint32_t i = 0; i < uip_clause.size(); i++) {
          if (var(uip_clause[i]).decision_level > var(uip_clause[max_i]).decision_level)
              max_i = i;
      }
      std::swap(uip_clause[max_i], uip_clause[0]);
      return var(uip_clause[0]).decision_level;
  }
}

uint32_t Counter::find_lev_to_set(int32_t other_lev) {
  if (uip_clause.empty()) return 0;
  if (uip_clause.size() == 1) return 0;
  int32_t max_lev = 0;
  bool updated = false;
  uint32_t switch_to = 0;
  for (uint32_t i = 0; i < uip_clause.size(); i++) {
    int32_t lev = var(uip_clause[i]).decision_level;
      if (lev > max_lev && lev < other_lev) {
        max_lev = lev;
        updated = true;
        switch_to = i;
      }
  }
#ifdef VERBOSE_DEBUG
  cout << "max_lev: " << max_lev << " other_lev: " << other_lev << " updated: " << (int)updated << endl;
#endif
  assert(updated);
  std::swap(uip_clause[1], uip_clause[switch_to]);
  return max_lev;
}

void Counter::print_trail(bool check_entail, bool check_anything) const
{
  cout << "Current trail :" << endl;
  for(uint32_t i = 0; i < trail.size(); i++) {
    const auto l = trail[i];
    cout << "lit " << std::setw(6) << l
      << " lev: " << std::setw(4) << var(l).decision_level
      << " ante: " << std::setw(5) << std::left << var(l).ante
    << " val: " << std::setw(8) << lit_val_str(l)
    << " trail pos: " << std::setw(4) << i
    << " sublevel: "  << std::setw(3) << var(l).sublevel << endl;
  }
  if (check_anything) check_trail(check_entail);
}

void Counter::go_back_to(int32_t backj) {
  VERBOSE_DEBUG_DO(cout << "going back to lev: " << backj << " dec level now: " << decision_stack_.get_decision_level() << endl);
  while(decision_stack_.get_decision_level() > backj) {
    VERBOSE_DEBUG_DO(cout << "at dec lit: " << top_dec_lit() << endl);
    VERBOSE_DEBUG_DO(print_comp_stack_info());
    decision_stack_.top().mark_branch_unsat();
    decision_stack_.top().zero_out_all_sol(); //not sure it's needed
    comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
    reactivate_comps_and_backtrack_trail(false);
    decision_stack_.pop_back();
    // WOW, if this is ALL solutions, we get wrong count on NICE.cnf
    decision_stack_.top().zero_out_branch_sol();
    comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
    comp_manager_->cleanRemainingComponentsOf(decision_stack_.top());
  }
  VERBOSE_DEBUG_DO(print_comp_stack_info());
  VERBOSE_DEBUG_DO(cout << "DONE backw cleaning" << endl);
}

void Counter::check_trail([[maybe_unused]] bool check_entail) const {
  vector<uint32_t> num_decs_at_level(decision_stack_.get_decision_level()+1, 0);
  for(const auto& t: trail) {
    int32_t lev = var(t).decision_level;
    if (lev > decision_stack_.get_decision_level()) {
      cout << "Too high decision level enqueued." << endl;
      assert(false);
    }
    if (var(t).ante.isDecision() && lev > 0) {
      num_decs_at_level.at(lev)++;
      if (num_decs_at_level.at(lev) >= 2) {
        cout << "Two or more of decs at level: " << lev << endl;
        assert(false);
      }
    }
    if (val(t) != T_TRI) {
      assert(false && "Trail is wrong, trail[val] is not TRUE");
    }
    bool entailment_fail = false;
#ifdef CHECK_TRAIL_ENTAILMENT
    if (check_entail) {
      // Check entailment
      // No need to check if we are flipping and immediately backtracking
      if (!var(t).ante.isDecision()) {
        CMSat::SATSolver s2;
        CMSat::copy_solver_to_solver(sat_solver, &s2);
        vector<CMSat::Lit> cl;
        cl.push_back(CMSat::Lit(t.var()-1, t.sign())); //add opposite of implied
        s2.add_clause(cl);
        int32_t this_lev = var(t).decision_level;
        for(const auto& t2: trail) {
          if (var(t2).ante.isDecision() &&
              var(t2).decision_level <= this_lev &&
              var(t2).decision_level != 0) {
            cl.clear();
            cl.push_back(CMSat::Lit(t2.var()-1, !t2.sign())); // add all decisions (non-negated)
            s2.add_clause(cl);
          }
        }
        auto ret = s2.solve();
        if (ret != CMSat::l_False) {
          cout << "Not implied by decisions: " << t << endl;
          entailment_fail = true;
        }
      }
    }
#endif
    assert(!entailment_fail);
  }
}

bool Counter::is_implied(const vector<Lit>& cl) {
    assert(sat_solver);
    vector<CMSat::Lit> lits;
    for(const auto& l: cl) {
      lits.push_back(CMSat::Lit(l.var()-1, l.sign()));
    }
    VERBOSE_DEBUG_DO(cout << "to check lits: " << lits << endl);
    auto ret = sat_solver->solve(&lits);
    VERBOSE_DEBUG_DO(cout << "Ret: " << ret << endl);
    return ret == CMSat::l_False;
}


void Counter::check_implied(const vector<Lit>& cl) {
  bool implied = is_implied(cl);
  if (!implied) {
    cout << "ERROR, not implied" << endl;
    cout << "last dec lit: " << top_dec_lit() << endl;
    print_comp_stack_info();
    print_conflict_info();
    assert(false);
  }
}

retStateT Counter::resolveConflict() {
  VERBOSE_DEBUG_DO(cout << "****** RECORD START" << endl);
  VERBOSE_DEBUG_DO(print_trail());
  recordLastUIPCauses();
  assert(uip_clause.front() != NOT_A_LIT);
  VERBOSE_DEBUG_DO(cout << "*RECORD FINISHED*" << endl);
  act_inc *= 1.0/config_.act_exp;

  if (stats.conflicts-stats.uip_not_added+stats.saved_uip_used > last_reduceDB_conflicts+10000) {
    reduceDB();
    if (stats.cls_deleted_since_compaction > 50000) alloc->consolidate(this);
    last_reduceDB_conflicts = stats.conflicts-stats.uip_not_added+stats.saved_uip_used;
  }
  VERBOSE_DEBUG_DO(print_conflict_info());

  stats.conflicts++;
  assert(decision_stack_.top().remaining_comps_ofs() <= comp_manager_->comp_stack_size());
  decision_stack_.top().zero_out_branch_sol();
  decision_stack_.top().mark_branch_unsat();

  VERBOSE_DEBUG_DO(cout << "backwards cleaning" << endl);
  VERBOSE_DEBUG_DO(print_comp_stack_info());
  int32_t backj = find_backtrack_level_of_learnt();
  int32_t lev_to_set = find_lev_to_set(backj);

  // This is DEFINITELY a decision, right?
  VERBOSE_PRINT("backj initially: " << backj);
  /* assert(variables_[decision_stack_.at(backj).var].ante == Antecedent(NOT_A_CLAUSE)); */
  lev_to_set = std::min(lev_to_set, backj);
  bool flipped = (
      uip_clause[0].neg().var() == decision_stack_.at(backj).var
       && lev_to_set+1 == backj);
  if (!flipped) {
    /* cout << "Not adding. backj: " << backj << " lev_to_set: " << lev_to_set << endl; */
    /* print_trail(false); */
    /* print_conflict_info(); */
    /* cout << "--------------" << endl << endl; */
    if (config_.do_save_uip && uip_clause.size() > 1 &&
        uip_clause[0].neg().var() == decision_stack_.at(backj).var) {
      if (saved_uip_cls.size() <= (uint32_t)lev_to_set)
        saved_uip_cls.resize(lev_to_set+1);
      saved_uip_cls[lev_to_set] = uip_clause;
      stats.saved_uip++;
    }

    stats.uip_not_added++;
    decision_stack_.top().resetRemainingComps();
    if (decision_stack_.top().is_right_branch()) {
      return BACKTRACK;
    } else {
      decision_stack_.top().change_to_right_branch();
      auto lit = top_dec_lit();
      reactivate_comps_and_backtrack_trail(false);
      setLiteral(lit.neg(), decision_stack_.get_decision_level(), Antecedent(NOT_A_CLAUSE));
      return RESOLVED;
    }
  }
  assert(flipped);
  VERBOSE_DEBUG_DO(cout << "after finding backj lev: " << backj << " lev_to_set: " << lev_to_set <<  endl);
  VERBOSE_DEBUG_DO(print_conflict_info());

  go_back_to(backj);
  VERBOSE_DEBUG_DO(print_conflict_info());
  VERBOSE_PRINT("decision_stack_.get_decision_level(): " << decision_stack_.get_decision_level());

  Antecedent ant(NOT_A_CLAUSE);
  assert(!uip_clause.empty());
  SLOW_DEBUG_DO(check_implied(uip_clause));
  if (decision_stack_.get_decision_level() > 0 &&
      top_dec_lit().neg() == uip_clause[0]) {
    VERBOSE_PRINT("FLIPPING. Setting reason the conflict cl");
    assert(var(uip_clause[0]).decision_level != -1);
    ant = addUIPConflictClause(uip_clause);
    var(top_dec_lit()).ante = ant;
  }
  VERBOSE_PRINT("Ant is :" << ant);
  VERBOSE_PRINT("AFTER conflict, setup: ");
  VERBOSE_DEBUG_DO(print_conflict_info());
  VERBOSE_PRINT("is right here? " << decision_stack_.top().is_right_branch());

  comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
  decision_stack_.top().zero_out_branch_sol();
  decision_stack_.top().mark_branch_unsat();
  decision_stack_.top().resetRemainingComps();

  if (decision_stack_.top().is_right_branch()) {
    var(uip_clause[0]).decision_level = lev_to_set; //TODO what to do with sublevel?
    var(uip_clause[0]).ante = ant;
    lit_values_[uip_clause[0]] = T_TRI;
    lit_values_[uip_clause[0].neg()] = F_TRI;
    trail[var(uip_clause[0]).sublevel] = uip_clause[0];
    qhead = std::min(qhead, var(uip_clause[0]).sublevel);

    reactivate_comps_and_backtrack_trail(false);
    bool ret = propagate();
    if (!ret) return GO_AGAIN;

#ifdef VERBOSE_DEBUG
    cout << "FLIPPED Returning from resolveConflict() with:";
    print_conflict_info();
    print_trail(false); // we re-written the level above, so entailment
                        // may fail. when backtracking it'll be fine, though
    cout << "We have already counted this LEFT branch, so we backtrack now." << endl;
#endif
    return BACKTRACK;
  }

  if (decision_stack_.get_decision_level() > 0) {
    assert(decision_stack_.top().remaining_comps_ofs() == comp_manager_->comp_stack_size());
  }

  decision_stack_.top().change_to_right_branch();
  reactivate_comps_and_backtrack_trail(false);
  setLiteral(uip_clause[0], lev_to_set, ant);

#ifdef VERBOSE_DEBUG
  cout << "Returning from resolveConflict() with:";
  print_conflict_info();
  print_trail();
#endif

  return RESOLVED; // will ALWAYS propagate afterwards.
}

bool Counter::clause_falsified(const vector<Lit>& cl) const {
  for(const auto&l: cl) {
    if (val(l) != F_TRI) return false;
  }
  return true;
}


bool Counter::clause_asserting(const vector<Lit>& cl) const {
  uint32_t num_unkn = 0;
  for(const auto&l: cl) {
    if (val(l) == T_TRI) return false;
    if (val(l) == X_TRI) {num_unkn++; if (num_unkn>=2) return false;}
  }
  return true;
}


bool Counter::clause_satisfied(const vector<Lit>& cl) const {
  for(const auto&l: cl) {
    if (val(l) == T_TRI) return true;
  }
  return false;
}

/* #define VERB_DEBUG_SAVED */

bool Counter::prop_and_probe() {
  VERBOSE_DEBUG_DO(cout << "in " << __FUNCTION__ << " now. " << endl);
  // the asserted literal has been set, so we start
  // bcp on that literal
  assert(trail.size() > 0 && "Mate added this, but it seems OK");

  bool bSucceeded;
prop_again:;
  bSucceeded = propagate();
  if (config_.do_save_uip && bSucceeded &&
      (int)saved_uip_cls.size() > decision_stack_.get_decision_level() &&
      !saved_uip_cls[decision_stack_.get_decision_level()].empty()) {
    auto& cl = saved_uip_cls[decision_stack_.get_decision_level()];
    std::sort(cl.begin(), cl.end(),
        [=](const Lit& a, const Lit& b) {
          if (val(a) == X_TRI || val(b) == X_TRI) {
            if (val(a) == val(b)) return false;
            if (val(a) == X_TRI) return true;
            if (val(b) == X_TRI) return false;
          }
          if (var(a).decision_level == var(b).decision_level) {
            return val(a) == T_TRI;
          }
          return var(a).decision_level > var(b).decision_level;
        });
#ifdef VERB_DEBUG_SAVED
    cout << "-----" << endl;
    cout << "here now. dec lev: " << decision_stack_.get_decision_level() << endl;
    print_dec_info();
    print_trail();
    cout << "CL: " << endl;
    print_cl(cl);
#endif
    // The below often MUST have two from the same decision level at the top. Otherwise,
    // we can have something like this:
    /* lit -15    lev: 28   ante: DEC             val: TRUE */
    /* lit -6     lev: 27   ante: CL:        4921 val: FALSE */
    /* lit -54    lev: 27   ante: CL:        2338 val: FALSE */
    /* lit -13    lev: 24   ante: DEC             val: FALSE */
    /* lit 118    lev: 15   ante: Lit:        -85 val: FALSE */
    // which would have propagated at 27. Not OK.

    if (clause_falsified(cl)) {
      auto ante = addUIPConflictClause(cl);
      if (ante.isAClause()) setConflictState(alloc->ptr(ante.asCl()));
      else setConflictState(cl[0], cl[1]);
#ifdef VERB_DEBUG_SAVED
      cout << "Falsified. Ante: " << ante << endl;
#endif
      stats.saved_uip_used++;
      saved_uip_cls[decision_stack_.get_decision_level()].clear();
      return false;
    } else if (var(cl[0]).decision_level == var(cl[1]).decision_level && clause_asserting(cl)) {
      auto ante = addUIPConflictClause(cl);
      setLiteral(cl[0], decision_stack_.get_decision_level(), ante);
#ifdef VERB_DEBUG_SAVED
      cout << "Asserting. Now cl: " << endl;
      print_cl(cl);
#endif
      stats.saved_uip_used++;
      saved_uip_cls[decision_stack_.get_decision_level()].clear();
      goto prop_again;
    } else if (var(cl[0]).decision_level == var(cl[1]).decision_level && clause_satisfied(cl)) {
      addUIPConflictClause(cl);
#ifdef VERB_DEBUG_SAVED
      cout << "Satisfied." << endl;
#endif
      stats.saved_uip_used++;
      saved_uip_cls[decision_stack_.get_decision_level()].clear();
    } else {
#ifdef VERB_DEBUG_SAVED
      cout << "Throwing this saved away." << endl;
#endif
      stats.saved_uip_thrown++;
    }
  }

  if (bSucceeded && config_.num_probe_multi > 0 && config_.failed_lit_probe_type != 0) {
    if (config_.failed_lit_probe_type == 2  &&
      (double)decision_stack_.size() >
        depth_q.getLongtTerm().avg()*config_.probe_only_after_ratio) {
      bSucceeded = failed_lit_probe();
    } else if (config_.failed_lit_probe_type == 1) {
      bSucceeded = failed_lit_probe();
    }
  }
  VERBOSE_DEBUG_DO(cout << "Returning " << bSucceeded << " from " << __FUNCTION__ << endl);
  return bSucceeded;
}

template<uint32_t start>
inline void Counter::get_maxlev_maxind(ClauseOfs ofs, int32_t& maxlev, uint32_t& maxind)
{
  Clause& cl = *alloc->ptr(ofs);
  for(auto i3 = start; i3 < cl.sz; i3++) {
    Lit l = cl[i3];
    int32_t nlev = var(l).decision_level;
    VERBOSE_DEBUG_DO(cout << "i3: " << i3 << " l : " << l << " var(l).decision_level: "
        << var(l).decision_level << " maxlev: " << maxlev << endl);
    if (nlev > maxlev) {maxlev = nlev; maxind = i3;}
  }
}

bool Counter::propagate() {
  confl = Antecedent(NOT_A_CLAUSE);
  VERBOSE_PRINT("qhead in propagate(): " << qhead << " trail sz: " << trail.size());
  for (; qhead < trail.size(); qhead++) {
    const Lit unLit = trail[qhead].neg();
    const int32_t lev = var(unLit).decision_level;
    VERBOSE_PRINT("&&Propagating: " << unLit.neg() << " qhead: " << qhead << " lev: " << lev);

    //Propagate bin clauses
    for (const auto& l : litWatchList(unLit).binary_links_) {
      if (val(l) == F_TRI) {
        setConflictState(unLit, l);
        VERBOSE_DEBUG_DO(cout << "Bin confl. otherlit: " << l << endl);
      } else if (val(l) == X_TRI) {
        setLiteral(l, lev, Antecedent(unLit));
        VERBOSE_DEBUG_DO(cout << "Bin prop: " << l << " lev: " << lev << endl);
      /* } else if (val(l) == T_TRI && var(l).decision_level > lev) { */
      /*   var(l).ante = Antecedent(unLit); */
      /*   VERBOSE_PRINT("Updated ante of " << l << " to: " << unLit); */
      /*   /1* var(l).decision_level = lev; *1/ */
      }
    }

    //Propagate long clauses
    auto& ws = litWatchList(unLit).watch_list_;

#if 0
    cout << "prop-> will go through norm cl:" << endl;
    for(const auto& w: ws) {
      cout << "norm cl offs: " << w.ofs << " cl: ";
      const auto ofs = w.ofs;
      for(Lit* c = beginOf(ofs); *c != NOT_A_LIT; c++) { cout << *c << " "; }
      cout << endl;
    }
    cout << " will do it now... " << endl;
#endif

    auto it2 = ws.begin();
    auto it = ws.begin();
    for (; it != ws.end(); it++) {
      if (isTrue(it->blckLit)) { *it2++ = *it; continue; }

      const auto ofs = it->ofs;
      Clause& c = *alloc->ptr(ofs);
      if (c[0] == unLit) { std::swap(c[0], c[1]); }

#ifdef VERBOSE_DEBUG
      cout << "Prop Norm cl: " << ofs << endl;
      for(const auto&l: c) {
        cout << "lit " << std::setw(6) << l
          << " lev: " << std::setw(4) << var(l).decision_level
          << " ante: " << std::setw(5) << std::left << var(l).ante
          << " val: " << lit_val_str(l) << endl;
      }
#endif

      assert(c[1] == unLit);
      if (isTrue(c[0])) {
        *it2++ = ClOffsBlckL(ofs, c[0]);
        continue;
      }

      uint32_t i = 2;
      for(; i < c.sz; i++) if (!isFalse(c[i])) break;
      // either we found a free or satisfied lit
      if (i != c.sz) {
        c[1] = c[i];
        c[i] = unLit;
        VERBOSE_PRINT("New watch for cl: " << c[1]);
        litWatchList(c[1]).addWatchLinkTo(ofs, c[0]);
      } else {
        *it2++ = *it;
        if (val(c[0]) == F_TRI) {
          VERBOSE_PRINT("Conflicting state from norm cl offs: " << ofs);
          if (lev != decision_stack_.get_decision_level()) {
            int32_t maxlev = lev;
            uint32_t maxind = 1;
            get_maxlev_maxind(ofs, maxlev, maxind);
            if (maxind != 1) {
              VERBOSE_PRINT("swapping. maxlev: " << maxlev << " maxind: " << maxind << " c[1]: " << c[1] << " c[maxind]: " << c[maxind]);
              std::swap(c[1], c[maxind]);
              it2--; // undo last watch
              litWatchList(c[1]).addWatchLinkTo(ofs, it->blckLit);
            }
          }
          setConflictState(&c);
          it++;
          break;
        } else {
          assert(val(c[0]) == X_TRI);
          VERBOSE_PRINT("prop long lev: " << lev << " dec_stack.get_lev : " << decision_stack_.get_decision_level());
          if (lev == decision_stack_.get_decision_level()) {
            setLiteral(c[0], lev, Antecedent(ofs));
            VERBOSE_PRINT("Norm long prop: " << c[0] << " lev: " << lev);
          } else {
            int32_t maxlev = lev;
            uint32_t maxind = 1;
            get_maxlev_maxind(ofs, maxlev, maxind);
            if (maxind != 1) {
                std::swap(c[1], c[maxind]);
                it2--; // undo last watch
                litWatchList(c[1]).addWatchLinkTo(ofs, it->blckLit);
            }
            setLiteral(c[0], maxlev, Antecedent(ofs));
            VERBOSE_DEBUG_DO(cout << "Weird long prop: " << c[0] << " lev: " << maxlev << endl);
          }
        }
      }
    }
    while(it != ws.end()) *it2++ = *it++;
    ws.resize(it2-ws.begin());
    if (confl != Antecedent(NOT_A_CLAUSE)) break;
  }
  SLOW_DEBUG_DO(if (confl == Antecedent(NOT_A_CLAUSE) && !check_watchlists()) {
      print_trail(false, false);assert(false);});
  SLOW_DEBUG_DO(if (confl == Antecedent(NOT_A_CLAUSE)) check_all_propagated());
  VERBOSE_PRINT("After propagate, qhead is: " << qhead);
  return confl == Antecedent(NOT_A_CLAUSE);
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

  uint32_t trail_ofs = trail_at_top();
  uint32_t num_curr_lits = 0;
  while (trail_ofs < trail.size()) {
    test_lits.clear();
    for (auto it = trail.begin() + trail_ofs; it != trail.end(); it++) {
      // Only going through the long, the binary clauses have set the variables already
      for (auto cl_ofs : occ_lists_[it->neg()]) {
        if (!isSatisfied(cl_ofs)) {
          for(const auto& l: *alloc->ptr(cl_ofs)) {
            if (isUnknown(l) && !viewed_lits[l.raw()]) {
              test_lits.push_back(l);
              print_debug("-> potential lit to test: " << l.neg());
              viewed_lits[l.raw()] = true;
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

  uint32_t trail_ofs = trail_at_top();
  uint32_t num_curr_lits = 0;
  while (trail_ofs < trail.size()) {
    test_lits.clear();
    for (auto it = trail.begin() + trail_ofs; it != trail.end(); it++) {
      // Only going through the long, the binary clauses have set the variables already
      for (const auto& off : occ_lists_[it->neg()]) {
        if (!isSatisfied(off)) {
          for (auto& l: *alloc->ptr(off)) {
            if (isUnknown(l) && !viewed_lits[l.var()]) {
              test_lits.push_back(l);
              print_debug("-> potential lit to test: " << l.neg());
              viewed_lits[l.var()] = true;
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
        setLiteral(l, decision_stack_.get_decision_level(), Antecedent(Lit(), true));
        bool bSucceeded = propagate();
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
  setLiteral(lit, decision_stack_.get_decision_level());
  assert(!hasAntecedent(lit));
  if (set == true) assert(toClear.empty());

  bool bSucceeded = propagate();
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
  assert(false && "qhead is not working for probing");
  if (!set) {
    for(const auto& v: toClear) tmp_seen[v] = 0;
    toClear.clear();
  }

  if (!bSucceeded) {
    stats.num_failed_literals_detected_++;
    print_debug("-> failed literal detected");
    sz = trail.size();
    setLiteral(lit.neg(), decision_stack_.get_decision_level(), Antecedent(Lit(), true));
    for(const auto& v: toClear) tmp_seen[v] = 0;
    toClear.clear();
    if (!propagate()) {
      print_debug("Failed literal probing END -- this comp/branch is UNSAT");
      return false;
    }
  }

  return true;
}

bool Counter::litRedundant(const Lit p, uint32_t abstract_levels) {
    VERBOSE_PRINT(__func__ << " called");

    analyze_stack.clear();
    analyze_stack.push_back(p);

    Lit tmpLit[2];
    uint32_t size;
    size_t top = toClear.size();
    while (!analyze_stack.empty()) {
      VERBOSE_PRINT("At point in litRedundant: " << analyze_stack.back());
      const auto reason = var(analyze_stack.back()).ante;
      assert(reason.isAnt());  //Must have a reason
      analyze_stack.pop_back();

      Lit* lits = NULL;
      size = 0;
      if (reason.isAClause()) {
        Clause* cl = alloc->ptr(reason.asCl());
        lits = cl->getData();
        size = cl->sz;
#ifdef VERBOSE_DEBUG
        cout << "CL offs: " << reason.asCl() << " in analyze_stack:" << endl;
        for(Lit* l = lits; *l != SENTINEL_LIT; l++) {
          cout << std::setw(5) << *l<< " lev: " << std::setw(3) << var(*l).decision_level
            << " ante: " << std::setw(8) << var(*l).ante
            << " val : " << std::setw(7) << lit_val_str(*l)
            << endl;
        }
#endif
      } else if (reason.isFake()) {
        assert(false && "not handled right now");
      } else if (!reason.isAClause()) {
        tmpLit[0] = NOT_A_LIT;
        tmpLit[1] = reason.asLit();
        size = 2;
#ifdef VERBOSE_DEBUG
        cout << "Bin cl in analyze_stack: " << endl;
        Lit* l = tmpLit;
        cout << std::setw(5) << *l<< " lev: " << std::setw(3) << var(*l).decision_level
          << " ante: " << std::setw(8) << var(*l).ante
          << " val : " << std::setw(7) << lit_val_str(*l)
          << endl;
#endif
        lits = tmpLit;
      } else {assert(false && "no such reason");}

      for (uint32_t i = 1; i < size; i++) {
        Lit p2 = lits[i];
        VERBOSE_PRINT("Examining lit " << p2 << " tmp_seen: " << tmp_seen[p2.var()]);
        if (!tmp_seen[p2.var()] && var(p2).decision_level > 0) {
          if (var(p2).ante.isAnt()
              && (abstractLevel(p2.var()) & abstract_levels) != 0
          ) {
              VERBOSE_PRINT("lit " << p2 << " OK");
              tmp_seen[p2.var()] = 1;
              analyze_stack.push_back(p2);
              toClear.push_back(p2.var());
          } else {
              VERBOSE_PRINT("lit " << p2 << " NOT OK");
              //Return to where we started before function executed
              for (size_t j = top; j < toClear.size(); j++) tmp_seen[toClear[j]] = 0;
              toClear.resize(top);
              return false;
          }
        }
      }
    }
    VERBOSE_PRINT("Returning OK from " << __func__);
    return true;
}

uint32_t Counter::abstractLevel(const uint32_t x) const
{
    return ((uint32_t)1) << (variables_[x].decision_level & 31);
}

void Counter::recursiveConfClauseMin()
{
    VERBOSE_DEBUG_DO(print_conflict_info());
  VERBOSE_PRINT("recursive ccmin now.");
  uint32_t abstract_level = 0;
  for (size_t i = 1; i < uip_clause.size(); i++) {
    //(maintain an abstraction of levels involved in conflict)
    abstract_level |= abstractLevel(uip_clause[i].var());
  }

  size_t i, j;
  for (i = j = 1; i < uip_clause.size(); i++) {
    if (var(uip_clause[i]).ante == Antecedent(NOT_A_CLAUSE)
      || !litRedundant(uip_clause[i], abstract_level)
    ) {
      VERBOSE_PRINT("ccmin -- keeping lit: " << uip_clause[i]);
      uip_clause[j++] = uip_clause[i];
    } else {
      VERBOSE_PRINT("ccmin -- NOT keeping lit: " << uip_clause[i]);
    }
  }
  uip_clause.resize(j);
}

void Counter::minimizeUIPClause() {
  stats.uip_cls++;
  stats.orig_uip_lits += uip_clause.size();
  recursiveConfClauseMin();
  for(const auto& c: toClear) tmp_seen[c] = 0;
  toClear.clear();

  SLOW_DEBUG_DO(check_implied(uip_clause));
  tmp_clause_minim.clear();
  for(const auto& l:uip_clause) tmp_clause_minim.push_back(l);

  stats.uip_lits_ccmin+=tmp_clause_minim.size();
  if (stats.rem_lits_tried <= (200ULL*1000ULL) ||
      (stats.rem_lits_tried > (200ULL*1000ULL) &&
      ((double)stats.rem_lits_with_bins/(double)stats.rem_lits_tried > 3)))
    minimize_uip_cl_with_bins(tmp_clause_minim);
  stats.final_cl_sz+=tmp_clause_minim.size();
  uip_clause.clear();
  for(const auto& l: tmp_clause_minim) uip_clause.push_back(l);
  SLOW_DEBUG_DO(check_implied(uip_clause));
}

int32_t Counter::get_confl_maxlev(const Lit p) const {
  Lit tmpLit[3];
  tmpLit[2] = SENTINEL_LIT;
  Lit* c;
  uint32_t size = 0;

  if (confl.isAClause()) {
    assert(confl.asCl() != NOT_A_CLAUSE);
    VERBOSE_PRINT("Conflicting CL offset: " << confl.asCl());
    Clause* cl = alloc->ptr(confl.asCl());
    c = cl->getData();
    size = cl->sz;
  } else if (confl.isFake()) {
    assert(false);
  } else {
    //Binary
    c = tmpLit;
    if (p == NOT_A_LIT) c[0] = conflLit;
    else c[0] = p;
    c[1] = confl.asLit();
    size = 2;
  }
  int32_t maxlev = -1;
  for(uint32_t i = 0; i < size; i ++) {
#ifdef VERBOSE_DEBUG
    cout << "confl cl[" << std::setw(5) << i << "]"
        << " lit: " << c[i]
        << " lev: " << std::setw(3) << var(c[i]).decision_level
        << " ante: " << std::setw(8) << var(c[i]).ante
        << " val : " << std::setw(7) << lit_val_str(c[i])
        << endl;
#endif
    if (var(c[i]).decision_level > maxlev) maxlev = var(c[i]).decision_level;
  }
  VERBOSE_DEBUG_DO(cout << "maxlev: " << maxlev << endl);
  VERBOSE_DEBUG_DO(print_dec_info());
  assert(var(c[0]).decision_level == var(c[1]).decision_level);
  return maxlev;
}

// Returns TRUE if it can generate a UIP. Otherwise, false
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

  SLOW_DEBUG_DO(for(const auto& t:tmp_seen) assert(t == 0););
  int32_t DL = var(top_dec_lit()).decision_level;
  VERBOSE_DEBUG_DO(cout << "orig DL: " << decision_stack_.get_decision_level() << endl);
  VERBOSE_DEBUG_DO(cout << "new DL : " << DL << endl);
  VERBOSE_DEBUG_DO(print_dec_info());
  int32_t maxlev = get_confl_maxlev(p);
  go_back_to(maxlev);
  DL = var(top_dec_lit()).decision_level;

  Lit tmpLit[3];
  tmpLit[2] = SENTINEL_LIT;
  Lit* c;
  uint32_t size;
  VERBOSE_DEBUG_DO(cout << "Doing loop:" << endl);
  int32_t index = trail.size()-1;
  uint32_t pathC = 0;
  do {
    if (confl.isAClause()) {
      // Long clause
      assert(confl.asCl() != NOT_A_CLAUSE);
      Clause& cl = *alloc->ptr(confl.asCl());
      if (cl.red && cl.lbd > lbd_cutoff) {
        cl.increaseScore();
        cl.update_lbd(calc_lbd(&cl));
      }
      if (p == NOT_A_LIT) std::swap(c[0], c[1]);
      c = cl.getData();
      size = cl.sz;
    } else if (confl.isFake()) {
      assert(false);
    } else {
      // Binary
      assert(!confl.isAClause());
      c = tmpLit;
      if (p == NOT_A_LIT) c[0] = conflLit;
      else c[0] = p;
      c[1] = confl.asLit();
      if (p == NOT_A_LIT && var(c[0]).decision_level < var(c[1]).decision_level)
        std::swap(c[0], c[1]);
      size = 2;
    }
    VERBOSE_DEBUG_DO(cout << "next cl: " << endl);
#ifdef VERBOSE_DEBUG
    for(uint32_t i = 0; i < size; i++) {
      Lit l = c[i];
      cout << std::setw(5) << l<< " lev: " << std::setw(3) << var(l).decision_level
        << " ante: " << std::setw(8) << var(l).ante
        << " val : " << std::setw(7) << lit_val_str(l)
        << endl;
    }
#endif
    int32_t nDecisionLevel = var(c[0]).decision_level;
    VERBOSE_DEBUG_DO(cout << "nDecisionLevel: " <<  nDecisionLevel << endl);
    if (p == NOT_A_LIT) assert(nDecisionLevel == DL);

    VERBOSE_DEBUG_DO(cout << "For loop." << endl);
    for(uint32_t j = ((p == NOT_A_LIT) ? 0 : 1); j < size ;j++) {
      Lit q = c[j];
      if (!tmp_seen[q.var()] && var(q).decision_level > 0){
        increaseActivity(q);
        tmp_seen[q.var()] = 1;
        toClear.push_back(q.var());
#ifdef VERBOSE_DEBUG
        cout << std::setw(5) << q << " lev: " << std::setw(3) << var(q).decision_level
          << " ante: " << std::setw(8) << var(q).ante
          << " val : " << std::setw(7) << lit_val_str(q)
          << endl;
#endif
        if (var(q).decision_level >= nDecisionLevel) {
          pathC++;
          VERBOSE_DEBUG_DO(cout << "pathc inc." << endl);
        } else {
          uip_clause.push_back(q);
          VERBOSE_DEBUG_DO(cout << "added to cl." << endl);
        }
      }
    }
    VERBOSE_DEBUG_DO(cout << "PathC: " << pathC << endl);

    do {
      while (!tmp_seen[trail[index--].var()]) { SLOW_DEBUG_DO(assert(index >= 0));};
      p = trail[index+1];
      assert(p != NOT_A_LIT);
#ifdef VERBOSE_DEBUG
      cout << "going back on trail: " << std::setw(5) << p<< " lev: " << std::setw(3) << var(p).decision_level
        << " ante: " << std::setw(8) << var(p).ante
        << " val : " << std::setw(7) << lit_val_str(p)
        << endl;
#endif
    } while(var(trail[index+1]).decision_level < nDecisionLevel);
    VERBOSE_DEBUG_DO(cout << "Next p: " << p << endl);
    confl = var(p).ante;
    tmp_seen[p.var()] = 0;
    pathC--;
  } while (pathC > 0);
  assert(pathC == 0);
  uip_clause[0] = p.neg();
#ifdef VERBOSE_DEBUG
  cout << "UIP cl: " << endl;
  for(const auto& l: uip_clause) {
      cout << std::setw(5) << l<< " lev: " << std::setw(3) << var(l).decision_level
        << " ante: " << std::setw(14) << var(l).ante
        << " val : " << std::setw(7) << lit_val_str(l)
        << endl;
  }
#endif
  SLOW_DEBUG_DO(check_implied(uip_clause));
  minimizeUIPClause();
  SLOW_DEBUG_DO(for(const auto& s: tmp_seen) assert(s == 0));
}

Counter::Counter(const CounterConfiguration& conf) : Instance(conf) {
  mtrand.seed(conf.seed);
}

Counter::~Counter() {
  delete comp_manager_;
}

bool Counter::check_watchlists() const {
  bool ret = true;
#if 0
  // All watchlists
  cout << "All watchlists: " << endl;
  for(uint32_t i = 2; i < (nVars()+1)*2; i++) {
    Lit lit = Lit(i/2, i%2);
    cout << "->Watchlist for lit " << lit << " (val: " << lit_val_str(lit) << ") " << endl;
    auto& ws = watches_[lit].watch_list_;
    for(const auto& w: ws) {
      const auto ofs = w.ofs;
      cout << "--> Cl ofs " << ofs << " lits: ";
      for(Lit const* c = beginOf(ofs); *c != NOT_A_LIT; c++) {
        cout << *c << " (val: " << lit_val_str(*c) << ") ";
      }
      cout << endl;
    }
  }
#endif

  // Also check that after propagation, if the clause is not satisfied,
  // it's NOT the case that prop queue contains
  // FALSE & UNK. Must be UNK & UNK
  for(uint32_t i = 2; i < (nVars()+1)*2; i++) {
    Lit lit = Lit(i/2, i%2);
    auto& ws = watches_[lit].watch_list_;
    for(const auto& w: ws) {
      const auto ofs = w.ofs;
      uint32_t num_unk = 0;
      bool sat = false;
      for(const auto& l: *alloc->ptr(ofs)) {
        if (isUnknown(l)) num_unk++;
        if (isTrue(l)) sat = true;
      }
      if (!sat && num_unk >=2 && !isUnknown(lit)) {
        cout << "ERROR, we are watching a FALSE: " << lit << ", but there are at least 2 UNK in cl offs: " << ofs << " clause: " << endl;
      for(const auto& l: *alloc->ptr(ofs)) {
          cout << l << " (val: " << lit_val_str(l)
            << " lev: " << var(l).decision_level << ") " << endl;
        }
        ret = false;
      }
    }
  }
  return ret;
}

#ifdef SLOW_DEBUG
void Instance::check_all_propagated() const {
  // Everything that should have propagated, propagated
  for(const auto& cl: debug_irred_cls) {
    Lit unk = NOT_A_LIT;
    uint32_t num_unknown = 0;
    bool satisfied = false;
    for(const auto& l: cl) {
      if (isTrue(l)) {satisfied = true; break;}
      if (isUnknown(l)) {num_unknown++; unk = l;}
      if (num_unknown > 1) break;
    }

    if (!satisfied && num_unknown == 1) {
      cout << "ERROR! Clause: " << cl << " should have propagated: " << unk << endl;
      assert(false);
    }
  }
}
#endif
