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

#include "counter.h"

#include <algorithm>
#include <complex>
#include <ios>
#include <iomanip>
#include <numeric>
#include <utility>
#include "common.h"
#include "comp_types/comp.h"
#include "cryptominisat5/cryptominisat.h"
#include "cryptominisat5/solvertypesmini.h"
#include "primitive_types.h"
#include "counter_config.h"
#include "stack.h"
#include "structures.h"
#include "time_mem.h"
#include "IFlowCutter.h"
#include "graph.hpp"

void Counter::simplePreProcess()
{
  for (auto lit : unit_clauses_) {
    assert(!isUnitClause(lit.neg()) && "Formula is not UNSAT, we ran CMS before");
    if (val(lit) == X_TRI) setLiteral(lit, 0);
    assert(val(lit) == T_TRI);
  }

  bool succeeded = propagate();
  release_assert(succeeded && "We ran CMS before, so it cannot be UNSAT");
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

void Counter::init_activity_scores() {
  for (auto l = Lit(1, false); l != watches_.end_lit(); l.inc()) {
    watches_[l].activity = watches_[l].binary_links_.size();
  }
  for(const auto& off: longIrredCls) {
    const auto& cl = *alloc->ptr(off);
    for(const auto& l: cl) watches_[l].activity++;
  }
}

void Counter::end_irred_cls()
{
  seen.clear();
  seen.resize(2*(nVars()+2), 0);
  delete comp_manager_;
  comp_manager_ = new ComponentManager(conf,stats, lit_values_, indep_support_end, this);
  comp_manager_->getrandomseedforclhash();

  // reset stats
  depth_q.clearAndResize(conf.first_restart);
  cache_miss_rate_q.clearAndResize(conf.first_restart);
  comp_size_q.clearAndResize(conf.first_restart);
  next_print_stat_cache = 2ULL*1000LL*1000LL;

  stats.maximum_cache_size_bytes_ = conf.maximum_cache_size_MB*1024*1024;
  init_decision_stack();
  simplePreProcess();
  ended_irred_cls = true;

  if (conf.verb) stats.printShortFormulaInfo();
  // This below will initialize the disjoint component analyzer (ana_)
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
  if (cl) longIrredCls.push_back(off);
  SLOW_DEBUG_DO(debug_irred_cls.push_back(lits));
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
    assert(cl->red);
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
  bool conditionOnCNF = indep_support_end > 3 && nVars() > 20 && nVars() <= conf.td_varlim;
  if (!conditionOnCNF) {
    verb_print(1, "skipping TD, too many/few vars. Setting branch to fallback");
    conf.branch_type = conf.branch_fallback_type;
    return;
  }

  Graph primal(nVars()+1);
  for(uint32_t i = 2; i < (nVars()+1)*2; i++) {
    Lit l(i/2, i%2);
    for(const auto& l2: watches_[l].binary_links_) {
      if ((!l2.red() || (l2.red() && conf.td_with_red_bins))
          && l < l2.lit()) {
        print_debug("v1: " << l.var());
        print_debug("v2: " << l2.var());
        primal.addEdge(l.var(), l2.lit().var());
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
  /*     density <= conf.td_denselim && */
  /*     edge_var_ratio <= conf.td_ratiolim; */

  /* if (!conditionOnPrimalGraph) { */
  /*   verb_print(1, "skipping td, primal graph is too large or dense." */
  /*       " Setting branch to fallback"); */
  /*   conf.branch_type = conf.branch_fallback_type; */
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

vector<CMSat::Lit> ganak_to_cms_cl(const vector<Lit>& cl) {
  vector<CMSat::Lit> cms_cl;
  for(const auto& l: cl) cms_cl.push_back(CMSat::Lit(l.var()-1, !l.sign()));
  return cms_cl;
}

vector<CMSat::Lit> ganak_to_cms_cl(const Lit& l) {
  vector<CMSat::Lit> cms_cl;
  cms_cl.push_back(CMSat::Lit(l.var()-1, !l.sign()));
  return cms_cl;
}

mpz_class Counter::check_norestart(const vector<Lit>& cube) {
  // Test
  CounterConfiguration conf2 = conf;
  conf2.do_restart = 0;
  conf2.verb = 0;
  vector<Lit> tmp;
  Counter* test_cnt = new Counter(conf2);
  test_cnt->new_vars(nVars());
  // Long cls
  for(const auto& off: longIrredCls) {
    const Clause& cl = *alloc->ptr(off);
    tmp.clear();
    for(const auto& l: cl) tmp.push_back(l);
    test_cnt->add_irred_cl(tmp);
  }
  // Bin cls
  for(uint32_t i = 2; i < (nVars()+1)*2; i++) {
    Lit l(i/2, i%2);
    for(const auto& l2: watches_[l].binary_links_) {
      if (l2.irred() && l < l2.lit()) {
        tmp.clear();
        tmp.push_back(l);
        tmp.push_back(l2.lit());
        test_cnt->add_irred_cl(tmp);
      }
    }
  }
  // Unit cls
  for(const auto& l: unit_clauses_) {
    tmp.clear();
    tmp.push_back(l);
    test_cnt->add_irred_cl(tmp);
  }
  // The cube
  for(const auto&l: cube) {
    tmp.clear();
    tmp.push_back(l.neg());
    test_cnt->add_irred_cl(tmp);
  }
  test_cnt->end_irred_cls();
  vector<Cube> ret;
  test_cnt->count(ret);
  assert(ret.size() == 1);
  delete test_cnt;
  return ret[0].val;
}

mpz_class Counter::outer_count(CMSat::SATSolver* _sat_solver) {
  mpz_class val = 0;
  sat_solver = _sat_solver;

  auto ret = sat_solver->solve();
  int32_t num_runs = 0;
  start_time = cpuTime();
  vector<Cube> cubes;
  if (conf.do_restart) {
    while(ret == CMSat::l_True) {
      vector<Cube> tmp_cubes;
      count(tmp_cubes);
      for(const auto&c: tmp_cubes) {val+=c.val; cubes.push_back(c);}
      num_runs++;
      cout << "Num runs: " << num_runs << endl;
      cout << "Cubes     : " << endl;
      for(const auto&c: tmp_cubes) cout << "-> c: " << c.cube << " val: " << c.val << endl;
#ifdef SLOW_DEBUG
      for(const auto& c: tmp_cubes) {
        auto check_cnt = check_norestart(c.cube);
        cout << "check cube: " << c.cube << " cnt here: " << check_cnt << " c.val: " << c.val << endl;
        assert(check_cnt == c.val);
      }
#endif

      // Add cubes to CMS
      for(const auto& c: tmp_cubes) sat_solver->add_clause(ganak_to_cms_cl(c.cube));
      ret = sat_solver->solve();
      if (ret == CMSat::l_False) break;

      // Add cubes to counter
      for(const auto& c: tmp_cubes) add_irred_cl(c.cube);
      end_irred_cls();
    }
  } else if (ret == CMSat::l_True) {
      count(cubes);
      assert(cubes.size() == 1);
      assert(cubes[0].cube.empty());
      val = cubes[0].val;
  }
  return val;
}

void Counter::count(vector<Cube>& ret_cubes) {
  release_assert(ret_cubes.empty());
  release_assert(ended_irred_cls && "ERROR *must* call end_irred_cls() before solve()");
  if (indep_support_end == std::numeric_limits<uint32_t>::max()) indep_support_end = nVars()+2;
  mini_cubes.clear();
  verb_print(1, "Sampling set size: " << indep_support_end-1);

  // Only compute TD decomposition once
  if (tdscore.empty() && conf.td_with_red_bins) {
    if (conf.branch_type == branch_t::sharptd ||
        conf.branch_type == branch_t::gpmc) td_decompose();
    verb_print(1, "branch type: " << conf.get_branch_type_str());
  }

  const auto exit_state = countSAT();
  if (conf.verb) stats.printShort(this, &comp_manager_->get_cache());
  if (exit_state == RESTART) {
    ret_cubes = mini_cubes;
  } else {
    assert(exit_state == SUCCESS);
    Cube c(vector<Lit>(), decision_stack_.top().getTotalModelCount());
    ret_cubes.push_back(c);
  }
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
  if (conf.verb) {
    verb_print(1, "GANAK time so far: " << (cpuTime() - start_time));
    stats.printShort(this, &comp_manager_->get_cache());
  }

  next_print_stat_cache = stats.num_cache_look_ups_ + (4ULL*1000LL*1000LL);
  next_print_stat_confl = stats.conflicts + (30LL*1000LL);
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

    // Here we can vivify I think
    if (conf.do_vivify) {
      vivify_clauses();
      bool ret = propagate();
      assert(ret);
    }
  }
  return SUCCESS;
}

bool Counter::standard_polarity(const uint32_t v) const {
    return watches_[Lit(v, true)].activity > watches_[Lit(v, false)].activity;
}

bool Counter::get_polarity(const uint32_t v) const
{
  bool polarity;
  if (conf.polar_type == 0) {
    if (var(Lit(v, false)).set_once) {
      polarity = var(Lit(v, false)).last_polarity;
    } else polarity = standard_polarity(v);
  } else if (conf.polar_type == 1) polarity = standard_polarity(v);
  else if (conf.polar_type == 4) polarity = !standard_polarity(v);
  else if (conf.polar_type == 2) polarity = false;
  else if (conf.polar_type == 3) polarity = true;
  else assert(false);
  return polarity;
}

bool Counter::decideLiteral() {
  print_debug("new decision level is about to be created, lev now: " << decision_stack_.get_decision_level() << " branch: " << decision_stack_.top().is_right_branch());
  decision_stack_.push_back(
    StackLevel(decision_stack_.top().currentRemainingComponent(),
               comp_manager_->comp_stack_size()));

  // The decision literal is now ready. Deal with it.
  uint32_t v = 0;
  isindependent = true;
  if (conf.branch_type == branch_t::gpmc) v = find_best_branch_gpmc(true);
  else v = find_best_branch(true);
  if (v == 0 && perform_projected_counting) {
    isindependent = false;
    if (conf.branch_type == branch_t::gpmc) v = find_best_branch_gpmc(false);
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

  if (conf.do_cache_score && best_var != 0) {
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

bool Counter::compute_cube(vector<Lit>& cube, mpz_class& cube_val) {
  assert(conf.do_restart);
  cube.clear();
  print_debug(COLWHT "-- " __func__ " BEGIN");
  print_debug_noendl(COLWHT "Decisions in the cube: ");

  // add decisions, components, and counts
  cube_val = decision_stack_.top().getTotalModelCount();
  bool error = false;
  for(uint32_t i = 1; i < decision_stack_.size()-1; i++) {
    const StackLevel& dec = decision_stack_[i];
    const Lit dec_lit = Lit(dec.var, val(dec.var) == T_TRI);
    // Add decisions
    cube.push_back(dec_lit.neg());
    print_debug_noendl(dec_lit.neg() << " ");
    if (dec.getTotalModelCount() > 0) cube_val *= dec.getTotalModelCount();
  }
  // Get a solution
  for(uint32_t i = 1; i < decision_stack_.size()-1; i++) {
    const StackLevel& dec = decision_stack_[i];
    const auto off_start = dec.remaining_comps_ofs();
    const auto off_end = dec.getUnprocessedComponentsEnd();
    // add all but the last component (it's the one we just counted)
    for(uint32_t i2 = off_start; i2 < off_end-1; i2++) {
      const auto& c = comp_manager_->at(i2);
      for(auto v = c->varsBegin(); *v != varsSENTINEL; v++) {
        cube.push_back(Lit(*v, sat_solver->get_model()[*v-1] == CMSat::l_True));
      }
    }
  }
  assert(!error);
  print_debug_noendl(endl);

#ifdef VERBOSE_DEBUG
  // Show decision stack's comps
  for(size_t i = 0; i < decision_stack_.size(); i++) {
    const auto& dst = decision_stack_.at(i);
    /* const auto dec_lit2 = trail[ds.trail_ofs()]; */
    /* cout << "dec_lit2: " << dec_lit2 << endl; */
    print_debug(COLWHT "decision_stack.at(" << i << "):"
      << " decision var: " << dst.var
      << " num unproc comps: " << dst.numUnprocessedComponents()
      << " unproc comps end: " << dst.getUnprocessedComponentsEnd()
      << " remain comps offs: " << dst.remaining_comps_ofs()
      << " count here: " << dst.getTotalModelCount()
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

  cout << COLWHT "cube so far. Size: " << cube.size() << " cube: ";
  for(const auto& l: cube) cout << l << " ";
  cout << endl;
  print_debug(COLWHT "cube's SOLE count: " << decision_stack_.top().getTotalModelCount());
  print_debug(COLWHT "cube's RECORDED count: " << cube_val);
  assert(!error);
#endif

  return true;
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

  if (!conf.do_restart) return false;
  bool restart = false;
  if (conf.restart_type == 0
      && comp_size_q.isvalid() &&
      comp_size_q.avg() < comp_size_q.getLongtTerm().avg()*conf.restart_cutoff_mult)
    restart = true;
  if (conf.restart_type == 1
      && cache_miss_rate_q.isvalid() &&
      cache_miss_rate_q.avg() > cache_miss_rate_q.getLongtTerm().avg()*0.95)
    restart = true;

  if (conf.restart_type == 2
      && depth_q.isvalid() &&
      depth_q.avg() > depth_q.getLongtTerm().avg()*(1.0/conf.restart_cutoff_mult))
    restart = true;

  if (conf.restart_type == 3 &&
      (stats.decisions-stats.last_restart_num_decisions) > conf.next_restart)
    restart = true;

  if (conf.restart_type == 4 && stats.cache_hits_misses_q.isvalid()
      && stats.cache_hits_misses_q.avg() <
      stats.cache_hits_misses_q.getLongtTerm().avg()*conf.restart_cutoff_mult)
      restart = true;

  if (conf.restart_type == 5 && stats.comp_size_times_depth_q.isvalid() &&
        stats.comp_size_times_depth_q.avg() >
          stats.comp_size_times_depth_q.getLongtTerm().avg()*(1.0/conf.restart_cutoff_mult))
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

    assert(mini_cubes.empty());
    while (decision_stack_.size() > 1) {
      cout
        << endl << "------" << endl
        << " lev: " << decision_stack_.get_decision_level()
        << " cnt: " << decision_stack_.top().getTotalModelCount()
        << endl;

      vector<Lit> cube;
      mpz_class cube_val;
      bool ret = compute_cube(cube, cube_val);
      cout << "->> cube here: " << cube << " val: " << cube_val << " ret: " << (int)ret << endl;
      if (ret && cube_val != 0) mini_cubes.push_back(Cube(cube, cube_val));
      comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
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
          cout << "ERROR: unit " << t << " is not correct!!" << endl;
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
      if (!var(t).ante.isNull() || var(t).decision_level == 0) continue;
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
        for(uint32_t i = 0; i < s2.nVars(); i++) {
          if (active.count(i+1)) {
            CMSat::Lit l = CMSat::Lit(i, s2.get_model()[i] == CMSat::l_True);
            cl.push_back(l);
          }
        }
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
    if (!decision_stack_.top().is_right_branch() && var(top_dec_lit()).ante.isNull()) {
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
    if (conf.do_cache_score) {
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
        << " right: " << parent_count_before_right << ")");

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

void Counter::print_conflict_info() const
{
  print_dec_info();
  cout << "UIP cl lits: " << endl;
  print_cl(uip_clause);
  cout << "uip_clause[0]: " << uip_clause[0] << endl;
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
    if (var(t).ante.isNull() && lev > 0) {
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
      if (!var(t).ante.isNull()) {
        CMSat::SATSolver s2;
        CMSat::copy_solver_to_solver(sat_solver, &s2);
        vector<CMSat::Lit> cl;
        cl.push_back(CMSat::Lit(t.var()-1, t.sign())); //add opposite of implied
        s2.add_clause(cl);
        int32_t this_lev = var(t).decision_level;
        for(const auto& t2: trail) {
          if (var(t2).ante.isNull() &&
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

Counter::ConflictData Counter::find_conflict_level() {
	ConflictData data;
  Lit* c;
  uint32_t size;
  fill_cl(confl, c, size, NOT_A_LIT);
	data.nHighestLevel = var(c[0]).decision_level;
  int32_t curr_dl = var(top_dec_lit()).decision_level;
	if (data.nHighestLevel == curr_dl && var(c[1]).decision_level == curr_dl) return data;

	int highestId = 0;
    data.bOnlyOneLitFromHighest = true;
	// find the largest decision level in the clause
	for (uint32_t nLitId = 1; nLitId < size; ++nLitId) {
		int32_t nLevel = var(c[nLitId]).decision_level;
		if (nLevel > data.nHighestLevel) {
			highestId = nLitId;
			data.nHighestLevel = nLevel;
			data.bOnlyOneLitFromHighest = true;
		} else if (nLevel == data.nHighestLevel && data.bOnlyOneLitFromHighest == true) {
			data.bOnlyOneLitFromHighest = false;
		}
	}

  // fixing clause & watchlist
  assert(false && "TODO below");
	/* if (highestId != 0) { */
	/* 	std::swap(c[0], c[highestId]); */
	/* 	if (highestId > 1 && size > 2) { */
      /* assert(confl.isAClause()); */
      /* ClauseOfs off = confl.asCl(); */
	/* 		remove(watches_[c[highestId].neg()].watch_list_, ClOffsBlckL(off, c[1])); */
	/* 		ws[~conflCls[0]].push(Watcher(cind, conflCls[1])); */
	/* 	} */
	/* } */
	return data;
}

retStateT Counter::resolveConflict() {
  VERBOSE_DEBUG_DO(cout << "****** RECORD START" << endl);
  VERBOSE_DEBUG_DO(print_trail());


  recordLastUIPCauses();
  assert(uip_clause.front() != NOT_A_LIT);
  VERBOSE_DEBUG_DO(cout << "*RECORD FINISHED*" << endl);
  act_inc *= 1.0/conf.act_exp;

  if (stats.conflicts-stats.uip_not_added+stats.saved_uip_used > last_reduceDB_conflicts+10000) {
    reduceDB();
    if (stats.cls_deleted_since_compaction > 30000 && alloc->consolidate(this)) {
        stats.cls_deleted_since_compaction = 0;
    }
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
    /* cout << "--------------" << endl << endl; */
    /* print_trail(false); */
    /* print_conflict_info(); */
    /* cout << "Not adding. backj: " << backj << " lev_to_set: " << lev_to_set */
    /*   << " current lev: " << decision_stack_.get_decision_level() << endl; */
    if (uip_clause.size() == 1) {
      if (!existsUnitClauseOf(uip_clause[0])) unit_clauses_.push_back(uip_clause[0]);
    }
    if (conf.do_save_uip && uip_clause.size() > 1) {
      if (saved_uip_cls.size() <= (uint32_t)lev_to_set) saved_uip_cls.resize(lev_to_set+1);
        // Latest seems better than smallest, so just upgrade to latest
        if (!saved_uip_cls[lev_to_set].empty()) stats.saved_uip_thrown++;
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
      setLiteral(lit.neg(), decision_stack_.get_decision_level(), Antecedent::fakeAnte());
      return RESOLVED;
    }
  }
  assert(flipped);
  VERBOSE_DEBUG_DO(cout << "after finding backj lev: " << backj << " lev_to_set: " << lev_to_set <<  endl);
  VERBOSE_DEBUG_DO(print_conflict_info());

  go_back_to(backj);
  VERBOSE_DEBUG_DO(print_conflict_info());
  VERBOSE_PRINT("decision_stack_.get_decision_level(): " << decision_stack_.get_decision_level());

  Antecedent ant;
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


template<class T> bool Counter::clause_satisfied(const T& cl) const {
  for(const auto&l: cl) if (val(l) == T_TRI) return true;
  return false;
}

/* #define VERB_DEBUG_SAVED */

Counter::SavedUIPRet Counter::deal_with_saved_uips() {
  auto DL = decision_stack_.get_decision_level();
  if ( (int)saved_uip_cls.size() <= DL || saved_uip_cls[DL].empty())
    return Counter::SavedUIPRet::cont;

  int32_t i = DL;
  /* for(; i < (int32_t)saved_uip_cls.size(); i++) { */
  /*   if (!saved_uip_cls[i].empty()) break; */
  /* } */
  /* if (i == (int32_t)saved_uip_cls.size()) return Counter::SavedUIPRet::cont; */
  auto& cl = saved_uip_cls[i];
  std::sort(cl.begin(), cl.end(),
    [=](const Lit& a, const Lit& b) {
      if (val(a) == X_TRI || val(b) == X_TRI) {
        if (val(a) == val(b)) return false;
        if (val(a) == X_TRI) return true;
        if (val(b) == X_TRI) return false;
      }
      if (var(a).decision_level == var(b).decision_level) {
        if(val(a) != val(b)) return val(a) == T_TRI;
        return false;
      }
      return var(a).decision_level > var(b).decision_level;
    });

  uint32_t num_unk = 0;
  uint32_t num_true = 0;
  uint32_t num_false = 0;
  for(const auto& l: cl) {
    switch(val(l)) {
      case X_TRI: num_unk++; break;
      case T_TRI: num_true++; break;
      case F_TRI: num_false++; break;
    }
  }

#ifdef VERB_DEBUG_SAVED
  cout << "-----" << endl; cout << "here now. dec lev: " << DL << endl;
  print_dec_info();
  print_trail();
#endif
  // The below often MUST have two from the same decision level at the top. Otherwise,
  // we can have something like this:
  /* lit -15    lev: 28   ante: DEC             val: TRUE */
  /* lit -6     lev: 27   ante: CL:        4921 val: FALSE */
  /* lit -54    lev: 27   ante: CL:        2338 val: FALSE */
  /* lit -13    lev: 24   ante: DEC             val: FALSE */
  /* lit 118    lev: 15   ante: Lit:        -85 val: FALSE */
  // which would have propagated at 27. Not OK.

  if (num_false == cl.size() &&
      var(cl[0]).decision_level == DL &&
      var(cl[1]).decision_level == DL) {
    auto ante = addUIPConflictClause(cl);
    if (ante.isAClause()) setConflictState(alloc->ptr(ante.asCl()));
    else setConflictState(cl[0], cl[1]);
#ifdef VERB_DEBUG_SAVED
    cout << "Falsified. Ante: " << ante << endl; print_cl(cl);
#endif
    stats.saved_uip_used++;
    stats.saved_uip_used_falsified++;
    saved_uip_cls[DL].clear();
    return Counter::SavedUIPRet::ret_false;
  } else if (var(cl[1]).decision_level == DL
      && num_unk == 1 && num_false == cl.size()-1) {
    auto ante = addUIPConflictClause(cl);
    setLiteral(cl[0], DL, ante);
#ifdef VERB_DEBUG_SAVED
    cout << "Asserting. Now cl: " << endl; print_cl(cl);
#endif
    stats.saved_uip_used++;
    stats.saved_uip_used_asserting++;
    saved_uip_cls[DL].clear();
    return Counter::SavedUIPRet::prop_again;
  } else if ((num_true >= 1 || num_unk >= 2) &&
      (var(cl[0]).decision_level == var(cl[1]).decision_level || num_unk > 0)) {
    addUIPConflictClause(cl);
#ifdef VERB_DEBUG_SAVED
    cout << "At least 2 unknown or Satisfied." << endl; print_cl(cl);
#endif
    stats.saved_uip_used++;
    stats.saved_uip_used_sat_or_unk++;
    saved_uip_cls[DL].clear();
  } else {
#ifdef VERB_DEBUG_SAVED
    cout << "Throwing this saved away." << endl; print_cl(cl);
#endif
    saved_uip_cls[DL].clear();
    stats.saved_uip_thrown++;
  }

  return Counter::SavedUIPRet::cont;
}

bool Counter::prop_and_probe() {
  VERBOSE_DEBUG_DO(cout << "in " << __FUNCTION__ << " now. " << endl);
  // the asserted literal has been set, so we start
  // bcp on that literal
  assert(trail.size() > 0 && "Mate added this, but it seems OK");

  bool bSucceeded;
prop_again:;
  bSucceeded = propagate();
  if (conf.do_save_uip && bSucceeded) {
    switch(deal_with_saved_uips()) {
      case Counter::SavedUIPRet::prop_again: goto prop_again;
      case Counter::SavedUIPRet::ret_false: return false;
      case Counter::SavedUIPRet::cont: ;
    }
  }
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
  confl = Antecedent();
  VERBOSE_PRINT("qhead in propagate(): " << qhead << " trail sz: " << trail.size());
  for (; qhead < trail.size(); qhead++) {
    const Lit unLit = trail[qhead].neg();
    const int32_t lev = var(unLit).decision_level;
    VERBOSE_PRINT("&&Propagating: " << unLit.neg() << " qhead: " << qhead << " lev: " << lev);

    //Propagate bin clauses
    for (const auto& bincl : litWatchList(unLit).binary_links_) {
      const auto& l = bincl.lit();
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
    if (!confl.isNull()) break;
  }
  VERY_SLOW_DEBUG_DO(if (confl.isNull() && !check_watchlists()) {
      print_trail(false, false);assert(false);});
  SLOW_DEBUG_DO(if (confl.isNull()) check_all_propagated());
  VERBOSE_PRINT("After propagate, qhead is: " << qhead);
  return confl.isNull();
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

bool Counter::litRedundant(Lit p, uint32_t abstract_levels) {
    VERBOSE_PRINT(__func__ << " called");

    analyze_stack.clear();
    analyze_stack.push_back(p);

    Lit* c = NULL;
    uint32_t size;
    size_t top = toClear.size();
    while (!analyze_stack.empty()) {
      VERBOSE_PRINT("At point in litRedundant: " << analyze_stack.back());
      const auto reason = var(analyze_stack.back()).ante;
      assert(reason.isAnt());  //Must have a reason
      p = analyze_stack.back();
      analyze_stack.pop_back();
      fill_cl(reason, c, size, p);

      for (uint32_t i = 1; i < size; i++) {
        VERBOSE_PRINT("at i: " << i);
        Lit p2 = c[i];
        VERBOSE_PRINT("Examining lit " << p2 << " seen: " << seen[p2.var()]);
        if (!seen[p2.var()] && var(p2).decision_level > 0) {
          if (var(p2).ante.isAnt()
              && (abstractLevel(p2.var()) & abstract_levels) != 0
          ) {
              VERBOSE_PRINT("lit " << p2 << " OK");
              seen[p2.var()] = 1;
              analyze_stack.push_back(p2);
              toClear.push_back(p2.var());
          } else {
              VERBOSE_PRINT("lit " << p2 << " NOT OK");
              //Return to where we started before function executed
              for (size_t j = top; j < toClear.size(); j++) seen[toClear[j]] = 0;
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
    if (var(uip_clause[i]).ante.isNull()
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
  for(const auto& c: toClear) seen[c] = 0;
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

void Counter::vivify_cls(vector<ClauseOfs>& cls) {
  stats.vivif_tried++;
  uint32_t j = 0;
  for(uint32_t i = 0; i < cls.size(); i++) {
    bool rem = false;
    auto& off = cls[i];
    if (v_tout > 0) {
      Clause& cl = *alloc->ptr(off);
      if (cl.vivifed == 0 && (
          !cl.red || (cl.red &&
            (cl.lbd <= lbd_cutoff || (cl.used && cl.total_used > 50)))))
        rem = vivify_cl(off);
    }
    if (!rem) cls[j++] = off;
  }

  // We didn't timeout, reset vivified flag.
  if (v_tout > 0) {
    for(const auto& off: cls) {
      alloc->ptr(off)->vivifed = 0;
    }
  }
  cls.resize(j);
}

void Counter::vivify_clauses() {
  if (last_confl_vivif + conf.vivif_every > stats.conflicts) return;
  vivif_g.seed(mtrand.randInt());
  double myTime = cpuTime();
  uint64_t last_vivif_lit_rem = stats.vivif_lit_rem;
  uint64_t last_vivif_cl_minim = stats.vivif_cl_minim;
  auto last_vivif_cl_tried = stats.vivif_tried_cl;

  // Sanity check here.
  last_confl_vivif = stats.conflicts;

  // Backup
  ws_pos.clear();
  for(const auto& off: longIrredCls) {
    const Clause& cl = *alloc->ptr(off);
    ws_pos[off] = std::make_pair(cl[0], cl[1]);
  }
  for(const auto& off: longRedCls) {
    const Clause& cl = *alloc->ptr(off);
    ws_pos[off] = std::make_pair(cl[0], cl[1]);
  }

  // Set ourselves up.
  v_lev = 0;
  v_levs.clear();
  v_levs.resize(nVars()+1, -1);
  v_values.clear();
  v_values.resize(nVars()+1, X_TRI);
  v_qhead = 0;

  // Set units up
  v_trail.clear();
  for(const auto& l: trail) {
    if (var(l).decision_level == 0) v_enqueue(l);
  }
  for(const auto& l: unit_clauses_) if (v_val(l) == X_TRI) v_enqueue(l);
  bool ret = v_propagate();
  assert(ret == true);
  verb_print(2, "[vivif] setup. T: " << (cpuTime()-myTime));

  // Vivify clauses
  v_tout = conf.vivif_mult*2LL*1000LL*1000LL;
  if (stats.vivif_tried % 3 == 0) vivify_cls(longIrredCls);
  bool tout_irred = (v_tout <= 0);
  verb_print(2, "[vivif] irred vivif remain: " << v_tout/1000 << "K T: " << (cpuTime()-myTime));

  v_tout = conf.vivif_mult*5LL*1000LL*1000LL;
  vivify_cls(longRedCls);
  verb_print(2, "[vivif] red vivif remain: " << v_tout/1000 << "K T: " << (cpuTime()-myTime));
  bool tout_red = (v_tout <= 0);

  // Restore
  for(auto& ws: watches_) ws.watch_list_.clear();
  for(const auto& off: longIrredCls) v_cl_repair(off);
  for(const auto& off: longRedCls) v_cl_repair(off);
  ws_pos.clear();
  verb_print(1, "vivif finished."
      << " cl tried: " << (stats.vivif_tried_cl - last_vivif_cl_tried)
      << " cl minim: " << (stats.vivif_cl_minim - last_vivif_cl_minim)
      << " lit rem: " << (stats.vivif_lit_rem - last_vivif_lit_rem)
      << " tout-irred: " << (int)tout_irred
      << " tout-red: " << (int)tout_red
      << " T: " << (cpuTime()-myTime));
}

void Counter::v_cl_repair(ClauseOfs off) {
  Clause& cl = *alloc->ptr(off);
  auto& offs = ws_pos[off];

  auto at = std::find(cl.begin(), cl.end(), offs.first);
  assert(at != cl.end());
  std::swap(cl[0], *at);

  at = std::find(cl.begin(), cl.end(), offs.second);
  assert(at != cl.end());
  std::swap(cl[1], *at);

  uint32_t val_f = 0;
  uint32_t val_u = 0;
  uint32_t val_t = 0;
  int32_t val_t_at = -1;
  int32_t t_dec_lev = -1;
  int32_t any_val_t_pos = -1;
  int32_t mindec_12 = std::min(var(cl[0]).decision_level, var(cl[1]).decision_level);
  for(uint32_t i = 0; i < cl.size(); i++) {
    const Lit l = cl[i];
    if (val(l) == T_TRI && var(l).decision_level <= mindec_12) {
      if (val_t_at == -1) {val_t_at = i;t_dec_lev = var(l).decision_level;}
      else if (t_dec_lev > var(l).decision_level) {
        val_t_at = i;t_dec_lev = var(l).decision_level;}
    }
    if (val(l) == T_TRI) {any_val_t_pos = i;}
    if (val(l) == T_TRI) {val_t++;}
    if (val(l) == F_TRI) {val_f++;}
    if (val(l) == X_TRI) {val_u++;}
  }

  // Not conflicting
  assert(!(val_u == 0 && val_t == 0));
  // Not propagating
  assert(!(val_u == 1 && val_t == 0));

  if (val_t_at != -1) {
    litWatchList(cl[0]).addWatchLinkTo(off, cl[val_t_at]);
    litWatchList(cl[1]).addWatchLinkTo(off, cl[val_t_at]);
    VERBOSE_PRINT("Vivified cl off: " << off);
    VERBOSE_DEBUG_DO(print_cl(cl));
    return;
  }

  // We removed the TRUE
  std::sort(cl.begin(), cl.end(),
    [=](const Lit& a, const Lit& b) {
      if (val(a) == X_TRI && val(b) != X_TRI) return true;
      if (val(b) == X_TRI && val(a) != X_TRI) return false;
      if (var(a).decision_level == var(b).decision_level) {
        if(val(a) != val(b)) return val(a) == T_TRI;
        return false;
      }
      return var(a).decision_level > var(b).decision_level;
    });

  VERBOSE_PRINT("Vivified cl off: " << off);
  VERBOSE_DEBUG_DO(print_cl(cl));
  int32_t pos = (any_val_t_pos == -1) ? cl.sz/2 : any_val_t_pos;
  litWatchList(cl[0]).addWatchLinkTo(off, cl[pos]);
  litWatchList(cl[1]).addWatchLinkTo(off, cl[pos]);
}

// We could have removed a TRUE. This may be an issue.
void Counter::v_fix_watch(Clause& cl, uint32_t i) {
  if (val(cl[i]) == X_TRI || val(cl[i]) == T_TRI) return;
  auto off = alloc->get_offset(&cl);
  litWatchList(cl[i]).removeWatchLinkTo(off);
  uint32_t i2 = 2;
  for(; i2 < cl.size(); i2++) if (val(cl[i2]) == X_TRI || val(cl[i2]) == T_TRI) break;
  /* print_cl(cl); */
  assert(i2 != cl.size());
  std::swap(cl[i], cl[i2]);
  litWatchList(cl[i]).addWatchLinkTo(off, cl[cl.sz/2]);
}

void Counter::v_new_lev() {
  assert(v_lev == 0);
  v_lev++;
  v_backtrack_to = v_trail.size();
}

void Counter::v_unset(const Lit l) {
  VERBOSE_PRINT("v-unset: " << l);
  assert(v_levs[l.var()] == 1);
  v_levs[l.var()] = -1;
  v_values[l] = X_TRI;
  v_values[l.neg()] = X_TRI;
}

void Counter::v_backtrack() {
  assert(v_lev == 1);
  for(uint32_t i = v_backtrack_to; i < v_trail.size(); i++) {
    const auto& l = v_trail[i];
    v_unset(l);
  }
  v_trail.resize(v_backtrack_to);
  v_lev = 0;
  v_qhead = v_trail.size();
}

template<class T> bool Counter::v_clause_satisfied(const T& cl) const {
  for(const auto&l : cl) {
    if (v_val(l) == T_TRI) return true;
  }
  return false;
}

template<class T> bool Counter::should_have_propagated_earlier(const T& cl) const {
  uint32_t num_t = 0;
  int32_t t_lev = -1;
  int32_t maxlev_f = -1;
  for(const auto&l: cl) {
    if (val(l) == T_TRI) {
      num_t++;
      if (num_t >= 2) return false;
      t_lev = var(l).decision_level;
    }
    if (val(l) == X_TRI) return false;
    if (val(l) == F_TRI) {
      maxlev_f = std::max(maxlev_f, var(l).decision_level);
    }
  }

  // Should have propagated at level maxlev_f -- but it only got set TRUE at t_lev!
  if (maxlev_f < t_lev) return true;
  return false;

}

// Returns TRUE if we can remove the clause
bool Counter::vivify_cl(const ClauseOfs off) {
  SLOW_DEBUG_DO(for(auto& l: seen) assert(l == 0));
  bool fun_ret = false;
  Clause& cl = *alloc->ptr(off);
  cl.vivifed = 1;
  stats.vivif_tried_cl++;

  /* cout << "orig CL: " << endl; v_print_cl(cl); */
  auto it = ws_pos.find(off);
  v_new_lev();
  v_tmp.clear();
  v_tmp2.clear();
  for(const auto&l: cl) v_tmp2.push_back(l);
  std::shuffle(v_tmp2.begin(), v_tmp2.end(), vivif_g);

  // Swap to 1st & 2nd the two original 1st & 2nd
  auto sw = std::find(v_tmp2.begin(), v_tmp2.end(), it->second.first);
  std::swap(*sw, v_tmp2[0]);
  sw = std::find(v_tmp2.begin(), v_tmp2.end(), it->second.second);
  std::swap(*sw, v_tmp2[1]);

  VERBOSE_PRINT("vivifying cl offs: " << off);
  VERBOSE_DEBUG_DO(print_cl(cl));
  for(uint32_t i = 0; i < v_tmp2.size(); i++) {
    const auto& l = v_tmp2[i];
    VERBOSE_PRINT("Vivif lit l: " << l << " val: " << val_str(v_val(l)));
    if (v_val(l) == T_TRI) {v_tmp.push_back(l);break;}
    if (v_val(l) == F_TRI) continue;
    v_tmp.push_back(l);
    v_enqueue(l.neg());
    bool ret = v_propagate();
    if (!ret) {
      VERBOSE_PRINT("vivif ret FALSE, exiting");
      break;
    }
  }
  v_backtrack();
  VERBOSE_DEBUG_DO(cout << "new vivified CL offs: " << off << endl; print_cl(v_tmp));
  uip_clause.clear();
  check_implied(v_tmp);
  for(const auto&l: v_tmp) seen[l.raw()] = 1;

  uint32_t removable = 0;
  assert(it != ws_pos.end());
  for(uint32_t i = 0; i < cl.sz; i ++) {
    const Lit l = cl[i];
    // watch 0 & 1 are never removable
    if (l != it->second.first && l != it->second.second) {
      seen[l.raw()] ^= 1;
      if (seen[l.raw()]) {
        removable++;
        toClear.push_back(l.raw());
      }
    } else {
      seen[l.raw()] = 0;
    }
  }
  if (removable != 0 &&
      // TODO once chronological backtracking works, we can have level-0 stuff. Not now.
      //      so we must skip this
      !propagating_cl(v_tmp) && !conflicting_cl(v_tmp) &&
      !should_have_propagated_earlier(v_tmp)) {
    litWatchList(cl[0]).removeWatchLinkTo(off);
    litWatchList(cl[1]).removeWatchLinkTo(off);
    VERBOSE_DEBUG_DO(cout << "orig CL: " << endl; v_print_cl(cl));
    stats.vivif_cl_minim++;
    stats.vivif_lit_rem += removable;
    auto it2 = std::remove_if(cl.begin(), cl.end(),
        [=](Lit l) -> bool { return seen[l.raw()] == 1; });

    cl.resize(it2-cl.begin());
    assert(cl.sz >= 2);
    VERBOSE_DEBUG_DO(cout << "vivified CL: " << endl; v_print_cl(cl));

    std::sort(cl.begin(), cl.end(), [=](const Lit l1, const Lit l2) {
        if (v_val(l1) != v_val(l2)) {
          if (v_val(l1) == X_TRI) return true;
          return false;
        }
        return false;
      });
    if (cl.sz == 2) {
      // Not propagating
      assert(!(val(cl[0]) == X_TRI && val(cl[1])==F_TRI));
      assert(!(val(cl[1]) == X_TRI && val(cl[0])==F_TRI));
      // Not conflicting
      assert(!(val(cl[0]) == F_TRI && val(cl[1])==F_TRI));

      add_bin_cl(cl[0], cl[1], cl.red);
      if (v_val(cl[0]) == X_TRI && v_val(cl[1]) == F_TRI) {
        v_enqueue(cl[0]);
      }
      for(uint32_t v = 1; v < variables_.size(); v++) {
        auto& vdat= variables_[v];
        if (vdat.ante.isAClause() && vdat.ante.asCl() == off) {
          assert(v == cl[0].var() || v == cl[1].var());
          Lit otherLit = (v == cl[0].var()) ? cl[1] : cl[0];
          vdat.ante = Antecedent(otherLit);
        }
      }
      markClauseDeleted(off);
      fun_ret = true;
    } else {
      litWatchList(cl[0]).addWatchLinkTo(off, cl[cl.sz/2]);
      litWatchList(cl[1]).addWatchLinkTo(off, cl[cl.sz/2]);
      if (!v_clause_satisfied(cl) && v_val(cl[0]) == X_TRI && v_val(cl[1]) != X_TRI) {
        assert(v_val(cl[1]) == F_TRI);
        v_enqueue(cl[0]);
      }
      cl.update_lbd(cl.sz); // we may be smaller than LBD
    }
    bool ret = v_propagate();
    assert(ret);
  } else {
    VERBOSE_DEBUG_DO(cout << "Can't vivify." << endl);
  }

  for(const auto& l: toClear) seen[l] = 0;
  toClear.clear();
  return fun_ret;
}

TriValue Counter::v_val(const Lit l) const {
  return v_values[l];
}

void Counter::v_enqueue(const Lit l) {
  VERBOSE_PRINT("v-enq: " << l << " lev: " << v_lev);
  assert(v_val(l) == X_TRI);
  v_levs[l.var()] = v_lev;
  v_trail.push_back(l);
  v_values[l] = T_TRI;
  v_values[l.neg()] = F_TRI;
}

bool Counter::v_propagate() {
  bool ret = true;
  for (; v_qhead < v_trail.size(); v_qhead++) {
    const Lit unLit = v_trail[v_qhead].neg();

    //Propagate bin clauses
    const auto& wsbin = litWatchList(unLit).binary_links_;
    for (const auto& bincl : wsbin) {
      const auto& l = bincl.lit();
      if (v_val(l) == F_TRI) {
        VERBOSE_PRINT("Conflict from bin.");
        return false;
      } else if (v_val(l) == X_TRI) {
        v_enqueue(l);
        VERBOSE_PRINT("Bin prop: " << l);
      }
    }

    //Propagate long clauses
    auto& ws = litWatchList(unLit).watch_list_;
    v_tout-=ws.size();

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
      if (v_val(it->blckLit) == T_TRI) { *it2++ = *it; continue; }

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
      if (v_val(c[0]) == T_TRI) {
        *it2++ = ClOffsBlckL(ofs, c[0]);
        continue;
      }

      uint32_t i = 2;
      for(; i < c.sz; i++) if (v_val(c[i]) != F_TRI) break;
      // either we found a free or satisfied lit
      if (i != c.sz) {
        c[1] = c[i];
        c[i] = unLit;
        VERBOSE_PRINT("New watch for cl: " << c[1]);
        litWatchList(c[1]).addWatchLinkTo(ofs, c[0]);
      } else {
        *it2++ = *it;
        if (v_val(c[0]) == F_TRI) {
          VERBOSE_PRINT("Conflicting state from norm cl offs: " << ofs);
          ret = false;
          it++;
          break;
        } else {
          assert(v_val(c[0]) == X_TRI);
          VERBOSE_PRINT("prop long");
          v_enqueue(c[0]);
        }
      }
    }
    while(it != ws.end()) *it2++ = *it++;
    ws.resize(it2-ws.begin());
    if (!ret) break;
  }
  VERBOSE_PRINT("After propagate, v_qhead is: " << v_qhead << " returning: " << ret);
  return ret;
}

void Counter::create_fake(Lit p, uint32_t& size, Lit*& c) const
{
    tmpLit.clear();
    tmpLit.push_back(p);
    for(int32_t i = var(p).decision_level; i > 0; i--) {
      auto const& d = decision_stack_[i];
      tmpLit.push_back(Lit(d.var, val(d.var) == F_TRI));
    }
    size = tmpLit.size();
    c = tmpLit.data();

#ifdef VERBOSE_DEBUG
    cout << "Fake cl: " << endl;
    for(const auto& l: tmpLit) {
      cout << std::setw(5) << l<< " lev: " << std::setw(3) << var(l).decision_level
        << " ante: " << std::setw(8) << var(l).ante
        << " val : " << std::setw(7) << lit_val_str(l)
        << endl;
    }
#endif
}

void Counter::fill_cl(const Antecedent& ante, Lit*& c, uint32_t& size, Lit p) const {
  if (ante.isAClause()) {
    Clause* cl = alloc->ptr(ante.asCl());
    c = cl->getData();
    size = cl->sz;
  } else if (ante.isFake()) {
    create_fake(p, size, c);
  } else if (ante.isALit()) {
    //Binary
    tmpLit.resize(2);
    c = tmpLit.data();
    if (p == NOT_A_LIT) c[0] = conflLit;
    else c[0] = p;
    c[1] = ante.asLit();
    size = 2;
  } else {assert(false);}
}

int32_t Counter::get_confl_maxlev(const Lit p) const {
  Lit* c;
  uint32_t size = 0;
  fill_cl(confl, c, size, p);

  int32_t maxlev = -1;
  for(uint32_t i = 0; i < size; i ++) {
#ifdef VERBOSE_DEBUG
    cout << "confl cl[" << std::setw(5) << i << "]"
        << " lit: " << std::setw(5) << c[i]
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

void Counter::recordLastUIPCauses() {
  assert(toClear.empty());

  uip_clause.clear();
  uip_clause.push_back(Lit(0, false));
  Lit p = NOT_A_LIT;

  SLOW_DEBUG_DO(for(const auto& t:seen) assert(t == 0););
  int32_t DL = var(top_dec_lit()).decision_level;
  VERBOSE_DEBUG_DO(cout << "orig DL: " << decision_stack_.get_decision_level() << endl);
  VERBOSE_DEBUG_DO(cout << "new DL : " << DL << endl);
  VERBOSE_DEBUG_DO(print_dec_info());
  int32_t maxlev = get_confl_maxlev(p);
  go_back_to(maxlev);
  DL = var(top_dec_lit()).decision_level;

  Lit* c;
  uint32_t size;
  VERBOSE_DEBUG_DO(cout << "Doing loop:" << endl);
  int32_t index = trail.size()-1;
  uint32_t pathC = 0;
  do {
    fill_cl(confl, c, size, p);
    if (confl.isAClause()) {
      Clause& cl = *alloc->ptr(confl.asCl());
      if (cl.red && cl.lbd > lbd_cutoff) {
        cl.increaseScore();
        cl.update_lbd(calc_lbd(cl));
      }
      if (p == NOT_A_LIT) std::swap(c[0], c[1]);
    } else if (confl.isALit()) {
      if (p == NOT_A_LIT && var(c[0]).decision_level < var(c[1]).decision_level)
        std::swap(c[0], c[1]);
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
      if (!seen[q.var()] && var(q).decision_level > 0){
        increaseActivity(q);
        seen[q.var()] = 1;
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
      while (!seen[trail[index--].var()]) { SLOW_DEBUG_DO(assert(index >= 0));};
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
    seen[p.var()] = 0;
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
  SLOW_DEBUG_DO(for(const auto& s: seen) assert(s == 0));
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
