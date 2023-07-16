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

#include "counter.hpp"

#include <algorithm>
#include <complex>
#include <ios>
#include <iomanip>
#include <numeric>
#include <utility>
#include "common.hpp"
#include "comp_types/comp.hpp"
#include "cryptominisat5/solvertypesmini.h"
#include "primitive_types.hpp"
#include "counter_config.hpp"
#include "stack.hpp"
#include "structures.hpp"
#include "time_mem.hpp"
#include "IFlowCutter.hpp"
#include "graph.hpp"

void Counter::simplePreProcess()
{
  for (auto lit : unit_clauses_) {
    assert(!existsUnitClauseOf(lit.neg()) && "Formula is not UNSAT, we ran CMS before");
    if (val(lit) == X_TRI) setLiteral(lit, 0);
    assert(val(lit) == T_TRI);
  }

  bool succeeded = propagate();
  release_assert(succeeded && "We ran CMS before, so it cannot be UNSAT");
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

  if (conf.verb) stats.printShortFormulaInfo(this);
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

void Counter::td_decompose() {
  bool conditionOnCNF = indep_support_end > 3 && nVars() > 20 && nVars() <= conf.td_varlim;
  if (!conditionOnCNF) {
    verb_print(1, "skipping TD, too many/few vars. Setting branch to fallback");
    conf.branch_type = conf.branch_fallback_type;
    return;
  }

  Graph primal(nVars()+1);
  all_lits(i) {
    Lit l(i/2, i%2);
    for(const auto& l2: watches_[l].binary_links_) {
      if ((!l2.red() || (l2.red() && conf.td_with_red_bins))
          && l < l2.lit()) {
        debug_print("v1: " << l.var());
        debug_print("v2: " << l2.lit().var());
        primal.addEdge(l.var(), l2.lit().var());
      }
    }
  }

  for(const auto& off: longIrredCls) {
    Clause& cl = *alloc->ptr(off);
    for(uint32_t i = 0; i < cl.sz; i++) {
      for(uint32_t i2 = i+1; i2 < cl.sz; i2++) {
        debug_print("v1: " << cl[i].var());
        debug_print("v2: " << cl[i2].var());
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

  // run FlowCutter
  verb_print(1, "FlowCutter is running...");
  IFlowCutter FC(nVars()+1, primal.numEdges(), conf.verb);
  FC.importGraph(primal);
  TreeDecomposition td = FC.constructTD();

  td.centroid(nVars()+1, conf.verb);
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

// Self-check count without restart with CMS only
mpz_class Counter::check_norestart_cms(const Cube& c) {
  cout << "Checking count with ourselves (no verb, no restart)" << endl;
  vector<Lit> tmp;
  CMSat::SATSolver test_solver;
  test_solver.new_vars(nVars());
  // Long cls
  for(const auto& off: longIrredCls) {
    const Clause& cl = *alloc->ptr(off);
    tmp.clear();
    for(const auto& l: cl) tmp.push_back(l);
    test_solver.add_clause(ganak_to_cms_cl(tmp));
  }
  // Bin cls
  all_lits(i) {
    Lit l(i/2, i%2);
    for(const auto& l2: watches_[l].binary_links_) {
      if (l2.irred() && l < l2.lit()) {
        tmp.clear();
        tmp.push_back(l);
        tmp.push_back(l2.lit());
        test_solver.add_clause(ganak_to_cms_cl(tmp));
      }
    }
  }
  // Unit cls
  for(const auto& l: unit_clauses_) {
    tmp.clear();
    tmp.push_back(l);
    test_solver.add_clause(ganak_to_cms_cl(tmp));
  }
  // The cube
  for(const auto&l: c.cnf) {
    tmp.clear();
    tmp.push_back(l.neg());
    test_solver.add_clause(ganak_to_cms_cl(tmp));
  }
  uint64_t num = 0;
  for(;; num++) {
    auto ret = test_solver.solve();
    if (ret == CMSat::l_False) break;
    vector<CMSat::Lit> ban;
    for(uint32_t i = 0; i < test_solver.nVars(); i++) {
      ban.push_back(CMSat::Lit(i, test_solver.get_model()[i] == CMSat::l_True));
    }
    test_solver.add_clause(ban);
  }
  return num;
}

// Self-check count without restart
mpz_class Counter::check_norestart(const Cube& c) {
  cout << "Checking count with ourselves (no verb, no restart)" << endl;
  CounterConfiguration conf2 = conf;
  conf2.do_restart = 0;
  conf2.verb = 0;
  vector<Lit> tmp;
  Counter test_cnt(conf2);
  test_cnt.new_vars(nVars());
  CMSat::SATSolver test_solver;
  test_solver.new_vars(nVars());
  // Long cls
  for(const auto& off: longIrredCls) {
    const Clause& cl = *alloc->ptr(off);
    tmp.clear();
    for(const auto& l: cl) tmp.push_back(l);
    test_cnt.add_irred_cl(tmp);
    test_solver.add_clause(ganak_to_cms_cl(tmp));
  }
  // Bin cls
  all_lits(i) {
    Lit l(i/2, i%2);
    for(const auto& l2: watches_[l].binary_links_) {
      if (l2.irred() && l < l2.lit()) {
        tmp.clear();
        tmp.push_back(l);
        tmp.push_back(l2.lit());
        test_cnt.add_irred_cl(tmp);
        test_solver.add_clause(ganak_to_cms_cl(tmp));
      }
    }
  }
  // Unit cls
  for(const auto& l: unit_clauses_) {
    tmp.clear();
    tmp.push_back(l);
    test_cnt.add_irred_cl(tmp);
    test_solver.add_clause(ganak_to_cms_cl(tmp));
  }
  // The cube
  for(const auto&l: c.cnf) {
    tmp.clear();
    tmp.push_back(l.neg());
    test_cnt.add_irred_cl(tmp);
    test_solver.add_clause(ganak_to_cms_cl(tmp));
  }
  test_cnt.end_irred_cls();
  vector<Cube> ret;
  return test_cnt.outer_count(&test_solver);
}

void Counter::disable_smaller_cube_if_overlap(uint32_t i, uint32_t i2, vector<Cube>& cubes) {
  if (cubes[i].cnf.size() < cubes[i2].cnf.size()) std::swap(i, i2);
  auto c1 = ganak_to_cms_cl(cubes[i].cnf);
  auto c2 = ganak_to_cms_cl(cubes[i2].cnf);
  std::set<CMSat::Lit> assumps;
  bool unsat = false;
  for(const auto& l: c1) {
    if (assumps.count(~l)) continue;
    if (assumps.count(l)) {unsat = true; break;}
    assumps.insert(~l);
  }
  auto sz = assumps.size();
  if (unsat) return;
  for(const auto& l: c2) {
    if (assumps.count(~l)) continue;
    if (assumps.count(l)) {unsat = true; break;}
    assumps.insert(~l);
  }
  if (unsat) return;

  bool overlap = false;
  if (assumps.size() == sz) {
    // We didn't add anything to it, so it's definitely SAT
    overlap = true;
  }
  if (!overlap) {
    vector<CMSat::Lit> ass;
    ass.insert(ass.begin(), assumps.begin(), assumps.end());
    auto ret = sat_solver->solve(&ass);
    if (ret != CMSat::l_False) overlap = true;
  }

  if (overlap) {
    cout << "Two cubes overlap." << endl;
    cout << "c1: " << c1 << endl;
    cout << "c2: " << c2   << endl;
    uint32_t to_disable = cubes[i].val > cubes[i2].val ? i2 : i;
    cubes[to_disable].enabled = false;
    cout << "Disabled cube " << cubes[to_disable] << endl;
  }
}

void Counter::print_and_check_cubes(vector<Cube>& cubes) {
  verb_print(1, "Num restarts: " << stats.num_restarts);
  verb_print(2, "cubes     : ");
  for(const auto&c: cubes) verb_print(2, "-> " << c);
#ifdef SLOW_DEBUG
  for(const auto& c: cubes) {
    auto check_cnt = check_norestart(c);
    cout << "check cube: " << c << " check_cnt: " << check_cnt << endl;
    assert(check_cnt == c.val);
  }
#endif
  cout << "Total num cubes: " << cubes.size() << endl;
}

void Counter::disable_cubes_if_overlap(vector<Cube>& cubes) {
  for(uint32_t i = 0; i < cubes.size(); i++) {
    if (!cubes[i].enabled) continue;
    for(uint32_t i2 = i+1; i2 < cubes.size(); i2++) {
      if (!cubes[i2].enabled) continue;
      disable_smaller_cube_if_overlap(i, i2, cubes);
    }
  }
}

mpz_class Counter::outer_count(CMSat::SATSolver* _sat_solver) {
  mpz_class val = 0;
  sat_solver = _sat_solver;

  auto ret = sat_solver->solve();
  start_time = cpuTime();
  while(ret == CMSat::l_True) {
    vector<Cube> cubes;
    count(cubes);
    print_and_check_cubes(cubes);
    disable_cubes_if_overlap(cubes);

    // Add cubes to count, cubes & CMS
    for(const auto&c: cubes) {
      if (!c.enabled) continue;
      val+=c.val;
      sat_solver->add_clause(ganak_to_cms_cl(c.cnf));
    }
    ret = sat_solver->solve();
    if (ret == CMSat::l_False) break;

    // Add cubes to counter
    for(auto it = cubes.rbegin(); it != cubes.rend(); it++) if (it->enabled) {
      vivify_cl_toplevel(it->cnf);
      add_irred_cl(it->cnf);
    }
    decision_stack_.clear();
    if (stats.num_restarts %2 == 0) {
      vivify_clauses(true, true);
      subsume_all();
      toplevel_full_probe();
    }
    end_irred_cls();
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
  if (exit_state == RESTART) {
    ret_cubes = mini_cubes;
  } else {
    if (conf.verb) stats.printShort(this, &comp_manager_->get_cache());
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
      << " (left: " << decision_stack_.at(dec_lev).get_left_model_count()
      << " right: " << decision_stack_.at(dec_lev).get_right_model_count() << ")"
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

bool Counter::chrono_work() {
  debug_print("--- CHRONO CHECK ----");
  VERBOSE_DEBUG_DO(print_trail());
  auto data = find_conflict_level(conflLit);
  if (data.bOnlyOneLitFromHighest) {
    debug_print("ChronoBT. going back to " << data.nHighestLevel-1 << " curlev: " << decision_level());
    go_back_to(data.nHighestLevel-1);
    VERBOSE_DEBUG_DO(print_trail());
    debug_print("now Dec lev: " << decision_level());
    return true;
  }
  return false;
}

int Counter::chrono_work_sat() {
  debug_print("SAT mode chrono check");
  VERBOSE_DEBUG_DO(print_trail());
  auto data = find_conflict_level(conflLit);
  if (data.bOnlyOneLitFromHighest) {
    debug_print("SAT mode ChronoBT. going back to " << data.nHighestLevel-1 << " curlev: " << decision_level());
    if (data.nHighestLevel-1 < sat_start_dec_level) return -1;
    go_back_to(data.nHighestLevel-1);
    VERBOSE_DEBUG_DO(print_trail());
    debug_print("now Dec lev: " << decision_level());
    return 1;
  }
  return 0;
}
SOLVER_StateT Counter::countSAT() {
  retStateT state = RESOLVED;

  while (true) {
    debug_print("var top of decision stack: " << decision_stack_.top().var);
    // NOTE: findNextRemainingComponentOf finds disjoint comps
    // we then solve them all with the decideLiteral & calling findNext.. again
    while (comp_manager_->findNextRemainingComponentOf(decision_stack_.top())) {
      // It's a component. It will ONLY fall into smaller pieces if we decide on a literal
      if (!decideLiteral()) {
        decision_stack_.top().nextUnprocessedComponent();
        continue;
      }
      print_stat_line();
      if (!isindependent) {
        debug_print("before SAT mode. cnt dec: " << decision_stack_.top().getTotalModelCount()
            << " left: " << decision_stack_.top().get_left_model_count()
            << " right: " << decision_stack_.top().get_right_model_count());
        bool ret = deal_with_independent();
        debug_print("after SAT mode. cnt dec: " << decision_stack_.top().getTotalModelCount()
            << " left: " << decision_stack_.top().get_left_model_count()
            << " right: " << decision_stack_.top().get_right_model_count());
        decision_stack_.push_back(StackLevel( decision_stack_.top().currentRemainingComponent(),
              comp_manager_->comp_stack_size()));
        if (ret) {
          decision_stack_.top().change_to_right_branch();
          decision_stack_.top().includeSolution(1);
          decision_stack_.top().var = 0;
        } else {
          decision_stack_.top().branch_found_unsat();
          decision_stack_.top().change_to_right_branch();
          decision_stack_.top().branch_found_unsat();
          bool ret2 = propagate(true);
          if (!ret2) goto resolve;
        }
        debug_print("after SAT mode. cnt of this comp: " << decision_stack_.top().getTotalModelCount()
          << " unproc comps end: " << decision_stack_.top().getUnprocessedComponentsEnd()
          << " remaining comps: " << decision_stack_.top().remaining_comps_ofs()
          << " has unproc: " << decision_stack_.top().hasUnprocessedComponents());
        assert(isindependent);
        state = BACKTRACK;
        break;
      }

      while (!propagate()) {
resolve:
        if (chrono_work()) continue;
        state = resolveConflict();
        while(state == GO_AGAIN) state = resolveConflict();
        if (state == BACKTRACK) break;
      }
      if (state == BACKTRACK) break;

      // we are in RESOLVED or PROCESS_COMPONENT state, continue.
      assert(state == PROCESS_COMPONENT || state == RESOLVED);
    }
    // we are here because there is no next component, or we had to backtrack

    state = backtrack();
    if (state == PROCESS_COMPONENT && restart_if_needed()) {return RESTART;}
    if (state == EXIT) return SUCCESS;

    while (!propagate()) {
      if (chrono_work()) continue;
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

bool Counter::get_polarity(const uint32_t v) const {
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
  VERBOSE_DEBUG_DO(print_all_levels());
  debug_print("new decision level is about to be created, lev now: " << decision_stack_.get_decision_level() << " branch: " << decision_stack_.top().is_right_branch());
  decision_stack_.push_back(
    StackLevel(decision_stack_.top().currentRemainingComponent(),
               comp_manager_->comp_stack_size()));

  // The decision literal is now ready. Deal with it.
  uint32_t v = 0;
  isindependent = true;
  if (conf.branch_type == branch_t::gpmc) v = find_best_branch_gpmc();
  else v = find_best_branch();
  if (v == 0 && perform_projected_counting) {
    decision_stack_.pop_back();
    isindependent = false;
    return true;
  }
  if (v == 0) {
    // we have set all remaining var(s) from a lower decision level.
    // so there is nothing to decide. Component has a single solution.
    debug_print("We have set ALL REMAINING VARS FROM LOWER LEVELS!!");
    decision_stack_.pop_back();
    return false;
  }
  assert(val(v) == X_TRI);

  decision_stack_.top().var = v;

  Lit lit = Lit(v, get_polarity(v));
  debug_print(COLYEL "decideLiteral() is deciding: " << lit << " dec level: "
      << decision_stack_.get_decision_level());
  setLiteral(lit, decision_level());
  stats.decisions++;
  if (stats.decisions % 128 == 0) comp_manager_->rescale_cache_scores();
  assert( decision_stack_.top().remaining_comps_ofs() <= comp_manager_->comp_stack_size());
  return true;
}

uint32_t Counter::find_best_branch_gpmc() {
  uint32_t maxv = 0;
  double max_score_a = -1;
  double max_score_f = -1;
  double max_score_td = -1;

  for (auto it = comp_manager_->getSuperComponentOf(decision_stack_.top()).varsBegin();
      *it != varsSENTINEL; it++) if (*it < indep_support_end) {
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

uint32_t Counter::find_best_branch() {
  vars_scores.clear();
  uint32_t best_var = 0;
  double best_var_score = -1;
  for (auto it = comp_manager_->getSuperComponentOf(decision_stack_.top()).varsBegin();
      *it != varsSENTINEL; it++) {
    if (val(*it) != X_TRI) continue;
    if (*it < indep_support_end) {
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
      if (*it < indep_support_end) {
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

// returns cube in `c`. Uses branch 0/1, i.e. LEFT/RIGHT branch
bool Counter::compute_cube(Cube& c, int branch) {
  assert(c.val == 0);
  assert(c.cnf.empty());
  assert(conf.do_restart);
  debug_print(COLWHT "-- " << __func__ << " BEGIN");

  c.val = decision_stack_.top().get_model_side(branch);
  debug_print("Own cnt: " << c.val);
  for(int32_t i = 0; i < decision_stack_.get_decision_level(); i++) {
    const StackLevel& dec = decision_stack_[i];
    const auto& mul = dec.getBranchSols();
    if (mul == 0) continue;
    else c.val*=mul;
  }
  debug_print("Mult cnt: " << c.val);
  if (c.val == 0) return false;

  const bool opposite_branch = branch != decision_stack_.top().is_right_branch();

  // Add decisions
  debug_print(COLWHT "Decisions in the c.cnf: ");
  for(const auto& l: trail) {
    if (!var(l).ante.isNull()) continue;
    if (var(l).decision_level == decision_stack_.get_decision_level() &&
        opposite_branch) {
      assert(l == top_dec_lit());
      c.cnf.push_back(l);
    } else {
      c.cnf.push_back(l.neg());
    }
    debug_print(l << " ");
  }

  // Get a solution
  vector<CMSat::Lit> ass;
  for(const auto&l: c.cnf) ass.push_back(CMSat::Lit(l.var()-1, l.sign()));
  auto solution = sat_solver->solve(&ass);
  debug_print("cube solution: " << solution);
  if (solution == CMSat::l_False) return false;

  // Add values for all components not yet counted
  for(int32_t i = 0; i <= decision_stack_.get_decision_level(); i++) {
    if (i == decision_stack_.get_decision_level() && opposite_branch) {
      // This has been fully counted, ALL components.
      continue;
    }
    const StackLevel& dec = decision_stack_[i];
    const auto off_start = dec.remaining_comps_ofs();
    const auto off_end = dec.getUnprocessedComponentsEnd();
    debug_print("lev: " << i << " off_start: " << off_start << " off_end: " << off_end);
    // add all but the last component (it's the one being counted lower down)
    int off_by_one = 1;
    if (i == decision_stack_.get_decision_level()) off_by_one = 0;
    for(uint32_t i2 = off_start; i2 < off_end-off_by_one; i2++) {
      const auto& comp = comp_manager_->at(i2);
      all_vars_in_comp(comp, v) {
        Lit l = Lit(*v, sat_solver->get_model()[*v-1] == CMSat::l_False);
        debug_print("Lit from comp: " << l);
        c.cnf.push_back(l);
      }
    }
  }

#ifdef VERBOSE_DEBUG
  // Show decision stack's comps
  for(int32_t i = 0; i <= decision_stack_.get_decision_level(); i++) {
    const auto& dst = decision_stack_.at(i);
    cout << COLWHT "decision_stack.at(" << i << "):"
      << " decision var: " << dst.var
      << " num unproc comps: " << dst.numUnprocessedComponents()
      << " unproc comps end: " << dst.getUnprocessedComponentsEnd()
      << " remain comps offs: " << dst.remaining_comps_ofs()
      << " total count here: " << dst.getTotalModelCount()
      << " left count here: " << dst.get_left_model_count()
      << " right count here: " << dst.get_right_model_count()
      << " branch: " << dst.is_right_branch() << endl;
    const auto off_start = dst.remaining_comps_ofs();
    const auto off_end = dst.getUnprocessedComponentsEnd();
    for(uint32_t i2 = off_start; i2 < off_end; i2++) {
      assert(i2 < comp_manager_->comp_stack_size());
      const auto& comp = comp_manager_->at(i2);
      cout << COLWHT "-> comp at: " << std::setw(3) << i2 << " ID: " << comp->id() << " -- vars : ";
      all_vars_in_comp(comp, v) cout << *v << " ";
      cout << COLDEF << endl;
    }
  }

  cout << COLORG "cube so far. Size: " << c.cnf.size() << " cube: ";
  for(const auto& l: c.cnf) cout << l << " ";
  cout << endl;
  cout << COLORG "cube's SOLE count: " << decision_stack_.top().get_model_side(branch) << endl;
  cout << COLORG "cube's RECORDED count: " << c.val << COLDEF << endl;
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

static double luby(double y, int x){
    // Find the finite subsequence that contains index 'x', and the
    // size of that subsequence:
    int size, seq;
    for (size = 1, seq = 0; size < x+1; seq++, size = 2*size+1);

    while (size-1 != x){
        size = (size-1)>>1;
        seq--;
        x = x % size;
    }

    return pow(y, seq);
}

bool Counter::restart_if_needed() {
  cache_miss_rate_q.push(stats.cache_miss_rate());
  depth_q.push(decision_stack_.size());
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

  if (conf.restart_type == 6 &&
      (stats.conflicts-stats.last_restart_num_conflicts) > conf.next_restart)
    restart = true;

  if (conf.restart_type == 7 &&
      (stats.conflicts-stats.last_restart_num_conflicts) >
        (luby(2, stats.num_restarts) * conf.first_restart))
    restart = true;

  if (conf.restart_type == 8 &&
      (stats.num_cached_comps_) > (1000*luby(2, stats.num_restarts) * conf.first_restart))
    restart = true;

  if (conf.restart_type == 4 && stats.cache_hits_misses_q.isvalid()
      && stats.cache_hits_misses_q.avg() <
      stats.cache_hits_misses_q.getLongtTerm().avg()*conf.restart_cutoff_mult)
      restart = true;

  if (conf.restart_type == 5 && stats.comp_size_times_depth_q.isvalid() &&
        stats.comp_size_times_depth_q.avg() >
          stats.comp_size_times_depth_q.getLongtTerm().avg()*(1.0/conf.restart_cutoff_mult))
      restart = true;

  if (!restart) return false;
  verb_print(1, "c  ************* Restarting.  **************");
  print_restart_data();
  verb_print(2, "Num decisions since last restart: "
    << stats.decisions-stats.last_restart_num_decisions
    << endl
    << "c o Num conflicts since last restart: "
    << stats.conflicts-stats.last_restart_num_conflicts
    << endl
    << "c o Num cache lookups since last restart: "
    << stats.num_cache_look_ups_-stats.last_restart_num_cache_look_ups);

  // Reset stats
  depth_q.clear();
  cache_miss_rate_q.clear();
  comp_size_q.clear();
  stats.cache_hits_misses_q.clear();
  stats.comp_size_times_depth_q.clear();
  stats.last_restart_num_conflicts = stats.conflicts;
  stats.last_restart_num_decisions = stats.decisions;
  stats.last_restart_num_cache_look_ups = stats.num_cache_look_ups_;

  assert(mini_cubes.empty());
  while (decision_stack_.size() > 1) {
    verb_print(2, COLBLBACK <<  COLCYN "--> Mini cube gen. "
      << " lev: " << decision_stack_.get_decision_level()
      << " left cnt: " << decision_stack_.top().get_left_model_count()
      << " right cnt: " << decision_stack_.top().get_right_model_count()
      << COLDEF);
    for(uint32_t i = 0; i < 2; i++) {
      if (decision_stack_.top().get_model_side(i) == 0) continue;
      verb_print(2, "->> branch: " << i << " doing compute_cube...");
      Cube cube;
      if (compute_cube(cube, i)) {
        mini_cubes.push_back(cube);
      } else { verb_print(2, "->> FALSE cube. "); }
    }
    comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
    reactivate_comps_and_backtrack_trail();
    decision_stack_.pop_back();
  }
  stats.num_restarts++;
  return true;
}

// Checks one-by-one using a SAT solver
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
      for(uint32_t i = 0; i < c->nVars(); i++) {
        uint32_t var = c->varsBegin()[i];
        if (var < indep_support_end) active.insert(var);
      }
    } else {
      assert(single_var < (int)indep_support_end && "NO IDEA if this check is needed! TODO");
      active.insert(single_var);
    }

#ifdef VERBOSE_DEBUG
    cout << "active: "; for(const auto&a: active) cout << a << " "; cout << endl;
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
        // Ban solution
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
  do {
    debug_print("[indep] top count here: " << decision_stack_.top().getTotalModelCount() << " dec lev: " << decision_level());
    if (decision_stack_.top().branch_found_unsat()) {
      comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
    } else if (decision_stack_.top().anotherCompProcessible()) {
      debug_print("[indep] Processing another comp at dec lev "
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
      debug_print("[indep] We have NOT explored the right branch (isSecondBranch==false). Let's do it!"
          << " -- dec lev: " << decision_stack_.get_decision_level());
      const Lit aLit = top_dec_lit();
      assert(decision_stack_.get_decision_level() > 0);
      CHECK_COUNT_DO(check_count(true));
      SLOW_DEBUG_DO(assert(decision_stack_.top().get_right_model_count() == 0));
      decision_stack_.top().change_to_right_branch();
      // could be the flipped that's FALSEified so that would
      // mean the watchlist is not "sane". We need to propagate the flipped var and
      // then it'll be fine
      /* SLOW_DEBUG_DO(check_all_propagated_conflicted()); */
      reactivate_comps_and_backtrack_trail(false);
      bool ret = propagate(true);
      assert(ret);
      debug_print("[indep] Flipping lit to: " << aLit.neg() << " val is: " << val_to_str(val(aLit)));
      if (val(aLit.neg()) == X_TRI) {
        setLiteral(aLit.neg(), decision_stack_.get_decision_level());
        VERBOSE_DEBUG_DO(print_trail());
        debug_print(COLORGBG "[indep] Backtrack finished -- we flipped the branch");
        return RESOLVED;
      } else {
        assert(val(aLit.neg()) == F_TRI && "Cannot be TRUE because that would mean that the branch we just explored was UNSAT and we should have detected that");
        decision_stack_.top().branch_found_unsat();
        continue;
      }
    }
    debug_print(COLORGBG "[indep] We have explored BOTH branches, actually BACKTRACKING."
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
      debug_print("[indep] Backtracking from lev 0, i.e. ending");
      break;
    }

    CHECK_COUNT_DO(check_count());
    reactivate_comps_and_backtrack_trail(false);
    assert(decision_level() >= 1);
#ifdef VERBOSE_DEBUG
    const auto parent_count_before = (decision_stack_.end() - 2)->getTotalModelCount();
    const auto parent_count_before_left = (decision_stack_.end() - 2)->get_left_model_count();
    const auto parent_count_before_right = (decision_stack_.end() - 2)->get_right_model_count();
#endif
    (decision_stack_.end() - 2)->includeSolution(decision_stack_.top().getTotalModelCount());
    debug_print("[indep] Backtracking from level " << decision_stack_.get_decision_level()
        << " count here is: " << decision_stack_.top().getTotalModelCount());
    decision_stack_.pop_back();

    // var == 0 means it's coming from a fake decision due to normal SAT solving
    assert(decision_stack_.top().var == 0 || decision_stack_.top().var < indep_support_end);
    auto& dst = decision_stack_.top();
    debug_print("[indep] -> Backtracked to level " << decision_stack_.get_decision_level()
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

size_t Counter::find_backtrack_level_of_learnt() {
  assert(!uip_clause.empty());
  uint32_t max_i = 0;
  for (uint32_t i = 0; i < uip_clause.size(); i++) {
    if (var(uip_clause[i]).decision_level > var(uip_clause[max_i]).decision_level)
      max_i = i;
  }
  std::swap(uip_clause[max_i], uip_clause[0]);
  return var(uip_clause[0]).decision_level;
}

uint32_t Counter::find_lev_to_set(const int32_t backj) {
  assert(!uip_clause.empty());
  if (uip_clause.size() == 1) return 0;
  int32_t lev_to_set = 0;
  bool updated = false;
  uint32_t switch_to = 0;
  for (uint32_t i = 0; i < uip_clause.size(); i++) {
    int32_t lev = var(uip_clause[i]).decision_level;
      if (lev > lev_to_set && lev < backj) {
        lev_to_set = lev;
        updated = true;
        switch_to = i;
      }
  }
  debug_print("lev_to_set: " << lev_to_set << " backj: " << backj << " updated: " << (int)updated);
  assert(updated);
  std::swap(uip_clause[1], uip_clause[switch_to]);
  return lev_to_set;
}

void Counter::print_trail(bool check_entail, bool check_anything) const {
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
  cout << "qhead: " << qhead << endl;
  if (check_anything) check_trail(check_entail);
}

void Counter::go_back_to(int32_t backj) {
  debug_print("going back to lev: " << backj << " dec level now: " << decision_stack_.get_decision_level());
  while(decision_stack_.get_decision_level() > backj) {
    debug_print("at dec lit: " << top_dec_lit() << " lev: " << decision_level() << " cnt:" <<  decision_stack_.top().getTotalModelCount());
    VERBOSE_DEBUG_DO(print_comp_stack_info());
    if (!sat_mode()) {
      decision_stack_.top().mark_branch_unsat();
      decision_stack_.top().zero_out_all_sol(); //not sure it's needed
      comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
    }
    reactivate_comps_and_backtrack_trail(false);
    decision_stack_.pop_back();
    if (!sat_mode()) {
      decision_stack_.top().zero_out_branch_sol();
      comp_manager_->removeAllCachePollutionsOf(decision_stack_.top());
      comp_manager_->cleanRemainingComponentsOf(decision_stack_.top());
    }
    VERBOSE_DEBUG_DO(cout << "now at dec lit: " << top_dec_lit() << " lev: " << decision_level() << " cnt:" <<  decision_stack_.top().getTotalModelCount() << endl);
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
    if (check_entail && sat_solver) {
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
          cout << "Not implied by decisions above or at its level (" << this_lev << "): " << t << endl;
          entailment_fail = true;
        }
      }
    }
#endif
    if (entailment_fail) {
      cout << "Entailment fail." << endl;
      print_trail(false, false);
      print_dec_info();
    }
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

void Counter::reduceDB_if_needed() {
  if (stats.conflicts > last_reduceDB_conflicts+10000) {
    reduceDB();
    if (stats.cls_deleted_since_compaction > 30000 && alloc->consolidate(this)) {
        stats.cls_deleted_since_compaction = 0;
    }
    last_reduceDB_conflicts = stats.conflicts;
  }
}

// Returns TRUE if we would go further back
bool Counter::resolveConflict_sat() {
  assert(sat_mode());
  debug_print("SAT mode conflict resolution");
  recordLastUIPCause();
  if (uip_clause.size() == 1 && !existsUnitClauseOf(uip_clause[0]))
    unit_clauses_.push_back(uip_clause[0]);

  assert(uip_clause.front() != NOT_A_LIT);
  act_inc *= 1.0/conf.act_exp;
  VERBOSE_DEBUG_DO(print_conflict_info());

  int32_t backj = find_backtrack_level_of_learnt();
  int32_t lev_to_set = find_lev_to_set(backj);
  // NOTE TODO: we don't actually attach this UIP clause
  // that would take us back further up
  if (backj-1 < sat_start_dec_level) {
    debug_print("SAT mode backtrack would go back too far, not attaching UIP cl");
    return false;
  }

  stats.conflicts++;
  debug_print("SAT mode backj: " << backj << " lev_to_set: " << lev_to_set);
  VERBOSE_DEBUG_DO(print_trail());
  VERBOSE_DEBUG_DO(print_conflict_info());
  debug_print("SAT mode backj: " << backj << " lev_to_set: " << lev_to_set
    << " current lev: " << decision_level());
  go_back_to(backj-1);
  auto ant = addUIPConflictClause(uip_clause);
  setLiteral(uip_clause[0], lev_to_set, ant);
  VERBOSE_DEBUG_DO(print_trail());
  return true;
}

retStateT Counter::resolveConflict() {
  VERBOSE_DEBUG_DO(cout << "****** RECORD START" << endl);
  VERBOSE_DEBUG_DO(print_trail());

  recordLastUIPCause();
  if (uip_clause.size() == 1 && !existsUnitClauseOf(uip_clause[0]))
    unit_clauses_.push_back(uip_clause[0]);

  assert(uip_clause.front() != NOT_A_LIT);
  VERBOSE_DEBUG_DO(cout << "*RECORD FINISHED*" << endl);
  act_inc *= 1.0/conf.act_exp;

  reduceDB_if_needed();
  VERBOSE_DEBUG_DO(print_conflict_info());

  stats.conflicts++;
  assert(decision_stack_.top().remaining_comps_ofs() <= comp_manager_->comp_stack_size());
  decision_stack_.top().zero_out_branch_sol();
  decision_stack_.top().mark_branch_unsat();

  VERBOSE_DEBUG_DO(cout << "backwards cleaning" << endl);
  VERBOSE_DEBUG_DO(print_comp_stack_info());
  int32_t backj = find_backtrack_level_of_learnt();
  int32_t lev_to_set = find_lev_to_set(backj);

  debug_print("backj: " << backj << " lev_to_set: " << lev_to_set);
  bool flipped_declit = (
      uip_clause[0].neg().var() == decision_stack_.at(backj).var
       && lev_to_set+1 == backj);

  if (!flipped_declit) {
    debug_print("---- NOT FLIPPED DECLIT ----------");
    VERBOSE_DEBUG_DO(print_trail());
    VERBOSE_DEBUG_DO(print_conflict_info());
    debug_print("Not flipped. backj: " << backj << " lev_to_set: " << lev_to_set
      << " current lev: " << decision_level());
    go_back_to(backj-1);
    auto ant = addUIPConflictClause(uip_clause);
    setLiteral(uip_clause[0], lev_to_set, ant);
    VERBOSE_DEBUG_DO(print_trail());
    return RESOLVED;
  }

  assert(flipped_declit);
  VERBOSE_DEBUG_DO(cout << "after finding backj lev: " << backj << " lev_to_set: " << lev_to_set <<  endl);
  VERBOSE_DEBUG_DO(print_conflict_info());

  go_back_to(backj);
  VERBOSE_DEBUG_DO(print_conflict_info());
  debug_print("decision_stack_.get_decision_level(): " << decision_stack_.get_decision_level());

  Antecedent ant;
  assert(!uip_clause.empty());
  SLOW_DEBUG_DO(check_implied(uip_clause));
  if (decision_stack_.get_decision_level() > 0 && top_dec_lit().neg() == uip_clause[0]) {
    debug_print("FLIPPING. Setting reason the conflict cl");
    assert(var(uip_clause[0]).decision_level != -1);
    ant = addUIPConflictClause(uip_clause);
    var(top_dec_lit()).ante = ant;
  }
  debug_print("Ant is :" << ant);
  debug_print("AFTER conflict, setup: ");
  VERBOSE_DEBUG_DO(print_conflict_info());
  debug_print("is right here? " << decision_stack_.top().is_right_branch());

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

bool Counter::propagate(bool out_of_order) {
  confl = Antecedent();
  debug_print("qhead in propagate(): " << qhead << " trail sz: " << trail.size());
  for (; qhead < trail.size(); qhead++) {
    const Lit unLit = trail[qhead].neg();
    const int32_t lev = var(unLit).decision_level;
    bool lev_at_declev = false;
    if (!out_of_order) {
      if (decision_stack_.size() <= 1) lev_at_declev = true;
      else if (var(top_dec_lit()).decision_level == lev) lev_at_declev = true;
    }
    debug_print("&&Propagating: " << unLit.neg() << " qhead: " << qhead << " lev: " << lev);

    //Propagate bin clauses
    for (const auto& bincl : litWatchList(unLit).binary_links_) {
      const auto& l = bincl.lit();
      if (val(l) == F_TRI) {
        setConflictState(unLit, l);
        VERBOSE_DEBUG_DO(cout << "Bin confl. otherlit: " << l << endl);
      } else if (val(l) == X_TRI) {
        setLiteral(l, lev, Antecedent(unLit));
        VERBOSE_DEBUG_DO(cout << "Bin prop: " << l << " lev: " << lev << endl);
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
        debug_print("New watch for cl: " << c[1]);
        litWatchList(c[1]).addWatchLinkTo(ofs, c[0]);
      } else {
        *it2++ = *it;
        if (val(c[0]) == F_TRI) {
          debug_print("Conflicting state from norm cl offs: " << ofs);
          setConflictState(&c);
          it++;
          break;
        } else {
          assert(val(c[0]) == X_TRI);
          debug_print("prop long lev: " << lev << " dec_stack.get_lev : " << decision_stack_.get_decision_level());
          if (lev_at_declev) {
            setLiteral(c[0], lev, Antecedent(ofs));
            debug_print("Norm long prop: " << c[0] << " lev: " << lev);
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
  VERY_SLOW_DEBUG_DO(
      if (confl.isNull()) check_trail();
      if (confl.isNull() && !check_watchlists()) {
      print_trail(false, false);
      assert(false);
    });
  SLOW_DEBUG_DO(if (confl.isNull()) check_all_propagated_conflicted());
  debug_print("After propagate, qhead is: " << qhead);
  return confl.isNull();
}

const DataAndStatistics& Counter::get_stats() const
{
  return stats;
}

bool Counter::litRedundant(Lit p, uint32_t abstract_levels) {
    debug_print(__func__ << " called");

    analyze_stack.clear();
    analyze_stack.push_back(p);

    Lit* c = NULL;
    uint32_t size;
    size_t top = toClear.size();
    while (!analyze_stack.empty()) {
      debug_print("At point in litRedundant: " << analyze_stack.back());
      const auto reason = var(analyze_stack.back()).ante;
      assert(reason.isAnt());  //Must have a reason
      p = analyze_stack.back();
      analyze_stack.pop_back();
      fill_cl(reason, c, size, p);

      for (uint32_t i = 1; i < size; i++) {
        debug_print("at i: " << i);
        Lit p2 = c[i];
        debug_print("Examining lit " << p2 << " seen: " << (int)seen[p2.var()]);
        if (!seen[p2.var()] && var(p2).decision_level > 0) {
          if (var(p2).ante.isAnt()
              && (abstractLevel(p2.var()) & abstract_levels) != 0
          ) {
              debug_print("lit " << p2 << " OK");
              seen[p2.var()] = 1;
              analyze_stack.push_back(p2);
              toClear.push_back(p2.var());
          } else {
              debug_print("lit " << p2 << " NOT OK");
              //Return to where we started before function executed
              for (size_t j = top; j < toClear.size(); j++) seen[toClear[j]] = 0;
              toClear.resize(top);
              return false;
          }
        }
      }
    }
    debug_print("Returning OK from " << __func__);
    return true;
}

uint32_t Counter::abstractLevel(const uint32_t x) const
{
    return ((uint32_t)1) << (variables_[x].decision_level & 31);
}

void Counter::recursiveConfClauseMin()
{
    VERBOSE_DEBUG_DO(print_conflict_info());
  debug_print("recursive ccmin now.");
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
      debug_print("ccmin -- keeping lit: " << uip_clause[i]);
      uip_clause[j++] = uip_clause[i];
    } else {
      debug_print("ccmin -- NOT keeping lit: " << uip_clause[i]);
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

void Counter::vivify_clauses(bool force, bool only_irred) {
  if (!force && last_confl_vivif + conf.vivif_every > stats.conflicts) return;
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
  for(const auto& l: trail) if (var(l).decision_level == 0) v_enqueue(l);
  for(const auto& l: unit_clauses_) if (v_val(l) == X_TRI) v_enqueue(l);
  bool ret = v_propagate();
  assert(ret == true);
  verb_print(2, "[vivif] setup. T: " << (cpuTime()-myTime));

  // Vivify clauses
  v_tout = conf.vivif_mult*1LL*1000LL*1000LL;
  vivify_cls(longIrredCls);
  bool tout_irred = (v_tout <= 0);
  verb_print(2, "[vivif] irred vivif remain: " << v_tout/1000 << "K T: " << (cpuTime()-myTime));

  bool tout_red = false;
  if (!only_irred) {
    v_tout = conf.vivif_mult*5LL*1000LL*1000LL;
    vivify_cls(longRedCls);
    verb_print(2, "[vivif] red vivif remain: " << v_tout/1000 << "K T: " << (cpuTime()-myTime));
    tout_red = (v_tout <= 0);
  }

  // Restore
  for(auto& ws: watches_) ws.watch_list_.clear();
  if (decision_stack_.size() != 0) {
    for(const auto& off: longIrredCls) v_cl_repair(off);
    for(const auto& off: longRedCls) v_cl_repair(off);
  } else {
    // Move all 0-level stuff to unit_clauses_
    for(const auto& l: v_trail) {
      if (val(l) == X_TRI) {
        setLiteral(l, 0);
        if (!existsUnitClauseOf(l)) unit_clauses_.push_back(l);
      }
      assert(val(l) != F_TRI); // it would be UNSAT
    }
    bool ret2 = propagate();
    assert(ret2);
    v_cl_toplevel_repair(longIrredCls);
    v_cl_toplevel_repair(longRedCls);
  }
  ws_pos.clear();
  verb_print(1, "vivif finished."
      << " cl tried: " << (stats.vivif_tried_cl - last_vivif_cl_tried)
      << " cl minim: " << (stats.vivif_cl_minim - last_vivif_cl_minim)
      << " lit rem: " << (stats.vivif_lit_rem - last_vivif_lit_rem)
      << " tout-irred: " << (int)tout_irred
      << " tout-red: " << (int)tout_red
      << " T: " << (cpuTime()-myTime));
}

template<class T> bool Counter::v_satisfied(const T& lits) {
  for(auto& l: lits) if (v_val(l) == T_TRI) return true;
  return false;
}

template<class T> bool Counter::v_unsat(const T& lits) {
  for(auto& l: lits) if (v_val(l) == T_TRI || v_val(l) == X_TRI) return false;
  return true;
}

void Counter::v_shrink(Clause& cl) {
  uint32_t j = 0;
  for(uint32_t i = 0; i < cl.size(); i++) {
    if (v_val(cl[i]) == F_TRI) continue;
    cl[j++] = cl[i];
  }
  cl.resize(j);
}

void Counter::v_cl_toplevel_repair(vector<ClauseOfs>& offs) {
  uint32_t j = 0;
  for(uint32_t i = 0; i < offs.size(); i++) {
    Clause* cl = alloc->ptr(offs[i]);
    assert(!v_unsat(*cl));
    if (v_satisfied(*cl)) {alloc->clauseFree(cl);continue;}
    v_shrink(*cl);
    assert(cl->size() >= 2);
    if (cl->size() == 2) {
      add_bin_cl((*cl)[0], (*cl)[1], cl->red);
      alloc->clauseFree(cl);
      continue;
    }
    attach_cl(offs[i], (*cl));
    offs[j++] = offs[i];
  }
  offs.resize(j);
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

  /* uint32_t val_f = 0; */
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
    /* if (val(l) == F_TRI) {val_f++;} */
    if (val(l) == X_TRI) {val_u++;}
  }

  // Not conflicting
  assert(!(val_u == 0 && val_t == 0));
  // Not propagating
  assert(!(val_u == 1 && val_t == 0));

  if (val_t_at != -1) {
    litWatchList(cl[0]).addWatchLinkTo(off, cl[val_t_at]);
    litWatchList(cl[1]).addWatchLinkTo(off, cl[val_t_at]);
    debug_print("Vivified cl off: " << off);
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

  debug_print("Vivified cl off: " << off);
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
  debug_print("v-unset: " << l);
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

  debug_print("vivifying cl offs: " << off);
  VERBOSE_DEBUG_DO(print_cl(cl));
  for(uint32_t i = 0; i < v_tmp2.size(); i++) {
    const auto& l = v_tmp2[i];
    debug_print("Vivif lit l: " << l << " val: " << val_to_str(v_val(l)));
    if (v_val(l) == T_TRI) {v_tmp.push_back(l);break;}
    if (v_val(l) == F_TRI) continue;
    v_tmp.push_back(l);
    v_enqueue(l.neg());
    bool ret = v_propagate();
    if (!ret) {
      debug_print("vivif ret FALSE, exiting");
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
      alloc->clauseFree(off);
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
  debug_print("v-enq: " << l << " lev: " << v_lev);
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
        debug_print("Conflict from bin.");
        return false;
      } else if (v_val(l) == X_TRI) {
        v_enqueue(l);
        debug_print("Bin prop: " << l);
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
        debug_print("New watch for cl: " << c[1]);
        litWatchList(c[1]).addWatchLinkTo(ofs, c[0]);
      } else {
        *it2++ = *it;
        if (v_val(c[0]) == F_TRI) {
          debug_print("Conflicting state from norm cl offs: " << ofs);
          ret = false;
          it++;
          break;
        } else {
          assert(v_val(c[0]) == X_TRI);
          debug_print("prop long");
          v_enqueue(c[0]);
        }
      }
    }
    while(it != ws.end()) *it2++ = *it++;
    ws.resize(it2-ws.begin());
    if (!ret) break;
  }
  debug_print("After propagate, v_qhead is: " << v_qhead << " returning: " << ret);
  return ret;
}

void Counter::vivify_cl_toplevel(vector<Lit>& cl) {
  /* cout << "orig CL: " << endl; print_cl(cl); */
  assert(decision_level() == 0);
  decision_stack_.push_back(
  StackLevel(decision_stack_.top().currentRemainingComponent(),
             comp_manager_->comp_stack_size()));
  v_tmp2 = cl;
  cl.clear();
  std::shuffle(v_tmp2.begin(), v_tmp2.end(), vivif_g);

  for(uint32_t i = 0; i < v_tmp2.size(); i++) {
    const auto& l = v_tmp2[i];
    debug_print("Vivif lit l: " << l << " val: " << val_to_str(val(l)));
    if (val(l) == T_TRI) {cl.push_back(l);break;}
    if (val(l) == F_TRI) continue;
    cl.push_back(l);
    setLiteral(l.neg(), 1);
    bool ret = propagate();
    if (!ret) {
      debug_print("vivif ret FALSE, exiting");
      break;
    }
  }
  reactivate_comps_and_backtrack_trail();
  decision_stack_.pop_back();
  /* cout << "toplevel vivified CL: " << endl; print_cl(cl); */
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
    c = cl->data();
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

Counter::ConflictData Counter::find_conflict_level(Lit p) {
  ConflictData data;
  Lit* c;
  uint32_t size;
  fill_cl(confl, c, size, p);
  VERBOSE_DEBUG_DO(cout << "CL in find_conflict_level: " << endl;print_cl(c, size));
  data.nHighestLevel = var(c[0]).decision_level;
  if (data.nHighestLevel == decision_level() && var(c[1]).decision_level == decision_level())
    return data;

  int highestId = 0;
  data.bOnlyOneLitFromHighest = true;
  // find the largest decision level in the clause
  for (uint32_t i = 1; i < size; ++i) {
    int32_t lev = var(c[i]).decision_level;
    if (lev > data.nHighestLevel) {
      highestId = i;
      data.nHighestLevel = lev;
      data.bOnlyOneLitFromHighest = true;
    } else if (lev == data.nHighestLevel && data.bOnlyOneLitFromHighest == true) {
      data.bOnlyOneLitFromHighest = false;
    }
  }

  // fixing clause & watchlist
  if (highestId != 1 && confl.isAClause()) {
    Clause& cl = *alloc->ptr(confl.asCl());
    std::swap(cl[1], cl[highestId]); // swap to position 1, since we'll swap 1&0 in recordLastUIPClauses
    debug_print("SWAPPED");
    VERBOSE_DEBUG_DO(print_cl(cl.data(), cl.size()));
    if (highestId > 1 && size > 2) {
      ClauseOfs off = confl.asCl();
      litWatchList(cl[highestId]).removeWatchLinkTo(off);
      litWatchList(c[1]).addWatchLinkTo(off, c[0]);
    }
  }
	return data;
}

void Counter::recordLastUIPCause() {
  assert(toClear.empty());

  uip_clause.clear();
  uip_clause.push_back(Lit(0, false));
  Lit p = NOT_A_LIT;

  SLOW_DEBUG_DO(for(const auto& t:seen) assert(t == 0););
  VERBOSE_DEBUG_DO(print_dec_info());
  int32_t nDecisionLevel;

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
    SLOW_DEBUG_DO(if (p == NOT_A_LIT) check_cl_unsat(c, size));

    VERBOSE_DEBUG_DO(cout << "next cl: " << endl;print_cl(c, size));
    if (p == NOT_A_LIT) nDecisionLevel = var(c[0]).decision_level;
    VERBOSE_DEBUG_DO(cout << "nDecisionLevel: " <<  nDecisionLevel << endl);

    VERBOSE_DEBUG_DO(cout << "For loop." << endl);
    for(uint32_t j = ((p == NOT_A_LIT) ? 0 : 1); j < size ;j++) {
      Lit q = c[j];
      if (!seen[q.var()] && var(q).decision_level > 0){
        increaseActivity(q);
        seen[q.var()] = 1;
        toClear.push_back(q.var());
#ifdef VERBOSE_DEBUG
        cout << std::setw(5) << q
          << " lev: " << std::setw(3) << var(q).decision_level
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
      cout << "going back on trail: " << std::setw(5) << p
        << " lev: " << std::setw(3) << var(p).decision_level
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
  VERBOSE_DEBUG_DO(cout << "UIP cl: " << endl; print_cl(uip_clause.data(), uip_clause.size()));
  SLOW_DEBUG_DO(check_implied(uip_clause));
  minimizeUIPClause();
  SLOW_DEBUG_DO(for(const auto& s: seen) assert(s == 0));
}

Counter::Counter(const CounterConfiguration& _conf) :
    Instance(_conf), order_heap(VarOrderLt(watches_))
{
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
  all_lits(i) {
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
  all_lits(i) {
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

  // Check that all clauses are attached 2x in the watchlist
  map<ClauseOfs, uint32_t> off_att_num;
  all_lits(i) {
    Lit lit = Lit(i/2, i%2);
    for(const auto& ws: watches_[lit].watch_list_) {
      if (off_att_num.find(ws.ofs) == off_att_num.end()) off_att_num[ws.ofs] = 1;
      else off_att_num[ws.ofs]++;
    }
  }
  auto check_attach = [&](ClauseOfs off) {
    if (off_att_num.find(off) == off_att_num.end()) {
      cout << "ERROR: Not found clause in watchlist." << endl;
      print_cl(*alloc->ptr(off));
      ret = false;
    }
    if (off_att_num[off] !=2 ) {
      cout << "ERROR: Clause not attached 2 times. It's attached: " << off_att_num[off] << " times" << endl;
      print_cl(*alloc->ptr(off));
      ret = false;
    }
    off_att_num.erase(off);
  };
  for(const auto& off: longIrredCls) check_attach(off);
  for(const auto& off: longRedCls) check_attach(off);
  if (!off_att_num.empty()) {
    cout << "ERROR: The following clauses are attached but are NOT in longRed/longIrred clauses" << endl;
    for(const auto& p: off_att_num) {
      cout << "Offset: " << p.first << endl;
      print_cl(*alloc->ptr(p.first));
    }
  }
  return ret;
}

void Counter::attach_occ(vector<ClauseOfs>& cls) {
  for(const auto& off: cls) {
    clauses.push_back(off);
    Clause& cl = *alloc->ptr(off);
    std::sort(cl.begin(), cl.end());
    auto abs = calcAbstraction(cl);
    for(const auto& l: cl) occ[l.raw()].push_back(OffAbs(off, abs));
  }
  cls.clear();
}

void Counter::backw_susume_cl(ClauseOfs off) {
  Clause& cl = *alloc->ptr(off);
  uint32_t abs = calcAbstraction(cl);
  uint32_t smallest = numeric_limits<uint32_t>::max();
  uint32_t smallest_at = 0;
  for(uint32_t i = 0; i < cl.size(); i++) {
    Lit l = cl[i];
    if (occ[l.raw()].size() < smallest) {
      smallest = occ[l.raw()].size();
      smallest_at = i;
    }
  }

  for(const auto& check: occ[cl[smallest_at].raw()]) {
    if (off == check.off) continue;
    if (!subsetAbst(abs, check.abs)) continue;
    Clause& check_cl = *alloc->ptr(check.off);
    if (check_cl.freed) continue;
    if (subset(cl, check_cl)) {
      if (cl.red && !check_cl.red) cl.red = false;
      if (cl.red && check_cl.red) {
        cl.used |= check_cl.used;
        cl.total_used += check_cl.total_used;
      }
      debug_print( "Subsumed cl: " << check_cl << endl
                << "->by cl    : " << cl);
      alloc->clauseFree(&check_cl);
      stats.subsumed_cls++;
    }
  }
}

void Counter::backw_susume_cl_with_bin(BinClSub& cl) {
  uint32_t abs = calcAbstraction(cl);
  uint32_t smallest = numeric_limits<uint32_t>::max();
  uint32_t smallest_at = 0;
  for(uint32_t i = 0; i < cl.size(); i++) {
    Lit l = cl[i];
    if (occ[l.raw()].size() < smallest) {
      smallest = occ[l.raw()].size();
      smallest_at = i;
    }
  }

  for(const auto& check: occ[cl[smallest_at].raw()]) {
    if (!subsetAbst(abs, check.abs)) continue;
    Clause& check_cl = *alloc->ptr(check.off);
    if (check_cl.freed) continue;
    if (subset(cl, check_cl)) {
      if (cl.red && !check_cl.red) cl.red = false;
      debug_print( "Subsumed cl: " << check_cl << endl
                << "->by cl    : " << cl);
      alloc->clauseFree(&check_cl);
      stats.subsumed_cls++;
    }
  }
}

void Counter::toplevel_full_probe() {
  SLOW_DEBUG_DO(for(auto& l: seen) assert(l == 0));
  assert(toClear.empty());

  double myTime = cpuTime();
  auto old_probe = stats.toplevel_probe_fail;
  auto old_bprop = stats.toplevel_bothprop_fail;
  stats.toplevel_probe_runs++;
  assert(decision_stack_.size() == 0);
  // 0 dec level.
  decision_stack_.push_back(StackLevel(1,2));

  for(uint32_t i = 1; i <= nVars(); i++) {
    Lit l = Lit(i, 0);
    if (val(l) != X_TRI) continue;

    assert(decision_level() == 0);
    decision_stack_.push_back(StackLevel(1,2));
    decision_stack_.back().var = l.var();
    setLiteral(l, 1);
    uint32_t trail_before = trail.size();
    bool ret = propagate();
    if (ret) {
      for(uint32_t i2 = trail_before; i2 < trail.size(); i2++) {
        Lit l2 = trail[i2];
        seen[l2.raw()] = 1;
        toClear.push_back(l2.raw());
      }
    }
    reactivate_comps_and_backtrack_trail();
    decision_stack_.pop_back();
    if (!ret) {
      setLiteral(l.neg(), 0);
      ret = propagate();
      assert(ret && "we are never UNSAT");
      stats.toplevel_probe_fail++;
      continue;
    }

    // Negation
    assert(decision_level() == 0);
    decision_stack_.push_back(StackLevel(1,2));
    decision_stack_.back().var = l.var();
    setLiteral(l.neg(), 1);

    trail_before = trail.size();
    ret = propagate();
    if (ret) {
      for(uint32_t i2 = trail_before; i2 < trail.size(); i2++) {
        Lit l2 = trail[i2];
        if (seen[l2.raw()] == 1) {
          bothprop_toset.push_back(l2);
          stats.toplevel_bothprop_fail++;
        }
      }
    }
    reactivate_comps_and_backtrack_trail();
    decision_stack_.pop_back();
    if (!ret) {
      for(const auto& x: toClear) seen[x] = 0;
      toClear.clear();
      setLiteral(l, 0);
      ret = propagate();
      assert(ret && "we are never UNSAT");
      stats.toplevel_probe_fail++;
      continue;
    }
    for(const auto& x: toClear) seen[x] = 0;
    toClear.clear();
    for(const auto& x: bothprop_toset) setLiteral(x, 0);
    bothprop_toset.clear();
    ret = propagate();
    assert(ret && "we are never UNSAT");
  }
  decision_stack_.clear();
  verb_print(1, "toplevel "
      << " probe f: " << (old_probe - stats.toplevel_probe_fail)
      << " bprop f: " << (old_bprop - stats.toplevel_bothprop_fail)
      << " T: " << (cpuTime()-myTime));
}

void Counter::subsume_all() {
  assert(decision_stack_.size() == 0);
  assert(occ.empty());
  assert(clauses.empty());

  // setup
  double myTime = cpuTime();
  auto old_subsumed_cls = stats.subsumed_cls;
  auto old_subsumed_bin_cls = stats.subsumed_bin_cls;
  stats.subsume_runs++;
  occ.resize((nVars()+1)*2);
  attach_occ(longIrredCls);
  attach_occ(longRedCls);
  for(auto& ws: watches_) ws.watch_list_.clear();

  // Binary clauses
  vector<BinClSub> bin_cls;
  all_lits(i) {
    Lit lit = Lit(i/2, i%2);
    for(const auto& l2: watches_[lit].binary_links_) {
      if (l2.lit() < lit) continue;
      assert(lit < l2.lit());
      bin_cls.push_back(BinClSub(lit, l2.lit(), l2.red()));
    }
    watches_[lit].binary_links_.clear();
  }
  std::sort(bin_cls.begin(), bin_cls.end());
  uint32_t j = 0;
  for(uint32_t i = 1; i < bin_cls.size(); i++) {
    if (bin_cls[i] == bin_cls[j]) {stats.subsumed_bin_cls++; continue;}
    if (bin_cls[i].lit[0] == bin_cls[j].lit[0]
       && bin_cls[i].lit[1] == bin_cls[j].lit[1]) {
      // ordering ensures IRRED is first
      stats.subsumed_bin_cls++;
      continue;
    }
    j++;
    bin_cls[j] = bin_cls[i];
  }
  j++;
  if (!bin_cls.empty()) bin_cls.resize(j);
  for(auto& b: bin_cls) backw_susume_cl_with_bin(b);

  // Long clauses
  std::shuffle(clauses.begin(), clauses.end(), std::default_random_engine(mtrand.randInt()));
  for(const auto& off: clauses) {
    Clause* cl = alloc->ptr(off);
    if (cl->freed) continue;
    backw_susume_cl(off);
  }

  // Cleanup
  for(const auto& b: bin_cls) add_bin_cl(b[0], b[1], b.red);
  for(const auto& off: clauses) {
    Clause& cl = *alloc->ptr(off);
    if (cl.freed) continue;
    if (cl.red) longRedCls.push_back(off);
    else longIrredCls.push_back(off);
    attach_cl(off, cl);
  }
  occ.clear();
  clauses.clear();
  verb_print(1, "subs "
      << " bin-cls: " << stats.subsumed_bin_cls - old_subsumed_bin_cls
      << " long-cls: " << stats.subsumed_cls - old_subsumed_cls
      << " T: " << (cpuTime() - myTime))
}

// SAT or UNSAT
bool Counter::deal_with_independent() {
  assert(!isindependent);
  assert(order_heap.empty());
  assert(decision_stack_.size() > 0);
  assert(!sat_mode());
  stats.sat_called++;
  auto conflicts_before = stats.conflicts;

  debug_print("Entering SAT mode. Declev: " << decision_level());
  sat_start_dec_level = decision_level();
  bool sat = false;

  // Create dummy decision level in order for getSuperComponentOf work correctly.
  decision_stack_.push_back(StackLevel( decision_stack_.top().currentRemainingComponent(),
        comp_manager_->comp_stack_size()));
  for (auto it = comp_manager_->getSuperComponentOf(decision_stack_.top()).varsBegin();
      *it != varsSENTINEL; it++) {
    if (val(*it) != X_TRI) continue;
    if (*it < indep_support_end) {
      assert(false && "Only non-indep remains");
    } else {
      order_heap.insert(*it);
    }
  }
  decision_stack_.pop_back();

  // the SAT loop
  while(true) {
    uint32_t d;
    do {
      if (order_heap.empty()) {d = 0; break;}
      d = order_heap.removeMin();
    } while (val(d) != X_TRI);
    if (d == 0) {
      debug_print("SAT mode found a solution");
      SLOW_DEBUG_DO(check_sat_solution());
      sat = true;
      break;
    }
    assert(val(d) == X_TRI);
    Lit l(d, var(d).last_polarity);
    decision_stack_.push_back(StackLevel(1,2));
    decision_stack_.back().var = l.var();
    setLiteral(l, decision_level());
prop:
    bool ret = propagate();
    if (!ret) {
      print_stat_line();
      int x = chrono_work_sat();
      if (x == -1) {
        debug_print("SAT mode found UNSAT -- chrono BT");
        sat = false;
        break;
      }
      if (x == 1) goto prop;
      if (resolveConflict_sat()) {
        goto prop;
      } else {
        debug_print("SAT mode found UNSAT");
        sat = false;
        break;
      }
    } else {
      // TODO restart sometimes
    }
  }

  go_back_to(sat_start_dec_level);
  assert(decision_level() == sat_start_dec_level);
  order_heap.clear();
  sat_start_dec_level = -1;
  isindependent = true;
  debug_print("Exiting SAT mode. Declev: " << decision_level() << " sat: " << (int)sat);
  if (sat) stats.sat_found_sat++;
  else stats.sat_found_unsat++;
  stats.sat_conflicts += stats.conflicts-conflicts_before;

  return sat;
}

void Counter::check_sat_solution() const {
  assert(sat_mode());
  bool ok = true;

  for(const auto& off: longIrredCls) {
    Clause& cl = *alloc->ptr(off);
    if (clause_falsified(cl)) {
      ok = false;
      cout << "ERROR: SAT mode found a solution that falsifies a clause." << endl;
      print_cl(cl);
    }
  }

  for(const auto& off: longRedCls) {
    Clause& cl = *alloc->ptr(off);
    if (clause_falsified(cl)) {
      ok = false;
      cout << "ERROR: SAT mode found a solution that falsifies a clause." << endl;
      print_cl(cl);
    }
  }

  assert(ok);
}

#ifdef SLOW_DEBUG
template<class T> void Counter::check_cl_propagated_conflicted(T& cl) const {
  Lit unk = NOT_A_LIT;
  uint32_t num_unknown = 0;
  bool satisfied = false;
  for(const auto& l: cl) {
    if (isTrue(l)) {satisfied = true; break;}
    if (isUnknown(l)) {num_unknown++; unk = l;}
    if (num_unknown > 1) break;
  }

  if (!satisfied && num_unknown == 1) {
    cout << "ERROR! Clause should have propagated: " << unk << endl;
    print_cl(cl);
    assert(false);
  }
  if (!satisfied && num_unknown == 0) {
    cout << "ERROR! Clause should have conflicted" << endl;
    print_cl(cl);
    assert(false);
  }
}

void Counter::check_all_propagated_conflicted() const {
  // Everything that should have propagated, propagated
  for(const auto& cl: debug_irred_cls) check_cl_propagated_conflicted(cl);
  for(const auto& off: longIrredCls) {
    const Clause& cl = *alloc->ptr(off);
    check_cl_propagated_conflicted(cl);
  }
  for(const auto& off: longRedCls) {
    const Clause& cl = *alloc->ptr(off);
    check_cl_propagated_conflicted(cl);
  }
}
#endif
