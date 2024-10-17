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
#include <cstdint>
#include <ios>
#include <iomanip>
#include <limits>
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
#ifdef BUDDY_ENABLED
#include "bdd.h"
#endif
#include "mpreal.h"
#include "approxmc.h"
#include <thread>
#include <algorithm>

using std::setw;
using std::is_same;
using std::setprecision;
using std::setw;

template<typename T>
vector<uint32_t> Counter<T>::common_indep_code(const set<uint32_t>& indeps) {
  if (!num_vars_set) {
    cout << "ERROR: new_vars() MUST be called before setting indep support" << endl;
    exit(-1);
  }
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
      cout << "ERROR: independent support MUST start from variable 1 and be consecutive, e.g. 1,2,3,4,5. It cannot skip any variables. You skipped variable: " << i+1 << endl;
      exit(-1);
    }
  }

  return tmp;
}

template<typename T>
void Counter<T>::set_optional_indep_support(const set<uint32_t> &indeps) {
  auto tmp = common_indep_code(indeps);
  if (tmp.size() +1 < indep_support_end) {
    cout << "ERROR: The optional indeps MUST contain ALL indeps, plus the optional ones" << endl;
    assert(false);
    exit(-1);
  }
  if (tmp.empty()) { opt_indep_support_end = 0; return; }
  opt_indep_support_end = tmp.back()+1;

  verb_print(1, "opt ind size: " << std::min<int>((int)opt_indep_support_end-1, 0) << " ind size: " << std::min<int>((int)indep_support_end-1, 0)
    << " nvars: " << nVars());

  if (conf.verb) {
    cout << "c o indep/optional/none distribution: ";
    for(uint32_t i = 0; i <= nVars(); i++) {
      if (i < opt_indep_support_end) {
        if (i < indep_support_end) cout << "I";
        else cout << "O";
      } else cout << "N";
    }
    cout << endl;
  }
}

template<typename T>
void Counter<T>::set_indep_support(const set<uint32_t> &indeps) {
  opt_indep_support_end = nVars()+1;
  auto tmp = common_indep_code(indeps);
  if (tmp.empty()) {
    indep_support_end = 0;
    opt_indep_support_end = 0;
    return;
  }
  indep_support_end = tmp.back()+1;
  opt_indep_support_end = tmp.back()+1;

  verb_print(1, "ind size: " << indep_support_end-1 << " nvars: " << nVars());
}

// Returns false if the clause is auto-satisfied
template<typename T>
bool Counter<T>::remove_duplicates(vector<Lit>& lits) {
  if (lits.size() <= 1) return true;
  std::sort(lits.begin(), lits.end());
  uint32_t j = 1;
  Lit last_lit = lits[0];
  for(uint32_t i = 1; i < lits.size(); i++) {
    if (lits[i] == last_lit) continue;
    if (lits[i] == last_lit.neg()) return false;
    last_lit = lits[i];
    lits[j++] = lits[i];
  }
  lits.resize(j);
  return true;
}

template<typename T>
void Counter<T>::compute_score(TWD::TreeDecomposition& tdec, bool print) {
  const auto& bags = tdec.Bags();
  td_width = tdec.width();
  const auto& adj = tdec.get_adj_list();
#if 0
  for(uint32_t i = 0; i < bags.size(); i++) {
    const auto& b = bags[i];
    cout << "bag id:" << i << endl;
    for(const auto& bb: b) { cout << bb << " "; }
    cout << endl;
  }
  for(uint32_t i = 0; i < adj.size(); i++) {
    const auto& a = adj[i];
    for(const auto& nn: a) cout << i << " " << nn << endl;
  }
#endif
  std::vector<int> dists = tdec.distanceFromCentroid(opt_indep_support_end);
  if (dists.empty()) {
      verb_print(1, "All projected vars in the same bag, ignoring TD");
      return;
  } else {
    int max_dst = 0;
    for(int i=1; i < (int)opt_indep_support_end; i++)
      max_dst = std::max(max_dst, dists[i]);
    if (max_dst == 0) {
      verb_print(1, "All projected vars are the same distance, ignoring TD");
      return;
    }
  }
  sspp::TreeDecomposition dec(bags.size(), opt_indep_support_end);
  for(uint32_t i = 0; i < bags.size();i++) dec.SetBag(i+1, bags[i]);
  for(uint32_t i = 0; i < adj.size(); i++)
    for(const auto& nn: adj[i]) dec.AddEdge(i+1, nn+1);

  // We use 1-indexing, ignore index 0
  auto ord = dec.GetOrd();
  assert(ord.size() == opt_indep_support_end);
  int max_ord = 0;
  int min_ord = std::numeric_limits<int>::max();
  for (uint32_t i = 1; i < opt_indep_support_end; i++) {
    max_ord = std::max(max_ord, ord[i]);
    min_ord = std::min(min_ord, ord[i]);
  }
  max_ord -= min_ord;
  assert(max_ord >= 1);

  // calc td weight
  double rt = 0;
  if (td_width > 0) {
    // Larger the better
    rt = (double)opt_indep_support_end/(double)td_width;
    if (rt*conf.td_exp_mult > 20) td_weight = conf.td_maxweight;
    else td_weight = exp(rt*conf.td_exp_mult)/conf.td_divider;
  } else td_weight = conf.td_maxweight;
  if (conf.do_check_td_vs_ind && (int)indep_support_end < td_width) td_weight = 0.1;
  td_weight = std::min(td_weight, conf.td_maxweight);
  td_weight = std::max(td_weight, conf.td_minweight);
  if (!conf.do_td_weight) td_weight = 1;
  if (print) {
    verb_print(1,
        "TD weight: " << td_weight
        << " opt_end: " << opt_indep_support_end
        << " rt/width(=rt): " << rt
        << " rt*conf.td_exp_mult: " << rt*conf.td_exp_mult
        << " exp(rt*conf.td_exp_mult)/conf.td_divider: "
        << exp(rt*conf.td_exp_mult)/conf.td_divider);
  }

  // Calc td score
  for (uint32_t i = 1; i < opt_indep_support_end; i++) {
    // Normalize
    double val = max_ord - (ord[i]-min_ord);
    val /= (double)max_ord;
    assert(val > -0.01 && val < 1.01);

    assert(i < tdscore.size());
    tdscore[i] = val;
  }

  for(uint32_t i = 1; i < opt_indep_support_end; i++) {
      verb_print(2, "TD var: " << i << " tdscore: " << tdscore[i]);
  }
}

template<typename T>
TWD::TreeDecomposition Counter<T>::td_decompose_component(double mult) {
  auto const& sup_at = decisions.top().super_comp();
  const auto& c = comp_manager->at(sup_at);
  set<uint32_t> active;
  for(uint32_t i = 0; i < c->nVars(); i++) {
    uint32_t var = c->vars_begin()[i];
    active.insert(var);
  }

  TWD::Graph primal(nVars()+1);
  all_lits(i) {
    Lit l(i/2, i%2 == 0);
    for(const auto& l2: watches[l].binaries) {
      if (!l2.red() && l < l2.lit() && val(l) == X_TRI && val(l2.lit()) == X_TRI
          && active.count(l.var()) && active.count(l2.lit().var())) {
        /* debug_print("bin cl: " << l.var() << " " << l2.lit().var()); */
        primal.addEdge(l.var(), l2.lit().var());
      }
    }
  }

  for(const auto& off: long_irred_cls) {
    Clause& cl = *alloc->ptr(off);
    bool sat = false;
    for(Lit& l: cl) {
      if (val(l) == T_TRI) {sat = true; break;}
    }
    if (sat) continue;

    for(uint32_t i = 0; i < cl.sz; i++) {
      const Lit l = cl[i];
      if (!active.count(l.var())) continue;
      if (val(l) == F_TRI) continue;

      for(uint32_t i2 = i+1; i2 < cl.sz; i2++) {
        const Lit l2 = cl[i2];
        if (!active.count(l2.var())) continue;
        if (val(l2) == F_TRI) continue;

        debug_print("bin cl: " <<  l.var() << " " << l2.var());
        primal.addEdge(l.var(), l2.var());
      }
    }
  }

  // run FlowCutter
  verb_print(2, "[td-cmp] FlowCutter is running...");
  TWD::IFlowCutter fc(primal.numNodes(), primal.numEdges(), 0);
  fc.importGraph(primal);

  // Notice that this graph returned is VERY different
  TWD::TreeDecomposition td = fc.constructTD(conf.td_steps, conf.td_lookahead_iters * mult);
  td.centroid(primal.numNodes(), 0);
  verb_print(2, "[td-cmp] FlowCutter FINISHED, TD width: " << td.width());
  return td;
}

template<typename T>
void Counter<T>::td_decompose() {
  double my_time = cpu_time();
  if (indep_support_end <= 3 || nVars() <= 20 || nVars() > conf.td_varlim) {
    verb_print(1, "[td] too many/few vars, not running TD");
    return;
  }

  TWD::Graph primal(nVars()+1);
  all_lits(i) {
    Lit l(i/2, i%2 == 0);
    for(const auto& l2: watches[l].binaries) {
      if (!l2.red() && l < l2.lit()) {
        debug_print("bin cl: " << l.var() << " " << l2.lit().var());
        primal.addEdge(l.var(), l2.lit().var());
      }
    }
  }

  for(const auto& off: long_irred_cls) {
    Clause& cl = *alloc->ptr(off);
    for(uint32_t i = 0; i < cl.sz; i++) {
      for(uint32_t i2 = i+1; i2 < cl.sz; i2++) {
        debug_print("bin cl: " <<  cl[i].var() << " " << cl[i2].var());
        primal.addEdge(cl[i].var(), cl[i2].var());
      }
    }
  }
  for(uint32_t i = opt_indep_support_end; i < nVars()+1; i++) {
    primal.contract(i, conf.td_max_edges*100);
    if (primal.numEdges() > conf.td_max_edges*100 ) break;
  }

  const uint64_t n = (uint64_t)nVars()*(uint64_t)nVars();
  const double density = (double)primal.numEdges()/(double)n;
  const double edge_var_ratio = (double)primal.numEdges()/(double)nVars();
  verb_print(1, "[td] Primal graph  "
    << " nodes: " << primal.numNodes()
    << " edges: " <<  primal.numEdges()
    << " density: " << std::fixed << std::setprecision(3) << density
    << " edge/var: " << std::fixed << std::setprecision(3) << edge_var_ratio);
  if (primal.numEdges() > conf.td_max_edges) {
    verb_print(1, "[td] Too many edges, " << primal.numEdges() << " skipping TD");
    return;
  }
  if (density > 0.3) {
    verb_print(1, "[td] Density is too high, " << density << " skipping TD");
    return;
  }
  if (edge_var_ratio > 30) {
    verb_print(1, "[td] edge/var ratio is too high (" << edge_var_ratio  << "), not running TD");
    return;
  }
  TWD::Graph primal_alt(opt_indep_support_end);
  for(uint32_t i = 0 ; i < opt_indep_support_end; i++) {
    const auto& k = primal.get_adj_list()[i];
    for(const auto& i2: k) {
      if (i2 < (int)opt_indep_support_end)
        primal_alt.addEdge(i, i2);
    }
  }

  // run FlowCutter
  verb_print(2, "[td] FlowCutter is running...");
  TWD::IFlowCutter fc(primal_alt.numNodes(), primal_alt.numEdges(), conf.verb);
  fc.importGraph(primal_alt);

  // Notice that this graph returned is VERY different
  TWD::TreeDecomposition td = fc.constructTD(conf.td_steps, conf.td_iters);

  td.centroid(opt_indep_support_end, conf.verb);
  compute_score(td);
  verb_print(1, "[td] decompose time: " << cpu_time() - my_time);
}

// Self-check count without restart with CMS only
template<typename T>
T Counter<T>::check_count_norestart_cms(const Cube<T>& c) {
  verb_print(1, "Checking cube count with CMS (no verb, no restart)");
  vector<Lit> tmp;
  CMSat::SATSolver test_solver;
  test_solver.new_vars(nVars());
  // Long cls
  for(const auto& off: long_irred_cls) {
    const Clause& cl = *alloc->ptr(off);
    tmp.clear();
    for(const auto& l: cl) tmp.push_back(l);
    test_solver.add_clause(ganak_to_cms_cl(tmp));
  }
  // Bin cls
  all_lits(i) {
    Lit l(i/2, i%2);
    for(const auto& l2: watches[l].binaries) {
      if (l2.irred() && l < l2.lit()) {
        tmp.clear();
        tmp.push_back(l);
        tmp.push_back(l2.lit());
        test_solver.add_clause(ganak_to_cms_cl(tmp));
      }
    }
  }
  // Unit cls
  for(const auto& l: unit_cls) {
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
  T cnt = 0;
  while(true) {
    auto ret = test_solver.solve();
    if (ret == CMSat::l_False) break;
    vector<CMSat::Lit> ban;
    T this_cnt = 1;
    if constexpr (weighted) {
      for(uint32_t i = 0; i < opt_indep_support_end-1; i++) {
        Lit l(i+1, test_solver.get_model()[i] == CMSat::l_True);
        this_cnt *= get_weight(l);
      }
    }
    cnt += this_cnt;
    for(uint32_t i = 0; i < indep_support_end-1; i++) {
      ban.push_back(CMSat::Lit(i, test_solver.get_model()[i] == CMSat::l_True));
    }
    test_solver.add_clause(ban);
  }
  return cnt;
}

// Self-check count without restart
template<typename T>
T Counter<T>::check_count_norestart(const Cube<T>& c) {
  verb_print(1, "Checking count with ourselves (no verb, no restart), CNF: " << c.cnf);
  CounterConfiguration conf2 = conf;
  conf2.do_restart = 0;
  conf2.verb = 0;
  conf2.do_buddy = 0;
  conf2.do_cube_check_count = 0;
  vector<Lit> tmp;
  Counter test_cnt(conf2);
  test_cnt.new_vars(nVars());
  set<uint32_t> tmp_indep;
  for(uint32_t i = 1; i < indep_support_end; i++) tmp_indep.insert(i);
  test_cnt.set_indep_support(tmp_indep);
  if constexpr (weighted) {
    all_lits(i) {
      Lit l(i/2, i%2 == 0);
      test_cnt.set_lit_weight(l, get_weight(l));
    }
  }
  // Long cls
  for(const auto& off: long_irred_cls) {
    const Clause& cl = *alloc->ptr(off);
    tmp.clear();
    for(const auto& l: cl) tmp.push_back(l);
    test_cnt.add_irred_cl(tmp);
  }
  // Bin cls
  all_lits(i) {
    Lit l(i/2, i%2);
    for(const auto& l2: watches[l].binaries) {
      if (l2.irred() && l < l2.lit()) {
        tmp.clear();
        tmp.push_back(l);
        tmp.push_back(l2.lit());
        test_cnt.add_irred_cl(tmp);
      }
    }
  }
  // Unit cls
  for(const auto& l: unit_cls) {
    tmp.clear();
    tmp.push_back(l);
    test_cnt.add_irred_cl(tmp);
  }
  // The cube
  for(const auto&l: c.cnf) {
    tmp.clear();
    tmp.push_back(l.neg());
    test_cnt.add_irred_cl(tmp);
  }
  test_cnt.end_irred_cls();
  vector<Cube<T>> ret;
  return test_cnt.outer_count();
}

template<typename T>
void Counter<T>::disable_smaller_cube_if_overlap(uint32_t i, uint32_t i2, vector<Cube<T>>& cubes) {
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
    verb_print(2, "Two cubes overlap.");
    verb_print(2, "c1: " << c1);
    verb_print(2, "c2: " << c2);
    uint32_t to_disable = cubes[i].cnt > cubes[i2].cnt ? i2 : i;
    cubes[to_disable].enabled = false;
    verb_print(2, "Disabled cube " << cubes[to_disable]);
  }
}

template<typename T>
void Counter<T>::print_and_check_cubes(vector<Cube<T>>& cubes) {
  verb_print(2, "cubes     : ");
  for(const auto&c: cubes) verb_print(2, "-> " << c);
  if (conf.do_cube_check_count) {
    for(const auto& c: cubes) {
      T check_cnt;
      if (conf.do_cube_check_count == 1) check_cnt = check_count_norestart(c);
      else check_cnt = check_count_norestart_cms(c);
      cout << "checking cube [ " << c << " ] ---- check_cnt: " << check_cnt << endl;
      if constexpr (weighted) {
        T diff = check_cnt - c.cnt;
        if (diff/check_cnt > 0.01 || diff/check_cnt < -0.01) assert(false);
      } else assert(check_cnt == c.cnt);
    }
  }
}

template<typename T>
void Counter<T>::disable_cubes_if_overlap(vector<Cube<T>>& cubes) {
  for(uint32_t i = 0; i < cubes.size(); i++) {
    if (!cubes[i].enabled) continue;
    for(uint32_t i2 = i+1; i2 < cubes.size(); i2++) {
      if (!cubes[i2].enabled) continue;
      if (!cubes[i].enabled) continue;
      disable_smaller_cube_if_overlap(i, i2, cubes);
    }
  }
}

template<typename T>
int Counter<T>::cube_try_extend_by_lit(const Lit torem, const Cube<T>& c) {
  verb_print(2, "[cube-ext] Trying to remove " << torem << " from cube " << c);

  // Prop all but torem
  for(const auto& l: c.cnf) {
    if (l == torem) continue;
    if (v_val(l.neg()) == T_TRI) continue;
    if (v_val(l.neg()) == F_TRI) return 0; // don't want to deal with this
    v_enqueue(l.neg());
  }
  bool ret = v_propagate();
  assert(ret);

  if (v_val(torem) == F_TRI) {
    verb_print(2, "[cube-ext] Cube  can have " << torem << " removed, but no count change.");
    return 1;
  }
  if (v_val(torem) != X_TRI) {
    verb_print(1, "[cube-ext] Weeeeirrrddd --- " << torem << " ?????");
    return 0;
  }

  // Check if torem doesn't occur anymore
  //
  // TODO: e.g a=b could be true because of the variables set, which could
  //        remove clauses that torem/~torem is inside. So we could do more...
  for(const auto& l: {torem, torem.neg()}) {
    for(const auto& ws: watches[l].binaries) {
      if (ws.red()) continue;
      if (v_val(ws.lit()) == X_TRI) return 0;
    }
    for(const auto& ws: occ[l.raw()]) {
      Clause& cl = *alloc->ptr(ws.off);
      bool good = false;
      for(const auto& cl_lit: cl) {
        if (v_val(cl_lit) == T_TRI) { ok = true; break;}
      }
      verb_print(2, "[cube-ext] Cube can't have " << torem << " removed");
      if (!good) return 0;
    }
  }
  verb_print(2, "[cube-ext] Cube  can have " << torem << " removed AND count doubled");
  return 100;
}

template<typename T>
bool Counter<T>::clash_cubes(const set<Lit>& c1, const set<Lit>& c2) const {
  for(const auto& l: c1) if (c2.count(l.neg())) return true;
  return false;
}

template<typename T>
void Counter<T>::symm_cubes(vector<Cube<T>>& cubes) {
  vector<Cube<T>> extra_cubes;
  for(const auto& c: cubes) {
    set<Lit> orig_cube(c.cnf.begin(), c.cnf.end());
    if (!c.enabled) continue;

    for(const auto& gen: generators) {
      set<Lit> symm_cube;
      vector<Lit> tmp;
      uint32_t mapped = 0;
      for(auto& l: orig_cube) {
        Lit l2 = l;
        if (gen.count(l) != 0) {
          mapped++;
          l2 = gen.find(l)->second;
          tmp.push_back(l);
        }
        symm_cube.insert(l2);
      }
      if (mapped <= 0) continue; // need at least 1 for clash
      if (symm_cube == orig_cube) continue; // same no clash
      /* if (!clash_cubes(symm_cube, orig_cube)) continue; // must clash */
      verb_print(2, "[rst-symm-map] mapped lits: " << tmp
        << " Old cube:" << orig_cube
        << " New cube:" << symm_cube);
      extra_cubes.push_back(Cube<T>(vector<Lit>(symm_cube.begin(), symm_cube.end()), c.cnt, true));
      stats.num_cubes_symm++;
    }
  }
  cubes.insert(cubes.end(), extra_cubes.begin(), extra_cubes.end());
}

template<typename T>
void Counter<T>::extend_cubes(vector<Cube<T>>& cubes) {
  verb_print(2, "[rst-cube-ext] Extending cubes.");
  assert(occ.empty());
  assert(occ_cls.empty());
  auto my_time = cpu_time();
  const auto before_ext = stats.cube_lit_extend;
  const auto before_rem = stats.cube_lit_rem;

  occ.resize((nVars()+1)*2);
  attach_occ(long_irred_cls, false);

  v_backup();
  vivif_setup();

  for(auto& c: cubes) {
    verb_print(2, "--> Working cube: " << c);
    bool go_again = true;
    while (go_again) {
      go_again = false;
      Cube c2 = c;
      for(const auto& l: c2.cnf) {
        v_new_lev();
        int ret = cube_try_extend_by_lit(l, c2);
        v_backtrack();

        if (ret != 0) {
          if (ret == 100) {
            verb_print(2, COLRED "Cube " << c << " can have " << l << " removed, with cnt change" << COLDEF);
            c.cnt *= 2;
            stats.cube_lit_extend++;
          } else stats.cube_lit_rem++;
          c.cnf.erase(std::find(c.cnf.begin(), c.cnf.end(), l));
          go_again = true;
          break;
        }
      }
    }
    verb_print(2, "--> Final cube:   " << c);
  }
  occ.clear();
  occ_cls.clear();
  v_restore();
  verb_print(2, "[rst-cube-ext] Extended cubes. lit-rem: "
      << setw(4) << stats.cube_lit_rem - before_rem
      << " lit-ext: " << setw(4) << stats.cube_lit_extend - before_ext
      << " T: " << (cpu_time() - my_time));
}

template<typename T>
uint32_t Counter<T>::disable_small_cubes(vector<Cube<T>>& cubes) {
  uint32_t disabled = 0;
  std::sort(cubes.begin(), cubes.end(), [](const Cube<T>& a, const Cube<T>& b) {
      return a.lbd < b.lbd;
  });
  uint32_t enabled_so_far = 0;
  for(auto & c : cubes) {
    if (!c.enabled) continue;
    if (enabled_so_far < conf.max_num_cubes_per_restart
        || c.lbd <= conf.lbd_cutoff_always_keep_cube
        || c.cnf.size() <= conf.lbd_cutoff_always_keep_cube) {
      enabled_so_far++;
      continue;
    } else {
      c.enabled = false;
      disabled++;
    }
  }
  return disabled;
}

template<>
mpz_class Counter<mpz_class>::do_appmc_count() {
  is_approximate = true;
  ApproxMC::AppMC appmc;
  appmc.new_vars(nVars());
  appmc.set_verbosity(std::max<int>(0, conf.verb));
  appmc.set_epsilon(conf.appmc_epsilon);
  appmc.set_delta(conf.delta);
  appmc.set_seed(conf.seed);

  vector<Lit> unit(1);
  for(const auto& l: unit_cls) {
      unit[0] = l;
      appmc.add_clause(ganak_to_cms_cl(unit));
  }

  for(const auto& off: long_irred_cls) {
    const Clause& c = *alloc->ptr(off);
    appmc.add_clause(ganak_to_cms_cl(c));
  }

  vector<Lit> bin(2);
  all_lits(lit_i) {
    Lit l(lit_i/2, lit_i%2);
    for(const auto& l2: watches[l].binaries) {
      if (l2.irred() && l < l2.lit()) {
        bin[0] = l;
        bin[1] = l2.lit();
        if (l2.irred()) appmc.add_clause(ganak_to_cms_cl(bin));
        else appmc.add_red_clause(ganak_to_cms_cl(bin));
      }
    }
  }
  for(const auto& off: long_red_cls) {
    const Clause& c = *alloc->ptr(off);
    appmc.add_red_clause(ganak_to_cms_cl(c));
  }
  vector<uint32_t> indep;
  for(int32_t i = 0; i < (int)indep_support_end-1; i++) {
    assert(i >= 0);
    indep.push_back(i);
  }
  appmc.set_sampl_vars(indep);
  ApproxMC::SolCount appmc_cnt = appmc.count();

  mpz_class num_sols(2);
  mpz_pow_ui(num_sols.get_mpz_t(), num_sols.get_mpz_t(), appmc_cnt.hashCount);
  num_sols *= appmc_cnt.cellSolCount;
  verb_print(1, "[appmc] ApproxMC count: " << num_sols);
  return num_sols;
}

class Timer {
    bool clear = false;
public:
    template<typename Function>
    void set_timeout(Function function, double delay) {
      this->clear = false;
      std::thread t([=, this]() {
          if(this->clear) return;
          std::this_thread::sleep_for(std::chrono::milliseconds((int)(delay*1000.0)));
          if(this->clear) return;
          function();
      });
      t.detach();
    }
    void stop() { this->clear = true; }
};

template<typename T>
T Counter<T>::outer_count() {
  if (!ok) return 0;
  T cnt = 0;
  if constexpr (!weighted) if (conf.appmc_timeout > 0) {
    double time_so_far = cpu_time();
    double set_timeout = std::max<double>(conf.appmc_timeout-time_so_far, 0);
    if (conf.appmc_timeout > 500 && set_timeout < 500) {
      double new_set_timeout = 300;
      verb_print(1, "[appmc] Too little time would be given to ganak: " << set_timeout
          << " adjusting to: " << new_set_timeout);
      set_timeout = new_set_timeout;
    }
    verb_print(1, "[appmc] timeout set to: " << set_timeout);
    Timer t;
    t.set_timeout([=, this]() {
        appmc_timeout_fired = true;
        verb_print(3, "**** ApproxMC timer fired ****");
      }, set_timeout);
  }

  verb_print(1, "Sampling set size: " << indep_support_end-1);
  verb_print(1, "Opt sampling set size: " << opt_indep_support_end-1);
  init_activity_scores();
  if (conf.verb) stats.print_short_formula_info(this);
  auto ret = sat_solver->solve();
  start_time = cpu_time();
  uint32_t next_rst_print = 0;
  bool done = false;
  while(ret == CMSat::l_True) {
    auto cubes = one_restart_count();
    if (cubes.size() == 1 && cubes[0].cnf.empty()) done = true;
    CHECK_PROPAGATED_DO(check_all_propagated_conflicted());
    stats.num_cubes_orig += cubes.size();

    // Extend, tighten, symm, disable cubes
    if (!done) extend_cubes(cubes);
    symm_cubes(cubes);
    print_and_check_cubes(cubes);
    disable_cubes_if_overlap(cubes);
    /* const auto disabled = disable_small_cubes(cubes); */
    /* verb_print(2, "[rst-cube] Disabled " << disabled << " cubes."); */

    // Add cubes to count, Ganak & CMS
    T cubes_cnt_this_rst = 0;
    for(const auto&c: cubes) {
      if (!c.enabled) continue;
      cnt+=c.cnt;
      cubes_cnt_this_rst += c.cnt;
      sat_solver->add_clause(ganak_to_cms_cl(c.cnf));
      stats.num_cubes_final++;
    }
    if (conf.verb >= 2 || stats.num_cache_look_ups > next_rst_print) {
      next_rst_print = stats.num_cache_look_ups + (1ULL*1000LL*1000LL);
      verb_print(1,"[rst-cube] Num restarts: " << stats.num_restarts
          << " orig cubes this rst: " << cubes.size()
          << " total orig cubes: " << stats.num_cubes_orig
          << " total final cubes: " << stats.num_cubes_final
          << " total so far: " << cnt
          << " this rst: " << cubes_cnt_this_rst);
    }

    ret = sat_solver->solve();
    if (ret == CMSat::l_False) {done = true; break;}

    // Add cubes to counter
    for(auto it = cubes.rbegin(); it != cubes.rend(); it++) if (it->enabled) {
      add_irred_cl(it->cnf);
      verb_print(2,  "[rst-cube] added cube CL to GANAK: " << it->cnf << " cnt: " << it->cnt);
    }
    decisions.clear();

    end_irred_cls();
    if (!done && conf.do_vivify && (stats.num_restarts % (conf.vivif_outer_every_n)) == (conf.vivif_outer_every_n-1)) {
      double my_time = cpu_time();
      vivify_all(true, true);
      subsume_all();
      toplevel_full_probe();
      verb_print(2, "[rst-vivif] Outer vivified/subsumed/probed all. T: " << (cpu_time() - my_time));
    }
    if (appmc_timeout_fired) break;
  }

  if (!done && ret == CMSat::l_True) {
    if constexpr (weighted) {
      cout << "ERROR: Not done, so we should be doing appmc, but it's weighted!!!" << endl;
      exit(-1);
    } else cnt += do_appmc_count();
  }
  return cnt;
}

template<typename T>
vector<Cube<T>> Counter<T>::one_restart_count() {
  release_assert(ended_irred_cls && "ERROR *must* call end_irred_cls() before solve()");
  if (indep_support_end == std::numeric_limits<uint32_t>::max()) {
    indep_support_end = nVars()+1;
    opt_indep_support_end = nVars()+1;
  }
  mini_cubes.clear();
  assert(opt_indep_support_end >= indep_support_end);

  if (tdscore.empty() && nVars() > 5 && conf.do_td) {
    tdscore.resize(nVars()+1, 0);
    td_decompose();
  }
  count_loop();
  if (conf.verb) stats.print_short(this, &comp_manager->get_cache());
  return mini_cubes;
}

template<typename T>
void Counter<T>::print_all_levels() {
  cout << COLORG "--- going through all decision levels now, printing comps --" << endl;
  uint32_t dec_lev = 0;
  for(const auto& s: decisions) {
    auto const& sup_at = s.super_comp();
    cout << COLORG "super comp of dec_lev " << dec_lev
      << " is at comp_stack position: " << sup_at
      << " branch var here: " << decisions.at(dec_lev).var
      << " unproc'd comp end: " << decisions.at(dec_lev).get_unproc_comps_end()
      << " remaining comp ofs: " << decisions.at(dec_lev).remaining_comps_ofs()
      << " num unproc'd comps: " << decisions.at(dec_lev).num_unproc_comps()
      << " count: " << decisions.at(dec_lev).total_model_count()
      << " (left: " << decisions.at(dec_lev).left_model_count()
      << " right: " << decisions.at(dec_lev).right_model_count()
      << " active: " << (decisions.at(dec_lev).is_right_branch() ? "right" : "left") << "). -- ";

    const auto& c = comp_manager->at(sup_at);
    cout << COLORG "-> Vars in comp_manager->at(" << sup_at << ")."
      << " num vars: " << c->nVars() << " vars: ";
    for(uint32_t i = 0; i < c->nVars(); i++) cout << c->vars_begin()[i] << " ";
    cout << endl;
    dec_lev++;
  }

  cout << "Full comp stack now." << endl;
  for(uint32_t i = 1; i < comp_manager->get_comp_stack().size(); i++) {
    const auto& c = comp_manager->at(i);
    cout << COLORG "-> Vars in comp_manager->at(" << i << ")."
      << " num vars: " << c->nVars() << " vars: ";
    for(uint32_t i2 = 0; i2 < c->nVars(); i2++) cout << c->vars_begin()[i2] << " ";
    cout << endl;
  }
  cout << COLORG "--- Went through all levels now --" << COLDEF << endl;
}

template<typename T>
void Counter<T>::print_stat_line() {
  if (next_print_stat_cache > stats.num_cache_look_ups) return;
  if (next_print_stat_confl > stats.conflicts) return;
  if (conf.verb) stats.print_short(this, &comp_manager->get_cache());
  next_print_stat_cache = stats.num_cache_look_ups + (20ULL*1000LL*1000LL);
  next_print_stat_confl = stats.conflicts + 150LL*1000LL;
}

template<typename T>
bool Counter<T>::chrono_work() {
  debug_print("--- CHRONO CHECK ----");
  VERBOSE_DEBUG_DO(print_trail());
  auto data = find_conflict_level(confl_lit);
  if (data.bOnlyOneLitFromHighest) {
    debug_print(COLYEL2  << __func__ << " -- going back to " << data.nHighestLevel-1
        << " curlev: " << dec_level());
    go_back_to(data.nHighestLevel-1);
    VERBOSE_DEBUG_DO(print_trail());
    debug_print(COLYEL2  << __func__ << " went back -- now Declev: " << dec_level());
    return true;
  }
  return false;
}

template<typename T>
void Counter<T>::count_loop() {
  assert(mini_cubes.empty());
  RetState state = RESOLVED;

  while (true) {
    debug_print("var top of decision stack: " << decisions.top().var);
    // NOTE: find_next_remain_comp_of finds disjoint comps
    // we then solve them all with the decide_lit & calling findNext.. again
    while (comp_manager->find_next_remain_comp_of(decisions.top())) {
      if (!decide_lit()) {
        decisions.top().next_unproc_comp();
        continue;
      }

      if (!isindependent) {
        // The only decision we could make would be non-indep for this component.
        debug_print("before SAT mode. cnt dec: " << decisions.top().total_model_count()
            << " left: " << decisions.top().left_model_count()
            << " right: " << decisions.top().right_model_count());

        bool sat = use_sat_solver(state);
        debug_print("after SAT mode. cnt dec: " << decisions.top().total_model_count()
            << " left: " << decisions.top().left_model_count()
            << " right: " << decisions.top().right_model_count());
        if (sat) {
          state = BACKTRACK;
        } else {
          goto start11;
        }
        debug_print("after SAT mode. cnt of this comp: " << decisions.top().total_model_count()
          << " unproc comps end: " << decisions.top().get_unproc_comps_end()
          << " remaining comps: " << decisions.top().remaining_comps_ofs()
          << " has unproc: " << decisions.top().has_unproc_comps());
        assert(isindependent);

        // Now backtrack
        break;
      }

      if (conf.do_buddy && should_do_buddy_count()) {
        if (do_buddy_count()) {
          state = BACKTRACK;
          break;
        } else {
          goto start11;
        }
      }

      while (!propagate()) {
        start1:
        if (chrono_work()) continue; // will DEFINITELY conflict if TRUE
        state = resolve_conflict();
        start11:
        if (state == GO_AGAIN) goto start1;
        if (state == BACKTRACK) break;
      }
      if (state == BACKTRACK) break;
      if (state == RESOLVED && restart_if_needed()) goto end;

      // we are in RESOLVED or PROCESS_COMPONENT state, continue.
      // first time bug: old/out-ganak-6870225.pbs101-9/mc2023_track1_174.cnf.gz.out_ganak
      //       which is cc5fbdbad21351745d004b4b39e3f73244a42ecd
      //       against, say: 71ba4c4eaf19eb74ec9fb3a2e6701ebebe71b9e5
      //        --> seems like SAT solver, chrono work...
      if (state != PROCESS_COMPONENT && state != RESOLVED) cout << "ERROR: state: " << state << endl;
      assert(state == PROCESS_COMPONENT || state == RESOLVED);
    }
    // we are here because there is no next component, or we had to backtrack

    print_stat_line();
    state = backtrack();
    if (state == EXIT) goto end;

    while (!propagate()) {
      start2:
      if (chrono_work()) continue;
      state = resolve_conflict();
      if (state == GO_AGAIN) goto start2;
      if (state == BACKTRACK) {
        state = backtrack();
        if (state == EXIT) goto end;
      }
    }
    assert(state != GO_AGAIN);

    if (conf.do_vivify) {
      vivify_all();
      bool ret = propagate();
      assert(ret);
    }
  }

end:
  if (state == EXIT) {
    Cube<T> c(vector<Lit>(), decisions.top().total_model_count());
    debug_print("Exiting due to EXIT state, the cube count: " << c.cnt);
    mini_cubes.push_back(c);
  } else {/*restart*/}

  if constexpr (weighted) {
    T this_restart_multiplier = 1;
    for(uint32_t i = 1; i < opt_indep_support_end; i++)
      if (!is_unknown(i)) {
        Lit l(i, val(i) == T_TRI);
        this_restart_multiplier *= get_weight(l);
        debug_print("[cube-final] lit: " <<  l << " mul: " << get_weight(l));
      }
    debug_print("[cube-final] This restart multiplier: " << this_restart_multiplier);
    for (auto& c: mini_cubes) {
      c.cnt *= this_restart_multiplier;
      if (c.enabled) debug_print("[cube-final] cube: " << c);
    }
  }

  // We have to propagate at the end due to non-chrono BT
  bool ret = propagate();
  assert(ret && "never UNSAT");
}

template<typename T>
void Counter<T>::recomp_td_weight() {
  if (conf.td_lookahead != -1 && dec_level() < conf.td_lookahead+5) {
    auto td = td_decompose_component(3);
    compute_score(td, false);
  }
}

template<typename T>
bool Counter<T>::standard_polarity(const uint32_t v) const {
  if (watches[Lit(v, true)].activity == watches[Lit(v, false)].activity)
    return var(Lit(v, true)).last_polarity;
  return watches[Lit(v, true)].activity > watches[Lit(v, false)].activity;
}

template<typename T>
bool Counter<T>::get_polarity(const uint32_t v) const {
  bool polarity;
  switch (conf.polar_type) {
    case 0: polarity = standard_polarity(v); break;
    case 1: polarity = var(v).last_polarity; break;
    case 2: polarity = false; break;
    case 3: polarity = true; break;
    default: assert(false);
  }
  return polarity;
}

template<typename T>
bool Counter<T>::decide_lit() {
  recomp_td_weight();
  VERBOSE_DEBUG_DO(print_all_levels());
  debug_print("new decision level is about to be created, lev now: " << dec_level() << " branch: " << decisions.top().is_right_branch());
  decisions.push_back(
    StackLevel<T>(decisions.top().curr_remain_comp(),
               comp_manager->comp_stack_size()));

  // The decision literal is now ready. Deal with it.
  uint32_t v = 0;
  isindependent = true;
  switch (conf.decide) {
    case 0: v = find_best_branch(false); break;
    case 1: v = find_best_branch(true); break;
    default: assert(false);
  }
  if (v == 0) {
    decisions.pop_back();
    isindependent = false;
    return true;
  }
  assert(val(v) == X_TRI);

  decisions.top().var = v;

  Lit lit = Lit(v, get_polarity(v));
  /* cout << "decided on: " << std::setw(4) << lit.var() << " sign:" << lit.sign() <<  endl; */
  debug_print(COLYEL "decide_lit() is deciding: " << lit << " dec level: "
      << dec_level());
  set_lit(lit, dec_level());
  stats.decisions++;
  vsads_readjust();
  assert( decisions.top().remaining_comps_ofs() <= comp_manager->comp_stack_size());
  return true;
}

template<typename T>
double Counter<T>::var_act(const uint32_t v) const {
  auto w = (watches[Lit(v, false)].activity + watches[Lit(v, true)].activity);
  return w;
}

// The higher, the better. It is never below 0.
template<typename T>
double Counter<T>::score_of(const uint32_t v, bool ignore_td) const {
  bool print = false;
  /* if (stats.decisions % 40000 == 0) print = 1; */
  /* print = true; */
  /* print = false; */
  double act_score = 0;
  double td_score = 0;
  double freq_score = 0;
  vector<uint32_t> occ_cnt; // number of occurrences of a variable in the component

  // TODO Yash idea: let's cut this into activities and incidence
  if (!tdscore.empty() && !ignore_td) td_score = td_weight*tdscore[v];
  act_score = var_act(v)/3.0;
  freq_score = (double)comp_manager->freq_score_of(v)/25.0;
  double score = act_score+td_score+freq_score;
  if (print) cout << "v: " << std::setw(4) << v
    << std::setw(3) << " conflK: " << stats.conflicts/1000
    << std::setw(5) << " decK: " << stats.decisions/1000
    << std::setw(6) << " act_score: " << act_score/score
    << std::setw(6) << " freq_score: " << freq_score/score
    << std::setw(6) << " td_score: " << td_score/score
    << std::setw(6) << " total: " << score
    << std::setw(6) << endl;

  return score;
}

template<typename T>
double Counter<T>::td_lookahead_score(const uint32_t v, const uint32_t base_comp_tw) {
  double score = 0;
  auto my_time = cpu_time();

  int32_t w[2];
  int tdiff[2];
  for(bool b: {true, false}) {
    set_lit(Lit(v, b), dec_level());
    int tsz = trail.size();
    bool ret = propagate();
    if (!ret) {
      score = 1e5;
      reactivate_comps_and_backtrack_trail();
      return score;
    }
    tdiff[b] = trail.size()-tsz;
    if (tdiff[b] < 3) w[b] = base_comp_tw;
    else w[b] = td_decompose_component().width();
    reactivate_comps_and_backtrack_trail();
  }
  verb_print(1, "var: " << setw(4) << v << " w[0]: " << setw(4) << w[0]
    << " w[1]: " << setw(4) << w[1]
    << " trail diff: " << setw(3) << tdiff[0]
    << " trail diff: " << setw(3) << tdiff[1]
    << " T: " << (cpu_time()-my_time));
  /* return tdiff[0]*tdiff[1]; */
  return -1*std::max<int32_t>({w[0],w[1]})*w[0]*w[1];
}

template<typename T>
uint32_t Counter<T>::find_best_branch(bool ignore_td) {
  bool only_optional_indep = true;
  uint32_t best_var = 0;
  double best_var_score = -1e8;
  uint64_t* at;
  VERBOSE_DEBUG_DO(cout << "decision level: " << dec_level() << " var options: ");
  if constexpr (weighted) {
    if (vars_act_dec.size()  < (dec_level()+1) * (nVars()+1)) {
      uint64_t todo = (dec_level()+1)*(nVars()+1) - vars_act_dec.size();
      vars_act_dec.insert(vars_act_dec.end(), todo, 0);
    }
    at = vars_act_dec.data()+(nVars()+1)*dec_level();
    vars_act_dec_num++;
    at[0] = vars_act_dec_num;
    VERBOSE_DEBUG_DO(cout << "(at[0] = " << at[0] << ") ");
  }

  int32_t tw = 0;
  if (dec_level() < conf.td_lookahead && !conf.td_look_only_weight)
    tw = td_decompose_component().width();

  all_vars_in_comp(comp_manager->get_super_comp(decisions.top()), it) {
    const uint32_t v = *it;
    if (val(v) != X_TRI) continue;
    VERBOSE_DEBUG_DO(cout << v << " ");

    if (v < opt_indep_support_end) {
      if (v < indep_support_end) only_optional_indep = false;
      if constexpr (weighted) at[v] = vars_act_dec_num;
      double score;
      if (!conf.td_look_only_weight && dec_level() < conf.td_lookahead &&
          tw > conf.td_lookahead_tw_cutoff)
        score = td_lookahead_score(v, tw);
      else score = score_of(v, ignore_td) ;
      if (best_var == 0 || score > best_var_score) {
        best_var = v;
        best_var_score = score;
      }
    }
  }
  VERBOSE_DEBUG_DO(cout << endl);

  if (best_var != 0 && only_optional_indep) return 0;
  if (dec_level() < conf.td_lookahead && tw > conf.td_lookahead_tw_cutoff)
    verb_print(1, "best var: " << best_var << " score: " << best_var_score);
  return best_var;
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

template<typename T>
bool Counter<T>::restart_if_needed() {
  if (!appmc_timeout_fired && conf.max_num_rst > 0 && (int32_t)stats.num_restarts > conf.max_num_rst)
    return false;
  if (!appmc_timeout_fired && (!conf.do_restart || td_width < 60)) return false;

  bool restart = false;
  if (appmc_timeout_fired) {
    verb_print(1, "[rst] AppMC timeout fired. Restarting.");
    restart = true;
  }

  // Conflicts, luby
  if (conf.restart_type == 7) {
    verb_print(3, "[rst] Will restart at confl: " << luby(2, stats.num_restarts) * conf.first_restart
      << " now confl: " << stats.conflicts);
    if ((stats.conflicts-stats.last_restart_num_conflicts) >
        (luby(2, stats.num_restarts) * conf.first_restart)) {
      verb_print(2, "[rst] restarting. Next restart confl: "
          << (stats.conflicts + luby(2, stats.num_restarts+1) * conf.first_restart));
      restart = true;
    }
  }

  // Decisions, luby
  if (conf.restart_type == 8) {
      verb_print(3, "[rst] Will restart at decK: "
        << (stats.decisions + luby(2, stats.num_restarts) * conf.first_restart * 20)/1000.0
        << " now decK: " << stats.decisions/1000.0);

    if ((stats.decisions-stats.last_restart_num_decisions) >
      (luby(2, stats.num_restarts) * conf.first_restart * 20)) {
      verb_print(2, "[rst] restarting. Next restart decK: "
        << (stats.decisions + luby(2, stats.num_restarts+1) * conf.first_restart * 20)/1000);
      restart = true;
    }
  }

  // Comps, luby
  if (conf.restart_type == 9) {
    verb_print(3, "[rst] Will restart at comps: "
        << (stats.num_cached_comps + luby(2, stats.num_restarts) * conf.first_restart * 1000)
        << " now comps: " << stats.num_cached_comps);
    if ((stats.num_cached_comps) > (luby(2, stats.num_restarts) * conf.first_restart*1000)) {
      verb_print(2, "[rst] restarting. Next restart comps: "
        << (stats.num_cached_comps + luby(2, stats.num_restarts+1) * conf.first_restart * 1000));
      restart = true;
    }
  }

  if (!restart) return false;
  verb_print(2, "************* Restarting.  **************");
  verb_print(2, "Num decisions since last restart: "
    << stats.decisions-stats.last_restart_num_decisions << endl
    << "c o Num conflicts since last restart: "
    << stats.conflicts-stats.last_restart_num_conflicts << endl
    << "c o Num comps since last restart: "
    << stats.num_cache_look_ups-stats.last_restart_num_cache_look_ups);

  // Reset stats
  stats.last_restart_num_conflicts = stats.conflicts;
  stats.last_restart_num_decisions = stats.decisions;
  stats.last_restart_num_cache_look_ups = stats.num_cache_look_ups;

  assert(mini_cubes.empty());
  T tot_cnt = 0;
  while (dec_level() > 0) {
    verb_print(2, COLBLBACK <<  COLCYN "--> Mini cube gen. "
      << " lev: " << dec_level()
      << " left cnt: " << decisions.top().left_model_count()
      << " right cnt: " << decisions.top().right_model_count()
      << COLDEF);
    for(auto i: {0, 1}) {
      if (decisions.top().get_model_side(i) == 0) continue;
      verb_print(2, "->> branch: " << i << " doing compute_cube...");

      Cube<T> cube;
      if (compute_cube(cube, i)) {
        mini_cubes.push_back(cube);
        tot_cnt += cube.cnt;
        verb_print(2, "[mini-cube] rst: " << stats.num_restarts << " mini cube: " << cube);
      }
      else comp_manager->removeAllCachePollutionsOfIfExists(decisions.top());
    }
    reactivate_comps_and_backtrack_trail(false);
    bool ret = propagate(true);
    assert(ret);
    decisions.pop_back();
    VERY_SLOW_DEBUG_DO(if (!check_watchlists()) {print_trail(false, false);assert(false);});
  }
  verb_print(2, "[mini-cube] rst: " << stats.num_restarts << " tot cnt before overlap: " << tot_cnt);

  // Because of non-chrono backtrack, we need to propagate here:
  // zero decision level stuff now gets propagated at 0-level
  bool ret = propagate();
  assert(ret && "never UNSAT");
  CHECK_PROPAGATED_DO(check_all_propagated_conflicted());
  stats.num_restarts++;

  // Readjust
  if (conf.do_readjust_for_restart) {
    conf.decide = stats.num_restarts%2;
    /* conf.polar_type = (stats.num_restarts % 5 == 3) ? (stats.num_restarts%4) : 0; */
  }
  verb_print(2, "[rst] new config. decide: " << conf.decide
    << " polar_type: " << conf.polar_type);
  return true;
}

// returns cube in `c`. Uses LEFT/RIGHT branch
// if UNSAT that SAT solver figured out, returns false
//    in this case, the cache elements much be deleted (they are erroneous)
template<typename T>
bool Counter<T>::compute_cube(Cube<T>& c, const int side) {
  assert(c.cnt == 0);
  assert(c.cnf.empty());
  debug_print(COLWHT "-- " << __func__ << " BEGIN");

  c.cnt = decisions.top().get_model_side(side);
  debug_print("Own cnt: " << c.cnt);
  for(int32_t i = 0; i < dec_level(); i++) {
    const auto& dec = decisions[i];
    const auto& mul = dec.get_branch_sols(); // ACTIVE branch (i.e. currently counted one)
    if (mul == 0) continue;
    c.cnt*=mul;
  }
  debug_print("Mult cnt: " << c.cnt);
  if (c.cnt == 0) return false;

  const bool opposite_branch = side != decisions.top().is_right_branch();

  // Add decisions
  debug_print(COLWHT "Indep decisions in the c.cnf: ");
  for(const auto& l: trail) {
    if (!var(l).ante.isNull()) continue;
    if (l.var() >= opt_indep_support_end) continue;
    if (var(l).decision_level == dec_level() && opposite_branch) {
      assert(l == top_dec_lit());
      c.cnf.push_back(l);
    } else {
      c.cnf.push_back(l.neg());
    }
    debug_print_noendl(l << " ");
    c.lbd = calc_lbd(c.cnf);
  }
  debug_print_noendl(COLDEF << endl);

  // Get a solution
  vector<CMSat::Lit> ass; ass.reserve(c.cnf.size());
  for(const auto&l: c.cnf) ass.push_back(CMSat::Lit(l.var()-1, l.sign()));
  auto solution = sat_solver->solve(&ass);
  debug_print("cube solution: " << solution);
  if (solution == CMSat::l_False) return false;

  // Add values for all components not yet counted
  for(int32_t i = 0; i <= dec_level(); i++) {
    if (i == dec_level() && opposite_branch) {
      // This has been fully counted, ALL components.
      continue;
    }
    const StackLevel<T>& dec = decisions[i];
    const auto off_start = dec.remaining_comps_ofs();
    const auto off_end = dec.get_unproc_comps_end();
    debug_print("lev: " << i << " off_start: " << off_start << " off_end: " << off_end);
    // add all but the last component (it's the one being counted lower down)
    int off_by_one = 1;
    if (i == dec_level()) off_by_one = 0;
    for(uint32_t i2 = off_start; i2 < off_end-off_by_one; i2++) {
      const auto& comp = comp_manager->at(i2);
      all_vars_in_comp(*comp, v) {
        Lit l = Lit(*v, sat_solver->get_model()[*v-1] == CMSat::l_False);
        debug_print("Lit from comp: " << l);
        if (l.var() >= indep_support_end) continue;
        c.cnf.push_back(l);
      }
    }
  }

#ifdef VERBOSE_DEBUG
  // Show decision stack's comps
  for(int32_t i = 0; i <= dec_level(); i++) {
    const auto& dst = decisions.at(i);
    cout << COLWHT "decisions.at(" << i << "):"
      << " decision var: " << dst.var
      << " num unproc comps: " << dst.num_unproc_comps()
      << " unproc comps end: " << dst.get_unproc_comps_end()
      << " remain comps offs: " << dst.remaining_comps_ofs()
      << " total count here: " << dst.total_model_count()
      << " left count here: " << dst.left_model_count()
      << " right count here: " << dst.right_model_count()
      << " branch: " << dst.is_right_branch() << endl;
    const auto off_start = dst.remaining_comps_ofs();
    const auto off_end = dst.get_unproc_comps_end();
    for(uint32_t i2 = off_start; i2 < off_end; i2++) {
      assert(i2 < comp_manager->comp_stack_size());
      const auto& comp = comp_manager->at(i2);
      cout << COLWHT "-> comp at: " << std::setw(3) << i2 << " ID: " << comp->id() << " -- vars : ";
      all_vars_in_comp(*comp, v) cout << *v << " ";
      cout << COLDEF << endl;
    }
  }

  cout << COLORG "cube so far. Size: " << c.cnf.size() << " cube: ";
  for(const auto& l: c.cnf) cout << l << " ";
  cout << endl;
  cout << COLORG "cube's SOLE count: " << decisions.top().get_model_side(side) << endl;
  cout << COLORG "cube's RECORDED count: " << c.cnt << COLDEF << endl;
#endif
  return true;
}

// Checks one-by-one using a SAT solver
template<typename T>
T Counter<T>::check_count(const bool also_incl_curr_and_later_dec) {
    //let's get vars active
    set<uint32_t> active;

    const auto& s = decisions.top();
    auto const& sup_at = s.super_comp();
    const auto& c = comp_manager->at(sup_at);
#ifdef VERBOSE_DEBUG
    cout << "-> Checking count. also_incl_curr_and_later_dec: " COLRED << std::boolalpha <<  also_incl_curr_and_later_dec << COLDEF
      << " dec lev: " << dec_level() << " [ stats: decisions so far: " << stats.decisions << " confl so far: " << stats.conflicts << " ]" << endl;
    cout << "-> Vars in comp_manager->at(" << sup_at << ")."
      << " num vars: " << c->nVars() << " vars: ";
    for(uint32_t i = 0; i < c->nVars(); i++) cout << c->vars_begin()[i] << " ";
    cout << endl;
#endif

    T dec_w = 1;
    for(uint32_t i = 0; i < c->nVars(); i++) {
      uint32_t v = c->vars_begin()[i];
      if (v < opt_indep_support_end) {
        active.insert(v);
        if constexpr (weighted) if (val(v) != X_TRI && var(v).decision_level == dec_level()) {
            dec_w *= get_weight(Lit(v, val(v) == T_TRI));
            if (get_weight(Lit(v, val(v) == T_TRI)) != 1)
              debug_print(COLYEL "mult var: " << setw(4) << v << " val: " << setw(3) << val(v)
                << " weight: " << std::setw(9) << get_weight(Lit(v, val(v) == T_TRI)) << COLDEF
                << " dec_lev: " << var(v).decision_level);
        }
      }
    }

#ifdef VERBOSE_DEBUG
    cout << "active for count chk: "; for(const auto&a: active) cout << a << " "; cout << endl;
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
    debug_print("dec lev: " << dec_level());
    VERBOSE_DEBUG_DO(if (!trail.empty()) cout << "top dec lit: " << top_dec_lit() << endl;);
    CMSat::SATSolver s2;
    CMSat::copy_solver_to_solver(sat_solver, &s2);
    int32_t last_dec_lev = -1;
    for(const auto& t: trail) {
      last_dec_lev = std::max(last_dec_lev, var(t.var()).decision_level);
    }
    for(const auto& t: trail) {
      if (!also_incl_curr_and_later_dec) {
        if (var(t).decision_level >= dec_level()) continue;
      }
      // don't include propagations or lev0 stuff
      if (!var(t).ante.isNull() || var(t).decision_level == 0) continue;
      cl.clear();
      cl.push_back(CMSat::Lit(t.var()-1, !t.sign()));
      s2.add_clause(cl);
    }
    T cnt = 0;
    bool solution_exist = false;
    while(true) {
      auto ret = s2.solve();
      if (ret == CMSat::l_True) {
        solution_exist = true;
        if constexpr (!weighted) cnt++;
        else {
          T cube_cnt = 1;
          for(uint32_t i = 0; i < s2.nVars(); i++) {
            if (active.count(i+1)
                && (val(i+1) == X_TRI || var(i+1).decision_level >= dec_level())
                ) {
              cube_cnt *= get_weight(Lit(i+1, s2.get_model()[i] == CMSat::l_True));
            }
          }
          VERBOSE_DEBUG_DO(cout << cube_cnt << " + ";);
          cnt += cube_cnt;
        }

        // Ban solution
        cl.clear();
        for(uint32_t i = 0; i < s2.nVars(); i++) {
          if (active.count(i+1)) {
            CMSat::Lit l = CMSat::Lit(i, s2.get_model()[i] == CMSat::l_False);
            cl.push_back(~l);
          }
        }
        /* cout << "banning sol: " << cl << endl; */
        s2.add_clause(cl);
      } else if (ret == CMSat::l_False) break;
      else assert(false);
    }
    VERBOSE_DEBUG_DO(cout << endl);
    T after_mul = 0;
    if (!decisions.top().is_right_branch()) {
      after_mul += decisions.top().left_model_count()*dec_w;
      after_mul += decisions.top().right_model_count();
    } else {
      after_mul += decisions.top().left_model_count();
      after_mul += decisions.top().right_model_count()*dec_w;
    }
    debug_print("correct                            : " << std::setprecision(10) << cnt);
    debug_print("after_mul:                         : " << after_mul);
    debug_print("dec_w                              : " << dec_w);
    debug_print("active                             : " << (decisions.top().is_right_branch() ? "right" : "left"));
    debug_print("ds.top().left_model_count()    : " << decisions.top().left_model_count());
    debug_print("ds.top().right_model_count()   : " << decisions.top().right_model_count());

    // It can be that a subcomponent above is UNSAT, in that case, it'd be UNSAT
    // and the count cannot be checked
    if (solution_exist) {
      if constexpr (!weighted) assert(decisions.top().total_model_count() == cnt);
      else {
        bool okay = true;
        T diff = after_mul - cnt;
        if (diff != 0) {
          debug_print("OOps, diff              : " << diff << " diff ratio: " << diff/cnt);
          print_all_levels();
          okay = false;
        }
        assert(okay);
      }
    }
    cout << std::setprecision(3);
    return cnt;
}

template<typename T>
RetState Counter<T>::backtrack() {
  debug_print("in " << __FUNCTION__ << " now. Dec lev: " << dec_level());
  assert(decisions.top().remaining_comps_ofs() <= comp_manager->comp_stack_size());
  do {
#ifdef VERBOSE_DEBUG
    if (dec_level() > 0) {
      debug_print("[indep] top count here: " << decisions.top().total_model_count()
        << " left: " << decisions.top().left_model_count()
        << " right: " << decisions.top().right_model_count()
        << " is right: " << decisions.top().is_right_branch()
        << " dec lit: " << top_dec_lit()
        << " dec lev: " << dec_level());
    }
#endif
    if (decisions.top().branch_found_unsat()) {
      comp_manager->removeAllCachePollutionsOf(decisions.top());
    } else if (decisions.top().another_comp_possible()) {
      debug_print("[indep] Processing another comp at dec lev "
          << dec_level()
          << " instead of backtracking." << " Num unprocessed comps: "
          << decisions.top().num_unproc_comps()
          << " so far the count: " << decisions.top().total_model_count());
      return PROCESS_COMPONENT;
    }

    // We have NOT explored the other side and it hasn't been re-written to be
    // propagation.
    if (!decisions.top().is_right_branch() && var(top_dec_lit()).ante.isNull()) {
      debug_print("[indep] We have NOT explored the right branch (isSecondBranch==false). Let's do it!"
          << " -- dec lev: " << dec_level());
      const Lit lit = top_dec_lit();
      assert(dec_level() > 0);
      CHECK_COUNT_DO(check_count(true));
      SLOW_DEBUG_DO(assert(decisions.top().right_model_count() == 0));
      // could be the flipped that's FALSEified so that would
      // mean the watchlist is not "sane". We need to propagate the flipped var and
      // then it'll be fine
      /* CHECK_PROPAGATED_DO(check_all_propagated_conflicted()); */

      // NOTE: replacing a decision literal x with y when y->x binary clause exists does
      // not work, because we'll count (x, y) = 01 (left hand branch), and
      // 10 (right hand branch, setting y = 0, forcing x = 1), but not 11.
      reactivate_comps_and_backtrack_trail(false);
      decisions.top().change_to_right_branch();
      bool ret = propagate(true);
      assert(ret);
      debug_print("[indep] Flipping lit to: " << lit.neg() << " val is: " << val_to_str(val(lit)));
      if (val(lit.neg()) == X_TRI) {
        set_lit(lit.neg(), dec_level());
        VERBOSE_DEBUG_DO(print_trail());
        debug_print(COLORGBG "[indep] Backtrack finished -- we flipped the branch. "
            "count left: " << decisions.top().left_model_count()
            << " count right: " << decisions.top().right_model_count());
        return RESOLVED;
      } else {
        assert(val(lit.neg()) == F_TRI && "Cannot be TRUE because that would mean that the branch we just explored was UNSAT and we should have detected that");
        decisions.top().mark_branch_unsat();
        continue;
      }
    }
    debug_print(COLORGBG "[indep] We have explored BOTH branches, actually BACKTRACKING."
        << " -- dec lev: " << dec_level());
    // Backtrack from end, i.e. finished.
    if (dec_level() == 0) {
      debug_print(COLORGBG "[indep] Backtracking from lev 0, i.e. ending");
      CHECK_COUNT_DO(check_count());
      break;
    }

    CHECK_COUNT_DO(check_count());
    reactivate_comps_and_backtrack_trail(false);
    assert(dec_level() >= 1);
    if (conf.do_use_cache) {
#ifdef VERBOSE_DEBUG
      cout << "comp vars: ";
      all_vars_in_comp(comp_manager->get_super_comp(decisions.top()), it) {
        cout << *it << " ";
      }
      cout << endl;
#endif
      if constexpr (weighted) {
        T cnt = decisions.top().total_model_count();
        all_vars_in_comp(comp_manager->get_super_comp(decisions.top()), it) {
          if (val(*it) != X_TRI && var(*it).decision_level < dec_level()) {
            Lit l(*it, val(*it) == T_TRI);
            if (get_weight(l) != 1) {
              debug_print(COLYEL2 << "MULT STORE var: " << std::setw(3) << *it
                << " val: " << val_to_str(val(*it))
                << " dec lev: " << var(*it).decision_level);
              cnt *= get_weight(Lit(*it, val(*it) == T_TRI));
            }
          }
        }
        comp_manager->save_count(decisions.top().super_comp(), cnt);
      } else {
        comp_manager->save_count(decisions.top().super_comp(), decisions.top().total_model_count());
      }
    }

#ifdef VERBOSE_DEBUG
    const auto parent_count_before = (decisions.end() - 2)->total_model_count();
    const auto parent_count_before_left = (decisions.end() - 2)->left_model_count();
    const auto parent_count_before_right = (decisions.end() - 2)->right_model_count();
#endif
    (decisions.end() - 2)->include_solution(decisions.top().total_model_count());
    decisions.pop_back();

    // var == 0 means it's coming from a fake decision due to normal SAT solving
    assert(decisions.top().var == 0 || decisions.top().var < opt_indep_support_end);
    auto& dst = decisions.top();
    debug_print("[indep] -> Backtracked to level " << dec_level()
        // NOTE: -1 here because we have JUST processed the child
        //     ->> (see below next_unproc_comp() call)
        << " num unprocessed comps here: " << dst.num_unproc_comps()-1
        << " current count here: " << dst.total_model_count()
        << " branch: " << dst.is_right_branch()
        << " before including child it was: " <<  parent_count_before
        << " (left: " << parent_count_before_left
        << " right: " << parent_count_before_right
        << ")");

    // step to the next comp not yet processed
    dst.next_unproc_comp();

    assert(dst.remaining_comps_ofs() < comp_manager->comp_stack_size() + 1);
  } while (true);
  return EXIT;
}

template<typename T>
void Counter<T>::print_dec_info() const {
  cout << "dec lits: " << endl;
  for(uint32_t i = 1; i < decisions.size(); i ++) {
    uint32_t dvar = decisions[i].var;
    Lit l = Lit(dvar, val(dvar) == T_TRI);
    cout << "dec lev: " << std::setw(3) << i <<
      " lit: " << std::setw(6)
      << l
      << " is right: "
      << (int)decisions[i].is_right_branch()
      << " ante: " << std::setw(10) << var(l).ante
      << " lev: " << var(l).decision_level
      << endl;
  }
}

template<typename T>
void Counter<T>::print_conflict_info() const
{
  print_dec_info();
  cout << "UIP cl lits: " << endl;
  print_cl(uip_clause);
  cout << "uip_clause[0]: " << uip_clause[0] << endl;
}

template<typename T>
void Counter<T>::print_comp_stack_info() const {
    cout << "decisions.top().remaining_comps_ofs(): "
      << decisions.top().remaining_comps_ofs() << endl;
    cout << "comp_manager->comp_stack_size(): " <<
      comp_manager->comp_stack_size() << endl;
}

struct UIPFixer {
  UIPFixer(vector<VarData>& _vars) : vars(_vars) { }
  bool operator()(const Lit& a, const Lit& b) const {
    auto a_dec = vars[a.var()].decision_level;
    auto b_dec = vars[b.var()].decision_level;
    if (a_dec != b_dec) return a_dec > b_dec;
    auto a_ante = vars[a.var()].ante;
    /* auto b_ante = vars[b.var()].ante; */
    if (!a_ante.isAnt()) return true;
    return false;
  }
  vector<VarData>& vars;
};

template<typename T>
int32_t Counter<T>::find_backtrack_level_of_learnt() {
  assert(!uip_clause.empty());
  uint32_t max_i = 0;
  for (uint32_t i = 0; i < uip_clause.size(); i++) {
    if (var(uip_clause[i]).decision_level > var(uip_clause[max_i]).decision_level)
      max_i = i;
  }
  std::swap(uip_clause[max_i], uip_clause[0]);
  return var(uip_clause[0]).decision_level;
}

template<typename T>
int32_t Counter<T>::find_lev_to_set(const int32_t backj) {
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

template<typename T>
void Counter<T>::print_trail(bool check_entail, bool check_anything) const {
  cout << "Current trail :" << endl;
  for(uint32_t i = 0; i < trail.size(); i++) {
    const auto l = trail[i];
    cout << "lit " << std:: left << std::setw(6) << l
      << " lev: " << std::setw(4) << var(l).decision_level
      << " ante: " << std::setw(5) << std::left << var(l).ante
    << " val: " << std::setw(8) << lit_val_str(l)
    << " trail pos: " << std::setw(4) << i
    << " sublevel: "  << std::setw(3) << var(l).sublevel << endl;
  }
  cout << "qhead: " << qhead << endl;
  if (check_anything) check_trail(check_entail);
}

template<typename T>
void Counter<T>::go_back_to(int32_t backj) {
  debug_print("going back to lev: " << backj << " dec level now: " << dec_level());
  while(dec_level() > backj) {
    debug_print("at dec lit: " << top_dec_lit() << " lev: " << dec_level() << " cnt:" <<  decisions.top().total_model_count());
    VERBOSE_DEBUG_DO(print_comp_stack_info());
    decisions.top().mark_branch_unsat();
    decisions.top().zero_out_all_sol(); //not sure it's needed
    if (!sat_mode()) {
      comp_manager->removeAllCachePollutionsOf(decisions.top());
    }
    reactivate_comps_and_backtrack_trail(false);
    decisions.pop_back();
    decisions.top().zero_out_branch_sol();
    if (!sat_mode()) {
      comp_manager->removeAllCachePollutionsOf(decisions.top());
      comp_manager->clean_remain_comps_of(decisions.top());
    }
    VERBOSE_DEBUG_DO(cout << "now at dec lit: " << top_dec_lit() << " lev: " << dec_level() << " cnt:" <<  decisions.top().total_model_count() << endl);
  }
  VERBOSE_DEBUG_DO(print_comp_stack_info());
  VERBOSE_DEBUG_DO(cout << "DONE backw cleaning" << endl);
}

template<typename T>
void Counter<T>::check_trail([[maybe_unused]] bool check_entail) const {
  vector<uint32_t> num_decs_at_level(dec_level()+1, 0);
  bool entailment_fail = false;
  for(const auto& t: trail) {
    int32_t lev = var(t).decision_level;
    if (lev > dec_level()) {
      cout << "Too high decision level enqueued." << endl;
      assert(false);
    }
    if (var(t).ante.isNull() && lev > 0) {
      num_decs_at_level.at(lev)++;
      if (num_decs_at_level.at(lev) >= 2) {
        cout << "ERROR: Two or more of decs at level: " << lev << " trail follows." << endl;
        print_trail(false, false);
        assert(false);
      }
    }
    if (val(t) != T_TRI) {
      assert(false && "Trail is wrong, trail[val] is not TRUE");
    }
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
          cout << "Not implied by decisions above or at its level "
            << this_lev << " lit: " << t << " solver said: " << ret << endl;
          entailment_fail = true;
        }
      }
    }
#endif
  }
  if (entailment_fail) {
    cout << "Entailment fail." << endl;
    print_trail(false, false);
    print_dec_info();
  }
  assert(!entailment_fail);
}

template<typename T>
bool Counter<T>::is_implied(const vector<Lit>& cl) {
    assert(sat_solver);
    vector<CMSat::Lit> lits; lits.reserve(cl.size());
    for(const auto& l: cl) lits.push_back(CMSat::Lit(l.var()-1, l.sign()));
    debug_print("to check lits: " << lits);
    auto ret = sat_solver->solve(&lits);
    debug_print("Ret: " << ret);
    return ret == CMSat::l_False;
}


template<typename T>
void Counter<T>::check_implied(const vector<Lit>& cl) {
  bool implied = is_implied(cl);
  if (!implied) {
    cout << "ERROR, not implied" << endl;
    cout << "last dec lit: " << top_dec_lit() << endl;
    print_comp_stack_info();
    print_conflict_info();
    assert(false);
  }
}

template<typename T>
void Counter<T>::reduce_db_if_needed() {
  if (stats.conflicts > last_reducedb_confl+conf.reduce_db_everyN) {
    reduce_db();
    if (stats.cls_deleted_since_compaction > conf.consolidate_every_n && alloc->consolidate(this)) {
        stats.cls_deleted_since_compaction = 0;
    }
    last_reducedb_confl = stats.conflicts;
    last_reducedb_dec = stats.decisions;
  }
}

///out-ganak-7178163.pbs101-2/mc2023_track3_152.cnf.gz.out_d4
template<typename T>
RetState Counter<T>::resolve_conflict() {
  VERBOSE_DEBUG_DO(cout << "******" << __FUNCTION__<< " START" << endl);
  VERBOSE_DEBUG_DO(print_trail());

  create_uip_cl();
  if (uip_clause.size() == 1 && !exists_unit_cl_of(uip_clause[0]))
    unit_cls.push_back(uip_clause[0]);

  assert(uip_clause.front() != NOT_A_LIT);

  reduce_db_if_needed();
  VERBOSE_DEBUG_DO(print_conflict_info());

  stats.conflicts++;
  assert(decisions.top().remaining_comps_ofs() <= comp_manager->comp_stack_size());
  decisions.top().zero_out_branch_sol();
  decisions.top().mark_branch_unsat();

  VERBOSE_DEBUG_DO(cout << "backwards cleaning" << endl);
  VERBOSE_DEBUG_DO(print_comp_stack_info());
  int32_t backj = find_backtrack_level_of_learnt();
  int32_t lev_to_set = find_lev_to_set(backj);

  debug_print("backj: " << backj << " lev_to_set: " << lev_to_set);
  bool flipped_declit = (
      uip_clause[0].neg().var() == decisions.at(backj).var
           && lev_to_set+1 == backj);

  if (!flipped_declit || (sat_mode() && backj-1 >= sat_start_dec_level)) {
    debug_print("---- NOT FLIPPED DECLIT ----------");
    VERBOSE_DEBUG_DO(print_trail());
    VERBOSE_DEBUG_DO(print_conflict_info());
    debug_print("Not flipped. backj: " << backj << " lev_to_set: " << lev_to_set
      << " current lev: " << dec_level());
    go_back_to(backj-1);
    auto ant = add_uip_confl_cl(uip_clause);
    set_lit(uip_clause[0], lev_to_set, ant);
    VERBOSE_DEBUG_DO(print_trail());
    return RESOLVED;
  }

  assert(flipped_declit);
  VERBOSE_DEBUG_DO(cout << "after finding backj lev: " << backj << " lev_to_set: " << lev_to_set <<  endl);
  VERBOSE_DEBUG_DO(print_conflict_info());

  go_back_to(backj);
  VERBOSE_DEBUG_DO(print_conflict_info());
  debug_print("dec_level(): " << dec_level());

  Antecedent ant;
  assert(!uip_clause.empty());
  CHECK_IMPLIED_DO(check_implied(uip_clause));
  if (dec_level() > 0 && top_dec_lit().neg() == uip_clause[0]) {
    debug_print("FLIPPING. Setting reason the conflict cl");
    assert(var(uip_clause[0]).decision_level != -1);
    ant = add_uip_confl_cl(uip_clause);
    var(top_dec_lit()).ante = ant;
  }
  debug_print("Ant is :" << ant);
  debug_print("AFTER conflict, setup: ");
  VERBOSE_DEBUG_DO(print_conflict_info());
  debug_print("is right here? " << decisions.top().is_right_branch());

  decisions.top().zero_out_branch_sol();
  decisions.top().mark_branch_unsat();
  if (!sat_mode()) {
    comp_manager->removeAllCachePollutionsOf(decisions.top());
    decisions.top().reset_remain_comps();
  }

  if (decisions.top().is_right_branch()) {
    reactivate_comps_and_backtrack_trail(false);
    set_lit(uip_clause[0], lev_to_set, ant);
    qhead = std::min(qhead, var(uip_clause[0]).sublevel);
    if (!propagate()) return GO_AGAIN;

#ifdef VERBOSE_DEBUG
    cout << "FLIPPED Returning from resolveConflict() with:";
    print_conflict_info();
    print_trail(false); // we re-written the level above, so entailment
                        // may fail. when backtracking it'll be fine, though
    cout << "We have already counted this LEFT branch, so we backtrack now." << endl;
#endif
    return BACKTRACK;
  }

  if (dec_level() > 0 && !sat_mode()) {
    assert(decisions.top().remaining_comps_ofs() == comp_manager->comp_stack_size());
  }

  reactivate_comps_and_backtrack_trail(false);
  decisions.top().change_to_right_branch();
  set_lit(uip_clause[0], lev_to_set, ant);

#ifdef VERBOSE_DEBUG
  cout << "Returning from resolveConflict() with:";
  print_conflict_info();
  print_trail();
#endif

  return RESOLVED; // will ALWAYS propagate afterwards.
}

template<typename T>
inline void Counter<T>::get_maxlev_maxind(ClauseOfs ofs, int32_t& maxlev, uint32_t& maxind) {
  Clause& cl = *alloc->ptr(ofs);
  for(uint32_t i3 = 2; i3 < cl.sz; i3++) {
    Lit l = cl[i3];
    int32_t nlev = var(l).decision_level;
    VERBOSE_DEBUG_DO(cout << "i3: " << i3 << " l : " << l << " var(l).decision_level: "
        << var(l).decision_level << " maxlev: " << maxlev << endl);
    if (nlev > maxlev) {maxlev = nlev; maxind = i3;}
  }
}

template<typename T>
bool Counter<T>::propagate(bool out_of_order) {
  confl = Antecedent();
  debug_print("qhead in propagate(): " << qhead << " trail sz: " << trail.size() << " dec lev: " << dec_level() << " trail follows.");
  VERBOSE_DEBUG_DO(print_trail());
  for (; qhead < trail.size(); qhead++) {
    const Lit plit = trail[qhead].neg();
    const int32_t lev = var(plit).decision_level;
    bool lev_at_declev = false;

    if (!out_of_order) {
      if (decisions.size() <= 1) lev_at_declev = true;
      else if (var(top_dec_lit()).decision_level == lev) lev_at_declev = true;
    }
    debug_print("&&Propagating: " << plit.neg() << " qhead: " << qhead << " lev: " << lev);

    //Propagate bin clauses
    for (const auto& bincl : watches[plit].binaries) {
      const auto& l = bincl.lit();
      if (val(l) == F_TRI) {
        set_confl_state(plit, l);
        VERBOSE_DEBUG_DO(cout << "Bin confl. otherlit: " << l << endl);
      } else if (val(l) == X_TRI) {
        set_lit(l, lev, Antecedent(plit));
        VERBOSE_DEBUG_DO(cout << "Bin prop: " << l << " lev: " << lev << endl);
      }
    }

    //Propagate long clauses
    auto& ws = watches[plit].watch_list_;

#if 0
    cout << "prop-> will go through norm cl:" << endl;
    for(const auto& w: ws) {
      cout << "norm cl offsets: " << w.ofs << " cl: ";
      const auto ofs = w.ofs;
      for(Lit* c = beginOf(ofs); *c != NOT_A_LIT; c++) { cout << *c << " "; }
    }
    cout << "--> will do it now... " << endl;
#endif

    auto it2 = ws.begin();
    auto it = ws.begin();
    for (; it != ws.end(); it++) {
      if (is_true(it->blckLit)) { *it2++ = *it;
        debug_print("cl ofs: " << it->ofs << " blocked on lit: " << it->blckLit << " -> skipping");
        continue; }

      const auto ofs = it->ofs;
      Clause& c = *alloc->ptr(ofs);
      if (c[0] == plit) { std::swap(c[0], c[1]); }

#ifdef VERBOSE_DEBUG
      cout << "Prop Norm cl: " << ofs << " red: " << std::boolalpha << (bool)c.red << endl;
      for(const auto&l: c) {
        cout << "lit " << std::setw(6) << l
          << " lev: " << std::setw(4) << var(l).decision_level
          << " ante: " << std::setw(5) << std::left << var(l).ante
          << " val: " << lit_val_str(l) << endl;
      }
#endif

      assert(c[1] == plit);
      if (is_true(c[0])) {
        *it2++ = ClOffsBlckL(ofs, c[0]);
        continue;
      }

      uint32_t i = 2;
      for(; i < c.sz; i++) if (!is_false(c[i])) break;
      // either we found a free or satisfied lit
      if (i != c.sz) {
        c[1] = c[i];
        c[i] = plit;
        debug_print("New watch for cl: " << c[1]);
        watches[c[1]].add_cl(ofs, plit);
      } else {
        *it2++ = *it;
        if (val(c[0]) == F_TRI) {
          debug_print("Conflicting state from norm cl offs: " << ofs);
          set_confl_state(&c);
          it++;
          break;
        } else {
          assert(val(c[0]) == X_TRI);
          debug_print("prop long lev: " << lev << " dec_stack.get_lev : " << dec_level());
          if (lev_at_declev) {
            set_lit(c[0], lev, Antecedent(ofs));
            debug_print("Norm long prop: " << c[0] << " lev: " << lev);
          } else {
            int32_t maxlev = lev;
            uint32_t maxind = 1;
            get_maxlev_maxind(ofs, maxlev, maxind);
            if (maxind != 1) {
                std::swap(c[1], c[maxind]);
                it2--; // undo last watch
                watches[c[1]].add_cl(ofs, plit);
            }
            set_lit(c[0], maxlev, Antecedent(ofs));
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
  CHECK_PROPAGATED_DO(if (confl.isNull()) check_all_propagated_conflicted());
  debug_print("After propagate, qhead is: " << qhead << " conflict: " << std::boolalpha << !confl.isNull());
  return confl.isNull();
}

template<typename T>
const DataAndStatistics<T>& Counter<T>::get_stats() const {
  return stats; }

template<typename T>
bool Counter<T>::lit_redundant(Lit p, uint32_t abstract_levels) {
    debug_print(__func__ << " called");

    analyze_stack.clear();
    analyze_stack.push_back(p);

    Lit* c = nullptr;
    uint32_t size;
    size_t top = to_clear.size();
    while (!analyze_stack.empty()) {
      debug_print("At point in lit_redundant: " << analyze_stack.back());
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
              && (abst_level(p2.var()) & abstract_levels) != 0
          ) {
              debug_print("lit " << p2 << " OK");
              seen[p2.var()] = 1;
              analyze_stack.push_back(p2);
              to_clear.push_back(p2.var());
          } else {
              debug_print("lit " << p2 << " NOT OK");
              //Return to where we started before function executed
              for (size_t j = top; j < to_clear.size(); j++) seen[to_clear[j]] = 0;
              to_clear.resize(top);
              return false;
          }
        }
      }
    }
    debug_print("Returning OK from " << __func__);
    return true;
}

template<typename T>
uint32_t Counter<T>::abst_level(const uint32_t x) const {
  return ((uint32_t)1) << (var_data[x].decision_level & 31);
}

template<typename T>
void Counter<T>::recursive_cc_min()
{
  VERBOSE_DEBUG_DO(print_conflict_info());
  debug_print("recursive ccmin now.");
  uint32_t abstract_level = 0;
  for (size_t i = 1; i < uip_clause.size(); i++) {
    //(maintain an abstraction of levels involved in conflict)
    abstract_level |= abst_level(uip_clause[i].var());
  }

  size_t i, j;
  for (i = j = 1; i < uip_clause.size(); i++) {
    if (var(uip_clause[i]).ante.isNull()
      || !lit_redundant(uip_clause[i], abstract_level)
    ) {
      debug_print("ccmin -- keeping lit: " << uip_clause[i]);
      uip_clause[j++] = uip_clause[i];
    } else {
      debug_print("ccmin -- NOT keeping lit: " << uip_clause[i]);
    }
  }
  uip_clause.resize(j);
}

template<typename T>
void Counter<T>::minimize_uip_cl() {
  stats.uip_cls++;
  stats.orig_uip_lits += uip_clause.size();
  recursive_cc_min();
  for(const auto& c: to_clear) seen[c] = 0;
  to_clear.clear();

  CHECK_IMPLIED_DO(check_implied(uip_clause));
  tmp_cl_minim.clear();
  for(const auto& l:uip_clause) tmp_cl_minim.push_back(l);

  stats.uip_lits_ccmin+=tmp_cl_minim.size();
  if (stats.rem_lits_tried <= (200ULL*1000ULL) ||
      (stats.rem_lits_tried > (200ULL*1000ULL) &&
      ((double)stats.rem_lits_with_bins/(double)stats.rem_lits_tried > 3)))
    minimize_uip_cl_with_bins(tmp_cl_minim);
  stats.final_cl_sz+=tmp_cl_minim.size();
  uip_clause.clear();
  for(const auto& l: tmp_cl_minim) uip_clause.push_back(l);
  CHECK_IMPLIED_DO(check_implied(uip_clause));
}

template<typename T>
void Counter<T>::vivify_cls(vector<ClauseOfs>& cls) {
  stats.vivif_tried++;
  uint32_t j = 0;
  for(uint32_t i = 0; i < cls.size(); i++) {
    bool rem = false;
    auto& off = cls[i];
    if (v_tout > 0) {
      Clause& cl = *alloc->ptr(off);
      if (cl.vivifed == 0 &&
          (!cl.red || (cl.red && (cl.lbd <= lbd_cutoff || (cl.used && cl.total_used > conf.tot_used_cutoff_vivif)))))
        rem = vivify_cl(off);
    }
    if (!rem) cls[j++] = off;
  }

  // We didn't timeout, reset vivified flag.
  if (v_tout > 0) for(const auto& off: cls) alloc->ptr(off)->vivifed = 0;
  cls.resize(j);
}

template<typename T>
void Counter<T>::vivif_setup() {
  // Set up internals
  v_lev = 0;
  v_levs.clear();
  v_levs.resize(nVars()+1, -1);
  v_values.clear();
  v_values.resize(nVars()+1, X_TRI);
  v_qhead = 0;

  // Set up units
  v_trail.clear();
  for(const auto& l: trail) if (var(l).decision_level == 0) v_enqueue(l);
  for(const auto& l: unit_cls) if (v_val(l) == X_TRI) v_enqueue(l);
  bool ret = v_propagate();
  assert(ret);

}

template<typename T>
void Counter<T>::vivify_all(bool force, bool only_irred) {
  if (!force && last_confl_vivif + conf.vivif_every > stats.conflicts) return;

  CHECK_PROPAGATED_DO(check_all_propagated_conflicted());
  double my_time = cpu_time();
  uint64_t last_vivif_lit_rem = stats.vivif_lit_rem;
  uint64_t last_vivif_cl_minim = stats.vivif_cl_minim;
  auto last_vivif_cl_tried = stats.vivif_tried_cl;

  // Sanity check here.
  last_confl_vivif = stats.conflicts;

  // Backup 1st&2nd watch + block lit
  off_to_lit12.clear();
  for(const auto& off: long_irred_cls) {
    const Clause& cl = *alloc->ptr(off);
    bool curr_prop = currently_propagating_cl(cl);
    off_to_lit12[off] = SavedCl(cl[0], cl[1], curr_prop);
  }
  for(const auto& off: long_red_cls) {
    const Clause& cl = *alloc->ptr(off);
    bool curr_prop = currently_propagating_cl(cl);
    off_to_lit12[off] = SavedCl(cl[0], cl[1], curr_prop);
  }
  all_lits(i) {
    Lit lit(i/2, i%2);
    for(const auto& ws: watches[lit].watch_list_) {
      auto it = off_to_lit12.find(ws.ofs);
      assert(it != off_to_lit12.end());
      if (lit == it->second.first) it->second.blk1 = ws.blckLit;
      else if (lit == it->second.second) it->second.blk2 = ws.blckLit;
      else assert(false);
    }
  }

  vivif_setup();
  verb_print(2, "[vivif] setup. T: " << setprecision(2) << (cpu_time()-my_time));

  // Vivify clauses
  v_tout = conf.vivif_mult*2LL*1000LL*1000LL;
  if (force) v_tout *= 50;
  vivify_cls(long_irred_cls);
  bool tout_irred = (v_tout <= 0);
  verb_print(2, "[vivif] irred vivif remain: " << v_tout/1000 << "K T: " << (cpu_time()-my_time));

  bool tout_red = false;
  if (!only_irred) {
    v_tout = conf.vivif_mult*20LL*1000LL*1000LL;
    vivify_cls(long_red_cls);
    verb_print(2, "[vivif] red vivif remain: " << v_tout/1000 << "K T: " << (cpu_time()-my_time));
    tout_red = (v_tout <= 0);
  }

  // Restore
  for(auto& ws: watches) ws.watch_list_.clear();
  if (!decisions.empty()) {
    for(const auto& off: long_irred_cls) v_cl_repair(off);
    for(const auto& off: long_red_cls) v_cl_repair(off);
  } else {
    // Move all 0-level stuff to unit_cls
    for(const auto& l: v_trail) {
      if (val(l) == X_TRI) {
        set_lit(l, 0);
        if (!exists_unit_cl_of(l)) unit_cls.push_back(l);
      }
      assert(val(l) != F_TRI); // it would be UNSAT
    }
    bool ret2 = propagate();
    assert(ret2);
    v_cl_toplevel_repair(long_irred_cls);
    v_cl_toplevel_repair(long_red_cls);
  }
  off_to_lit12.clear();
  verb_print(2, "[vivif] finished."
      << " cl tried: " << (stats.vivif_tried_cl - last_vivif_cl_tried)
      << " cl minim: " << (stats.vivif_cl_minim - last_vivif_cl_minim)
      << " lit rem: " << (stats.vivif_lit_rem - last_vivif_lit_rem)
      << " force: " << (int)force
      << " tout-irred: " << (int)tout_irred
      << " tout-red: " << (int)tout_red
      << " T: " << (cpu_time()-my_time));
  CHECK_PROPAGATED_DO(check_all_propagated_conflicted());
}

template<typename T>
template<class T2>
bool Counter<T>::v_satisfied(const T2& lits) {
  for(auto& l: lits) if (v_val(l) == T_TRI) return true;
  return false;
}

template<typename T>
template<class T2>
bool Counter<T>::v_unsat(const T2& lits) {
  for(auto& l: lits) if (v_val(l) == T_TRI || v_val(l) == X_TRI) return false;
  return true;
}

template<typename T>
void Counter<T>::v_shrink(Clause& cl) {
  uint32_t j = 0;
  for(uint32_t i = 0; i < cl.size(); i++) {
    if (v_val(cl[i]) == F_TRI) continue;
    cl[j++] = cl[i];
  }
  cl.resize(j);
}

template<typename T>
void Counter<T>::v_cl_toplevel_repair(vector<ClauseOfs>& offs) {
  uint32_t j = 0;
  for(uint32_t i = 0; i < offs.size(); i++) {
    Clause* cl = alloc->ptr(offs[i]);
    assert(!v_unsat(*cl));
    if (v_satisfied(*cl)) {alloc->clause_free(cl);continue;}
    v_shrink(*cl);
    assert(cl->size() >= 2);
    if (cl->size() == 2) {
      add_bin_cl((*cl)[0], (*cl)[1], cl->red);
      alloc->clause_free(cl);
      continue;
    }
    attach_cl(offs[i], (*cl));
    offs[j++] = offs[i];
  }
  offs.resize(j);
}

template<typename T>
void Counter<T>::v_cl_repair(ClauseOfs off) {
  Clause& cl = *alloc->ptr(off);
  auto& offs = off_to_lit12[off];

  if (offs.currently_propagating) {
    // Move 1st & 2nd literal to position
    auto at = std::find(cl.begin(), cl.end(), offs.first);
    assert(at != cl.end());
    std::swap(cl[0], *at);
    at = std::find(cl.begin(), cl.end(), offs.second);
    assert(at != cl.end());
    std::swap(cl[1], *at);

    watches[cl[0]].add_cl(off, offs.blk1);
    watches[cl[1]].add_cl(off, offs.blk2);
    return;
  }

  std::sort(cl.begin(), cl.end(),
    [=, this](const Lit& a, const Lit& b) {
      // undef must be at the beginning.
      if (val(a) == X_TRI && val(b) != X_TRI) return true;
      if (val(b) == X_TRI && val(a) != X_TRI) return false;
      if (var(a).decision_level == 0) return false;
      if (var(b).decision_level == 0) return true;

      // Largest sublevel first
      return var(a).sublevel > var(b).sublevel;
    });

  int32_t t_at = -1;
  for(uint32_t i = 2; i < cl.size(); i++) {
    if (val(cl[i]) == T_TRI) {t_at = i; break;}
  }

  debug_print("Vivified cl off: " << off);
  VERBOSE_DEBUG_DO(print_cl(cl));
  Lit blk = (t_at == -1) ? cl[cl.sz/2] : cl[t_at];
  watches[cl[0]].add_cl(off, blk);
  watches[cl[1]].add_cl(off, blk);
}

// We could have removed a TRUE. This may be an issue.
template<typename T>
void Counter<T>::v_fix_watch(Clause& cl, uint32_t i) {
  if (val(cl[i]) == X_TRI || val(cl[i]) == T_TRI) return;
  auto off = alloc->get_offset(&cl);
  watches[cl[i]].del_c(off);
  uint32_t i2 = 2;
  for(; i2 < cl.size(); i2++) if (val(cl[i2]) == X_TRI || val(cl[i2]) == T_TRI) break;
  /* print_cl(cl); */
  assert(i2 != cl.size());
  std::swap(cl[i], cl[i2]);
  watches[cl[i]].add_cl(off, cl[cl.sz/2]);
}

template<typename T>
void Counter<T>::v_new_lev() {
  assert(v_lev == 0);
  v_lev++;
  v_backtrack_to = v_trail.size();
}

template<typename T>
void Counter<T>::v_unset_lit(const Lit l) {
  debug_print("v_unset_lit: " << l);
  assert(v_levs[l.var()] == 1);
  v_levs[l.var()] = -1;
  v_values[l] = X_TRI;
  v_values[l.neg()] = X_TRI;
}

template<typename T>
void Counter<T>::v_backtrack() {
  assert(v_lev == 1);
  for(uint32_t i = v_backtrack_to; i < v_trail.size(); i++) {
    const auto& l = v_trail[i];
    v_unset_lit(l);
  }
  v_trail.resize(v_backtrack_to);
  v_lev = 0;
  v_qhead = v_trail.size();
}

template<typename T>
template<class T2>
bool Counter<T>::v_cl_satisfied(const T2& cl) const {
  for(const auto&l : cl) {
    if (v_val(l) == T_TRI) return true;
  }
  return false;
}

template<typename T>
template<class T2>
bool Counter<T>::propagation_correctness_of_vivified(const T2& cl) const {
  uint32_t num_t = 0;
  int32_t t_lev = -1;
  int32_t maxlev_f = -1;
  // Check if it should_have_propagated_earlier
  for(const auto&l: cl) {
    if (val(l) == T_TRI) {
      num_t++;
      if (num_t >= 2) return true;
      t_lev = var(l).decision_level;
    }
    if (val(l) == X_TRI) return true;
    if (val(l) == F_TRI) {
      maxlev_f = std::max(maxlev_f, var(l).decision_level);
    }
  }

  // Should have propagated at level maxlev_f -- but it only got set TRUE at t_lev!
  if (maxlev_f < t_lev) return false;

  // Have to find a FALSE at the level the TRUE is at
  for(const auto&l: cl) {
    if (val(l) == F_TRI) {
      if (var(l).decision_level == t_lev) return true;
    }
  }
  return false;
}

// Returns TRUE if we can remove the clause
template<typename T>
bool Counter<T>::vivify_cl(const ClauseOfs off) {
  SLOW_DEBUG_DO(for(auto& l: seen) assert(l == 0));
  bool fun_ret = false;
  Clause& cl = *alloc->ptr(off);

  /* cout << "orig CL: " << endl; v_print_cl(cl); */
  auto it = off_to_lit12.find(off);
  if (it->second.currently_propagating) return false;
  v_tmp.clear();
  v_tmp2.clear();
  for(const auto&l: cl) v_tmp2.push_back(l);
  std::shuffle(v_tmp2.begin(), v_tmp2.end(), mtrand);

  // Swap to 1st & 2nd the two original 1st & 2nd
  auto sw = std::find(v_tmp2.begin(), v_tmp2.end(), it->second.first);
  std::swap(*sw, v_tmp2[0]);
  sw = std::find(v_tmp2.begin(), v_tmp2.end(), it->second.second);
  std::swap(*sw, v_tmp2[1]);
  if (v_val(v_tmp2[0]) != X_TRI || v_val(v_tmp2[1]) != X_TRI) return false;

  v_new_lev();
  cl.vivifed = 1;
  stats.vivif_tried_cl++;

  debug_print("vivifying cl offs: " << off);
  VERBOSE_DEBUG_DO(print_cl(cl));

  v_tmp.push_back(v_tmp2[0]);
  v_tmp.push_back(v_tmp2[1]);
  v_enqueue(v_tmp2[0].neg());
  v_enqueue(v_tmp2[1].neg());
  bool ret = v_propagate();
  if (ret) {
    for(uint32_t i = 2; i < v_tmp2.size(); i++) {
      const auto& l = v_tmp2[i];
      debug_print("Vivif lit l: " << l << " val: " << val_to_str(v_val(l)));
      if (v_val(l) == T_TRI) {v_tmp.push_back(l);break;}
      if (v_val(l) == F_TRI) continue;
      v_tmp.push_back(l);
      v_enqueue(l.neg());
      ret = v_propagate();
      if (!ret) {
        debug_print("vivif ret FALSE, exiting");
        break;
      }
    }
  }
  v_backtrack();
  VERBOSE_DEBUG_DO(cout << "new vivified CL offs: " << off << endl; print_cl(v_tmp));
  uip_clause.clear();
  CHECK_IMPLIED_DO(check_implied(v_tmp));
  for(const auto&l: v_tmp) seen[l.raw()] = 1;

  uint32_t removable = 0;
  for(const auto& l: v_tmp2) {
    if (seen[l.raw()] == 0) {
      removable++;
    } else {
      to_clear.push_back(l.raw());
    }
  }

  if (removable != 0 &&
      // TODO once chronological backtracking works, we can have level-0 stuff. Not now.
      //      so we must skip this
      !propagating_cl(v_tmp) && !conflicting_cl(v_tmp) &&
      propagation_correctness_of_vivified(v_tmp)) {
    watches[cl[0]].del_c(off);
    watches[cl[1]].del_c(off);
    VERBOSE_DEBUG_DO(cout << "orig CL: " << endl; v_print_cl(cl));
    stats.vivif_cl_minim++;
    stats.vivif_lit_rem += removable;
    for(uint32_t i = 0; i < v_tmp.size(); i++) cl[i] = v_tmp[i];
    cl.resize(v_tmp.size());
    assert(cl.sz >= 2);
    VERBOSE_DEBUG_DO(cout << "vivified CL: " << endl; v_print_cl(cl));

    std::sort(cl.begin(), cl.end(), [=, this](const Lit l1, const Lit l2) {
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
        // cannot propagate!
        assert(false);
      }
      for(uint32_t v = 1; v < var_data.size(); v++) {
        auto& vdat = var_data[v];
        if (vdat.ante.isAClause() && vdat.ante.as_cl() == off) {
          assert(v == cl[0].var() || v == cl[1].var());
          Lit other_lit = (v == cl[0].var()) ? cl[1] : cl[0];
          vdat.ante = Antecedent(other_lit);
        }
      }
      alloc->clause_free(off);
      fun_ret = true;
    } else {
      watches[cl[0]].add_cl(off, cl[cl.sz/2]);
      watches[cl[1]].add_cl(off, cl[cl.sz/2]);
      if (!v_cl_satisfied(cl) && v_val(cl[0]) == X_TRI && v_val(cl[1]) != X_TRI) {
        //cannot propagate!
        assert(v_val(cl[1]) == F_TRI);
        assert(false);
      }
      cl.update_lbd(cl.sz); // we may be smaller than LBD
    }
    ret = v_propagate();
    assert(ret);
  } else {
    VERBOSE_DEBUG_DO(cout << "Can't vivify." << endl);
  }

  for(const auto& l: to_clear) seen[l] = 0;
  to_clear.clear();
  return fun_ret;
}

template<typename T>
TriValue Counter<T>::v_val(const Lit l) const {
  return v_values[l];
}

template<typename T>
void Counter<T>::v_enqueue(const Lit l) {
  debug_print("v-enq: " << l << " lev: " << v_lev);
  assert(v_val(l) == X_TRI);
  v_levs[l.var()] = v_lev;
  v_trail.push_back(l);
  v_values[l] = T_TRI;
  v_values[l.neg()] = F_TRI;
}

template<typename T>
bool Counter<T>::v_propagate() {
  bool ret = true;
  for (; v_qhead < v_trail.size(); v_qhead++) {
    const Lit plit = v_trail[v_qhead].neg();

    //Propagate bin clauses
    const auto& wsbin = watches[plit].binaries;
    v_tout-=wsbin.size()/2;
    for (const auto& bincl : wsbin) {
      const auto& l = bincl.lit();
      if (v_val(l) == F_TRI) {
        debug_print("v Conflict from bin.");
        return false;
      } else if (v_val(l) == X_TRI) {
        v_enqueue(l);
        debug_print("v Bin prop: " << l);
      }
    }

    //Propagate long clauses
    auto& ws = watches[plit].watch_list_;
    v_tout-=ws.size()/2;

#ifdef VERBOSE_DEBUG
    cout << "v prop-> will go through norm cl:" << endl;
    for(const auto& w: ws) {
      cout << "norm cl offsets: " << w.ofs << " cl: ";
      const auto ofs = w.ofs;
      Clause& c = *alloc->ptr(ofs);
      cout << c << endl;
    }
    cout << "--> will do it now... " << endl;
#endif

    auto it2 = ws.begin();
    auto it = ws.begin();
    for (; it != ws.end(); it++) {
      if (v_val(it->blckLit) == T_TRI) { *it2++ = *it; continue; }

      const auto ofs = it->ofs;
      Clause& c = *alloc->ptr(ofs);
      if (c[0] == plit) { std::swap(c[0], c[1]); }
      v_tout--;

#ifdef VERBOSE_DEBUG
      cout << "v Prop Norm cl: " << ofs << endl;
      for(const auto&l: c) {
        cout << "lit " << std::setw(6) << l
          << " lev: " << std::setw(4) << var(l).decision_level
          << " ante: " << std::setw(5) << std::left << var(l).ante
          << " val: " << lit_val_str(l) << endl;
      }
#endif

      assert(c[1] == plit);
      if (v_val(c[0]) == T_TRI) {
        *it2++ = ClOffsBlckL(ofs, c[0]);
        continue;
      }

      uint32_t i = 2;
      for(; i < c.sz; i++) if (v_val(c[i]) != F_TRI) break;
      // either we found a free or satisfied lit
      if (i != c.sz) {
        c[1] = c[i];
        c[i] = plit;
        debug_print("v New watch for cl: " << c[1]);
        watches[c[1]].add_cl(ofs, plit);
      } else {
        *it2++ = *it;
        if (v_val(c[0]) == F_TRI) {
          debug_print("v Conflicting state from norm cl offs: " << ofs);
          ret = false;
          it++;
          break;
        } else {
          assert(v_val(c[0]) == X_TRI);
          debug_print("v prop long");
          v_enqueue(c[0]);
        }
      }
    }
    while(it != ws.end()) *it2++ = *it++;
    ws.resize(it2-ws.begin());
    if (!ret) break;
  }
  debug_print("v After propagate, v_qhead is: " << v_qhead << " returning: " << ret);
  return ret;
}

template<typename T>
void Counter<T>::fill_cl(const Antecedent& ante, Lit*& c, uint32_t& size, Lit p) const {
  if (ante.isAClause()) {
    Clause* cl = alloc->ptr(ante.as_cl());
    c = cl->data();
    size = cl->sz;
  } else if (ante.isALit()) {
    //Binary
    tmp_lit.resize(2);
    c = tmp_lit.data();
    if (p == NOT_A_LIT) c[0] = confl_lit;
    else c[0] = p;
    c[1] = ante.as_lit();
    size = 2;
  } else {assert(false && "Should never be a decision");}
}

template<typename T>
typename Counter<T>::ConflictData Counter<T>::find_conflict_level(Lit p) {
  ConflictData data;
  Lit* c;
  uint32_t size;
  fill_cl(confl, c, size, p);
  VERBOSE_DEBUG_DO(cout << "CL in find_conflict_level " << confl << " : " << endl;print_cl(c, size));
  data.nHighestLevel = var(c[0]).decision_level;
  if (data.nHighestLevel == dec_level() && var(c[1]).decision_level == dec_level())
    return data;

  int highest_id = 0;
  data.bOnlyOneLitFromHighest = true;
  // find the largest decision level in the clause
  for (uint32_t i = 1; i < size; ++i) {
    int32_t lev = var(c[i]).decision_level;
    if (lev > data.nHighestLevel) {
      highest_id = i;
      data.nHighestLevel = lev;
      data.bOnlyOneLitFromHighest = true;
    } else if (lev == data.nHighestLevel && data.bOnlyOneLitFromHighest == true) {
      data.bOnlyOneLitFromHighest = false;
    }
  }

  // fixing clause & watchlist
  if (highest_id != 0 && confl.isAClause()) {
    Clause& cl = *alloc->ptr(confl.as_cl());
    std::swap(cl[0], cl[highest_id]); // swap to position 1, since we'll swap 1&0 in recordLastUIPClauses
    debug_print("SWAPPED");
    VERBOSE_DEBUG_DO(print_cl(cl.data(), cl.size()));
    if (highest_id > 1 && size > 2) {
      ClauseOfs off = confl.as_cl();
      watches[cl[highest_id]].del_c(off);
      watches[c[0]].add_cl(off, c[1]);
    }
  }
  return data;
}

template<typename T>
void Counter<T>::create_uip_cl() {
  assert(to_clear.empty());

  uip_clause.clear();
  uip_clause.push_back(Lit(0, false));
  Lit p = NOT_A_LIT;

  SLOW_DEBUG_DO(for(const auto& t:seen) assert(t == 0););
  VERBOSE_DEBUG_DO(print_dec_info());
  int32_t n_dec_level = std::numeric_limits<int32_t>::min();

  Lit* c;
  uint32_t size;
  VERBOSE_DEBUG_DO(cout << "Doing loop:" << endl);
  int32_t index = trail.size()-1;
  uint32_t path_c = 0;
  do {
    fill_cl(confl, c, size, p);
    if (confl.isAClause()) {
      Clause& cl = *alloc->ptr(confl.as_cl());
      if (cl.red && cl.lbd > lbd_cutoff) {
        cl.set_used();
        /* cl.update_lbd(calc_lbd(cl)); */
      }
    }
    if (p == NOT_A_LIT) {
      if (var(c[0]).decision_level < var(c[1]).decision_level) std::swap(c[0], c[1]);
      n_dec_level = var(c[0]).decision_level;
      SLOW_DEBUG_DO(check_cl_unsat(c, size));
    }
    VERBOSE_DEBUG_DO(cout << "next cl: " << endl;print_cl(c, size));
    VERBOSE_DEBUG_DO(cout << "n_dec_level: " <<  n_dec_level << endl);

    VERBOSE_DEBUG_DO(cout << "For loop." << endl);
    for(uint32_t j = ((p == NOT_A_LIT) ? 0 : 1); j < size ;j++) {
      Lit q = c[j];
      if (!seen[q.var()] && var(q).decision_level > 0){
        inc_act(q);
        seen[q.var()] = 1;
        to_clear.push_back(q.var());
#ifdef VERBOSE_DEBUG
        cout << std::setw(5) << q
          << " lev: " << std::setw(3) << var(q).decision_level
          << " ante: " << std::setw(8) << var(q).ante
          << " val : " << std::setw(7) << lit_val_str(q)
          << endl;
#endif
        if (var(q).decision_level >= n_dec_level) {
          path_c++;
          VERBOSE_DEBUG_DO(cout << "pathc inc." << endl);
        } else {
          uip_clause.push_back(q);
          VERBOSE_DEBUG_DO(cout << "added to cl." << endl);
        }
      }
    }
    VERBOSE_DEBUG_DO(cout << "PathC: " << path_c << endl);

    do {
      while (!seen[trail[index--].var()]) { SLOW_DEBUG_DO(assert(index >= 0));}
      p = trail[index+1];
      assert(p != NOT_A_LIT);
#ifdef VERBOSE_DEBUG
      cout << "going back on trail: " << std::setw(5) << p
        << " lev: " << std::setw(3) << var(p).decision_level
        << " ante: " << std::setw(8) << var(p).ante
        << " val : " << std::setw(7) << lit_val_str(p)
        << endl;
#endif
    } while(var(trail[index+1]).decision_level < n_dec_level);
    VERBOSE_DEBUG_DO(cout << "Next p: " << p << endl);
    confl = var(p).ante;
    seen[p.var()] = 0;
    path_c--;
  } while (path_c > 0);
  assert(path_c == 0);
  uip_clause[0] = p.neg();
  VERBOSE_DEBUG_DO(cout << "UIP cl: " << endl; print_cl(uip_clause.data(), uip_clause.size()));
  CHECK_IMPLIED_DO(check_implied(uip_clause));
  minimize_uip_cl();
  SLOW_DEBUG_DO(for(const auto& s: seen) assert(s == 0));
  VERBOSE_DEBUG_DO(cout << __FUNCTION__ << " finished");
}

template<typename T>
bool Counter<T>::check_watchlists() const {
  bool ret = true;
#if 0
  // All watchlists
  cout << "All watchlists: " << endl;
  all_lits(i) {
    Lit lit = Lit(i/2, i%2);
    cout << "->Watchlist for lit " << lit << " (val: " << lit_val_str(lit) << ") " << endl;
    auto& ws = watches[lit].watch_list_;
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
    auto& ws = watches[lit].watch_list_;
    for(const auto& w: ws) {
      const auto ofs = w.ofs;
      uint32_t num_unk = 0;
      bool sat = false;
      for(const auto& l: *alloc->ptr(ofs)) {
        if (is_unknown(l)) num_unk++;
        if (is_true(l)) sat = true;
      }
      if (!sat && num_unk >=2 && !is_unknown(lit)) {
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
    for(const auto& ws: watches[lit].watch_list_) {
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
  for(const auto& off: long_irred_cls) check_attach(off);
  for(const auto& off: long_red_cls) check_attach(off);
  if (!off_att_num.empty()) {
    cout << "ERROR: The following clauses are attached but are NOT in longRed/longIrred clauses" << endl;
    for(const auto& p: off_att_num) {
      cout << "Offset: " << p.first << endl;
      print_cl(*alloc->ptr(p.first));
    }
  }
  return ret;
}

// BEWARE! It sorts clauses, hence invalidates a lot of invariants about
// propagations
template<typename T>
void Counter<T>::attach_occ(vector<ClauseOfs>& cls, bool sort_and_clear) {
  for(const auto& off: cls) {
    occ_cls.push_back(off);
    Clause& cl = *alloc->ptr(off);
    if (sort_and_clear) std::sort(cl.begin(), cl.end());
    auto abs = calc_abstr(cl);
    for(const auto& l: cl) {
      SLOW_DEBUG_DO(assert(l.var() <= nVars()));
      SLOW_DEBUG_DO(assert(occ.size() > l.raw()));
      occ[l.raw()].push_back(OffAbs(off, abs));
    }
  }
  if (sort_and_clear) cls.clear();
}

template<typename T>
void Counter<T>::backw_susume_cl(ClauseOfs off) {
  Clause& cl = *alloc->ptr(off);
  uint32_t abs = calc_abstr(cl);
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
    if (!subset_abstr(abs, check.abs)) continue;
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
      alloc->clause_free(&check_cl);
      stats.subsumed_long_red_cls+=cl.red;
      stats.subsumed_long_irred_cls+=!cl.red;
    }
  }
}

template<typename T>
void Counter<T>::backw_susume_cl_with_bin(BinClSub& cl) {
  uint32_t abs = calc_abstr(cl);
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
    if (!subset_abstr(abs, check.abs)) continue;
    Clause& check_cl = *alloc->ptr(check.off);
    if (check_cl.freed) continue;
    if (subset(cl, check_cl)) {
      if (cl.red && !check_cl.red) cl.red = false;
      debug_print( "Subsumed cl: " << check_cl << endl
                << "->by cl    : " << cl);
      alloc->clause_free(&check_cl);
      stats.subsumed_long_red_cls+=cl.red;
      stats.subsumed_long_irred_cls+=!cl.red;
    }
  }
}

template<typename T>
void Counter<T>::toplevel_full_probe() {
  SLOW_DEBUG_DO(for(auto& l: seen) assert(l == 0));
  assert(to_clear.empty());
  assert(bothprop_toset.empty());

  double my_time = cpu_time();
  auto old_probe = stats.toplevel_probe_fail;
  auto old_bprop = stats.toplevel_bothprop_fail;
  stats.toplevel_probe_runs++;
  assert(dec_level() == 0);

  SLOW_DEBUG_DO(for(const auto&c: seen) assert(c == 0));
  for(uint32_t i = 1; i <= nVars(); i++) {
    Lit l = Lit(i, 0);
    if (val(l) != X_TRI) continue;

    decisions.push_back(StackLevel<T>(1,2));
    decisions.back().var = l.var();
    set_lit(l, 1);
    uint32_t trail_before = trail.size();
    bool ret = propagate();
    if (ret) {
      for(uint32_t i2 = trail_before; i2 < trail.size(); i2++) {
        Lit l2 = trail[i2];
        seen[l2.raw()] = 1;
        to_clear.push_back(l2.raw());
      }
    }
    reactivate_comps_and_backtrack_trail();
    decisions.pop_back();
    if (!ret) {
      clear_toclear_seen();
      set_lit(l.neg(), 0);
      ret = propagate();
      assert(ret && "we are never UNSAT");
      stats.toplevel_probe_fail++;
      continue;
    }

    // Negation
    assert(dec_level() == 0);
    decisions.push_back(StackLevel<T>(1,2));
    decisions.back().var = l.var();
    set_lit(l.neg(), 1);

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
    decisions.pop_back();
    if (!ret) {
      clear_toclear_seen();
      set_lit(l, 0);
      ret = propagate();
      assert(ret && "we are never UNSAT");
      stats.toplevel_probe_fail++;
      continue;
    }

    clear_toclear_seen();
    for(const auto& x: bothprop_toset) {
      set_lit(x, 0);
    }
    bothprop_toset.clear();
    ret = propagate();
    assert(ret && "we are never UNSAT");
  }

  SLOW_DEBUG_DO(for(const auto&c: seen) assert(c == 0));
  verb_print(2, "[top-probe] "
      << " failed: " << (stats.toplevel_probe_fail - old_probe)
      << " bprop: " << (stats.toplevel_bothprop_fail - old_bprop)
      << " T: " << (cpu_time()-my_time));
}

template<typename T>
void Counter<T>::subsume_all() {
  assert(dec_level() == 0);
  assert(occ.empty());
  assert(occ_cls.empty());

  // setup
  double my_time = cpu_time();
  auto old_subsumed_long_irred_cls = stats.subsumed_long_irred_cls;
  auto old_subsumed_long_red_cls = stats.subsumed_long_red_cls;
  auto old_subsumed_bin_irred_cls = stats.subsumed_bin_irred_cls;
  auto old_subsumed_bin_red_cls = stats.subsumed_bin_red_cls;
  stats.subsume_runs++;
  occ.resize((nVars()+1)*2);
  attach_occ(long_irred_cls, true); // beware-- sorts the clauses, invalidates prop invariants
  attach_occ(long_red_cls, true);

  // Detach everything
  for(auto& ws: watches) ws.watch_list_.clear();

  // Binary clauses
  vector<BinClSub> bin_cls;
  all_lits(i) {
    Lit lit = Lit(i/2, i%2);
    for(const auto& l2: watches[lit].binaries) {
      if (l2.lit() < lit) continue;
      assert(lit < l2.lit());
      bin_cls.push_back(BinClSub(lit, l2.lit(), l2.red()));
    }
    watches[lit].binaries.clear();
  }
  std::sort(bin_cls.begin(), bin_cls.end());
  uint32_t j = 0;
  for(uint32_t i = 1; i < bin_cls.size(); i++) {
    if (bin_cls[i] == bin_cls[j]) {
      stats.subsumed_bin_red_cls+= bin_cls[i].red;
      stats.subsumed_bin_irred_cls+= !bin_cls[i].red;
      continue;
  }
    if (bin_cls[i].lit[0] == bin_cls[j].lit[0]
       && bin_cls[i].lit[1] == bin_cls[j].lit[1]) {
      // ordering ensures IRRED is first
      stats.subsumed_bin_red_cls+= bin_cls[i].red;
      stats.subsumed_bin_irred_cls+= !bin_cls[i].red;
      continue;
    }
    j++;
    bin_cls[j] = bin_cls[i];
  }
  j++;
  if (!bin_cls.empty()) bin_cls.resize(j);
  for(auto& b: bin_cls) backw_susume_cl_with_bin(b);

  // Long clauses
  std::shuffle(occ_cls.begin(), occ_cls.end(), mtrand);
  for(const auto& off: occ_cls) {
    Clause* cl = alloc->ptr(off);
    if (cl->freed) continue;
    backw_susume_cl(off);
  }

  // Cleanup
  for(const auto& b: bin_cls) add_bin_cl(b[0], b[1], b.red);
  for(const auto& off: occ_cls) {
    Clause& cl = *alloc->ptr(off);
    if (cl.freed) continue;
    if (cl.red) long_red_cls.push_back(off);
    else long_irred_cls.push_back(off);

    std::sort(cl.begin(), cl.end(),
      [=, this](const Lit& a, const Lit& b) {
        // undef must be at the beginning.
        if (val(a) == X_TRI && val(b) != X_TRI) return true;
        if (val(b) == X_TRI && val(a) != X_TRI) return false;

        // Undef first as long as it's the same declevel
        if (var(a).decision_level == var(b).decision_level) {
          if(val(a) != val(b)) return val(a) == X_TRI;
          return false;
        }

        // Largest declevel first
        return var(a).decision_level > var(b).decision_level;
      });
    attach_cl(off, cl);
  }
  occ.clear();
  occ_cls.clear();
  verb_print(2, "[sub] "
      << " bin-irred-cls: " << stats.subsumed_bin_irred_cls - old_subsumed_bin_irred_cls
      << " bin-red-cls: " << stats.subsumed_bin_red_cls - old_subsumed_bin_red_cls
      << " long-irred-cls: " << stats.subsumed_long_irred_cls - old_subsumed_long_irred_cls
      << " long-red-cls: " << stats.subsumed_long_red_cls - old_subsumed_long_red_cls
      << " T: " << (cpu_time() - my_time));
}

template<typename T>
void Counter<T>::vsads_readjust() {
  if (stats.decisions % conf.vsads_readjust_every == 0)
    for(auto& w: watches) w.activity *= 0.5;
}

// At this point, the problem is either SAT or UNSAT, we only care about 1 or 0,
// because ONLY non-independent variables remain
template<typename T>
bool Counter<T>::use_sat_solver(RetState& state) {
  assert(!isindependent);
  assert(order_heap.empty());
  assert(!decisions.empty());
  assert(!sat_mode());
  stats.sat_called++;
  auto conflicts_before = stats.conflicts;

  debug_print("Entering SAT mode. Declev: " << dec_level() << " trail follows.");
  VERBOSE_DEBUG_DO(print_trail());
  bool sat = false;
  decisions.push_back(StackLevel<T>(decisions.top().curr_remain_comp(),
        comp_manager->comp_stack_size()));
  sat_start_dec_level = dec_level();

  // Fill up order heap
  all_vars_in_comp(comp_manager->get_super_comp(decisions.top()), it) {
    debug_print("checking var to put in order_heap: " << *it << " val: "
      << val_to_str(val(*it)));
    if (val(*it) != X_TRI) continue;
    if (*it < opt_indep_support_end)
      assert(*it >= indep_support_end && "only optional indep or non-indep remains");
    order_heap.insert(*it);
  }
  debug_print("Order heap size: " << order_heap.size());
  decisions.top().var = 0;
  auto old_sublev = trail.size();

  // the SAT loop
  auto orig_confl = stats.conflicts;
  auto last_restart = 0;
  uint32_t num_rst = 0;
  while(true) {
    uint32_t d;
    do {
      if (order_heap.empty()) {d = 0; break;}
      d = order_heap.removeMin();
    } while (val(d) != X_TRI);
    if (d == 0) {
      debug_print("SAT mode found a solution. dec lev: " << dec_level());
      SLOW_DEBUG_DO(check_sat_solution());
      sat = true;
      break;
    }
    stats.decisions++;
    vsads_readjust();
    assert(val(d) == X_TRI);
    Lit l(d, var(d).last_polarity);
    if (decisions.top().var != 0) {
      decisions.push_back(StackLevel<T>(1,2));
    }
    decisions.back().var = l.var();
    set_lit(l, dec_level());

    while (!propagate()) {
      start1:
      if (chrono_work()) continue;
      state = resolve_conflict();
      if (state == GO_AGAIN) goto start1;
      if (state == BACKTRACK) break;
    }
    if (state == BACKTRACK) goto end;
    assert(state != GO_AGAIN);
    if (dec_level() < sat_start_dec_level) { goto end; }
    const auto sat_confl = stats.conflicts -orig_confl;
    if (sat_confl-last_restart >= luby(2, num_rst)*conf.sat_restart_mult) {
      debug_print("SAT restarting!");
      last_restart = sat_confl;
      go_back_to(sat_start_dec_level);
      decisions.top().reset();
      if (!propagate()) goto start1;
      stats.sat_rst++;
      num_rst++;
      continue;
    }
  }

  if (true) {
    state = RESOLVED;
    T cnt = 1;
    if constexpr (weighted) {
      all_vars_in_comp(comp_manager->get_super_comp(decisions.at(sat_start_dec_level)), it) {
        uint32_t v = *it;
        if (v >= opt_indep_support_end) continue;
        sat_solution[v] = val(v);
      }
    }
    go_back_to(sat_start_dec_level);
    bool ret = propagate();
    assert(ret);
    assert(dec_level() == sat_start_dec_level);

    //  We need to multiply here, because some things may get re-propagated, and that will
    //  be unset, which would affect the weight calculated. Yes, chrono-bt is hard.
    if constexpr (weighted) {
      all_vars_in_comp(comp_manager->get_super_comp(decisions.at(sat_start_dec_level)), it) {
        uint32_t v = *it;
        if (v >= opt_indep_support_end) continue;
        debug_print(COLYEL "SAT solver -- mult var: " << setw(4) << v << " val: " << setw(3) << sat_solution[v]
          << " weight: " << setw(9) << get_weight(Lit(v, sat_solution[v] == T_TRI)) << COLDEF
          << " dec_lev: " << setw(5) << var(v).decision_level << " sat_start_dec_level: " << sat_start_dec_level);
        if (var(v).decision_level != INVALID_DL && var(v).decision_level <= sat_start_dec_level) continue;
        cnt *= get_weight(Lit(v, sat_solution[v] == T_TRI));
      }
      debug_print(COLYEL "SAT cnt will be: " << cnt);
    }
    decisions.top().var = 0;
    var(0).sublevel = old_sublev; // hack not to re-propagate everything.
    decisions.top().reset();
    decisions.top().change_to_right_branch();
    decisions.top().include_solution(cnt);
    if constexpr (!weighted) assert(decisions.top().total_model_count() == 1);
  }

end:
  assert(state != GO_AGAIN);
  order_heap.clear();
  sat_start_dec_level = -1;
  isindependent = true;
  debug_print("Exiting SAT mode. Declev: " << dec_level() << " sat: " << (int)sat
      << " trail below.");
  VERBOSE_DEBUG_DO(print_trail());
  if (sat) stats.sat_found_sat++;
  else stats.sat_found_unsat++;
  stats.sat_conflicts += stats.conflicts-conflicts_before;

  return sat;
}

template<typename T>
void Counter<T>::check_sat_solution() const {
  assert(sat_mode());
  bool good = true;

  for(const auto& off: long_irred_cls) {
    Clause& cl = *alloc->ptr(off);
    if (clause_falsified(cl)) {
      good = false;
      cout << "ERROR: SAT mode found a solution that falsifies a clause." << endl;
      print_cl(cl);
    }
  }

  for(const auto& off: long_red_cls) {
    Clause& cl = *alloc->ptr(off);
    if (clause_falsified(cl)) {
      good = false;
      cout << "ERROR: SAT mode found a solution that falsifies a clause." << endl;
      print_cl(cl);
    }
  }

  assert(good);
}

#ifdef BUDDY_ENABLED
#define mybdd_add(a,l) \
  do { \
  if (!(l).sign()) tmp |= bdd_ithvar(vmap_rev[(l).var()]); \
  else tmp |= bdd_nithvar(vmap_rev[(l).var()]); \
  } while(0)


template<typename T>
bool Counter<T>::should_do_buddy_count() const {
  auto d = decisions.at(dec_level()-1);
  const Comp& c = comp_manager->get_super_comp(d);
  if (c.nVars() >= 62 || c.nVars() <= 3 || c.num_long_cls() > conf.buddy_max_cls) {
    /* cout << "vars: " << c->nVars() << " cls: " << c->numBinCls() + c->num_long_cls() << endl; */
    return false;
  }
  return true;
}

template<typename T>
bool Counter<T>::do_buddy_count() {
  stats.buddy_called++;
  uint64_t cnt = buddy_count();

  debug_print("Buddy count: " << cnt);
  if (cnt > 0) {
    decisions.top().reset();
    decisions.top().change_to_right_branch();
    decisions.top().include_solution(cnt);
  } else {
    decisions.top().branch_found_unsat();
    decisions.top().change_to_right_branch();
    decisions.top().branch_found_unsat();
  }
  return cnt > 0;
}

// TODO Yash's ideas:
// * merge the BDDs in a tree-like manner
// * Mate: tune bdd_setcacheratio
// * need to use double bdd_satcountlnset(BDD r, BDD varset) to do projected counting
//   --> NOTE: double needs to be changed to int64_t
template<typename T>
uint64_t Counter<T>::buddy_count() {
  const Lit top_lit = trail.back();
  const uint32_t top_var = top_lit.var();
  const auto& s = decisions.top();
  auto const& sup_at = s.super_comp(); //TODO bad -- it doesn't take into account
                                       //that it could have already fallen into pieces
                                       //at current level
  const auto& c = comp_manager->at(sup_at);
  vmap.clear();
  vmap_rev.resize(nVars()+1);

  // variable mapping
  uint32_t proj_end = 63;
  bool proj = false;
  for(uint32_t i = 0; i < c->nVars(); i++) {
    const uint32_t var = c->vars_begin()[i];
    assert(var == top_var || val(var) == X_TRI);
    if (var >= indep_support_end && !proj) {
      proj_end = i;
      proj = true;
    }
    vmap.push_back(var);
    vmap_rev[vmap[i]] = i;
  }
  VERBOSE_DEBUG_DO(cout << "Vars in BDD: "; for(const auto& v: vmap) cout << v << " "; cout << endl);
  debug_print("proj_end: " << proj_end << " indep_support_end: " << indep_support_end);

  // The final built bdd
  auto bdd = bdd_true();

  // Long clauses
  uint32_t actual_long = 0;
  const auto& ana = comp_manager->get_ana();
  for (auto it_cl = c->cls_begin(); *it_cl != sentinel; it_cl++) {
    auto idx = *it_cl;
    debug_print("IDX: " << idx);
    Lit const* cl = ana.get_idx_to_cl(idx);
    VERBOSE_DEBUG_DO(cout << "Long cl." << endl;
      for(Lit const* l = cl; *l != SENTINEL_LIT; l++) cout << *l << " ";
      cout << endl);

    auto tmp = bdd_false();
    for(Lit const* l = cl; *l != SENTINEL_LIT; l++) {
      assert(l->var() == top_var || val(*l) != T_TRI);
      if (l->var() != top_var && val(*l) != X_TRI) continue;
      mybdd_add(tmp, *l);
    }
    bdd &= tmp;
    actual_long++;
  }

  // Binary clauses
  uint32_t actual_bin = 0;
  for(const auto& v: vmap) for(uint32_t i = 0; i < 2; i++) {
    Lit l(v, i);
    if (v != top_var && val(l) != X_TRI) continue;
    for(const auto& ws: watches[l].binaries) {
      if (!ws.irred() || ws.lit() < l) continue;
      if (ws.lit().var() != top_var && val(ws.lit()) == T_TRI) continue;
      SLOW_DEBUG_DO(assert(
            ws.lit().var() == top_var ||
            val(ws.lit()) == X_TRI)); // otherwise would have propagated/conflicted

      auto tmp = bdd_false();
      mybdd_add(tmp, l);
      auto l2 = ws.lit();
      mybdd_add(tmp, l2);
      bdd &= tmp;
      actual_bin++;

      debug_print("bin cl: " << l << " " << l2 << " 0");
    }
  }
  stats.buddy_num_bin_cls += actual_bin;
  stats.buddy_num_long_cls += actual_long;
  stats.buddy_num_vars += vmap.size();
  stats.buddy_max_bin_cls = std::max<uint64_t>(actual_bin, stats.buddy_max_bin_cls);
  stats.buddy_max_long_cls = std::max<uint64_t>(actual_long, stats.buddy_max_long_cls);
  stats.buddy_max_num_vars = std::max<uint64_t>(vmap.size(), stats.buddy_max_num_vars);

  VERBOSE_DEBUG_DO(
  if (actual_long != c->num_long_cls()) {
    cout << "WARN: numlong: " << c->num_long_cls() << " actual long: " << actual_long << endl;
  });
  /* cout << "bin cls: " << actual_bin << " long cls: " << actual_long << " vars: " << vmap.size() << endl; */

  assert(c->num_long_cls() == actual_long);

#ifdef VERBOSE_DEBUG
  std::stringstream fname;
  fname << "bdd-" << stats.buddy_called << ".dot";
  bdd_fnprintdot(fname.str().c_str(), bdd, proj_end);
  debug_print("BDD written to: " << fname.str());
#endif

  // Trick: when all is projected, we can set 64, and then subtract. This is
  // more generic. Otherwise, bdd cache will not match if we set the correct value
  // This way, we are more generic on non-projected, and can still use it on projected
  uint64_t cnt;
  if (proj_end == 63)
    cnt = bdd_satcount_i64(bdd, proj_end)>>(63-vmap.size());
  else
    cnt = bdd_satcount_i64(bdd, proj_end);
  VERBOSE_DEBUG_DO(
  cout << "cnt: " << cnt << endl;
  cout << "num bin cls: " << actual_bin << endl;
  cout << "num long cls: " << actual_long << endl;
  cout << "----------------------------------------------" << endl);

  return cnt;
}
#else
template<typename T>
bool Counter<T>::should_do_buddy_count() const {
  cout << "ERROR: you must recompile with buddy enabled for BDD counting to work" << endl;
  exit(-1);
  return false;
}
template<typename T>
bool Counter<T>::do_buddy_count() {
  cout << "ERROR: you must recompile with buddy enabled for BDD counting to work" << endl;
  exit(-1);
}
#endif

template<typename T>
template<typename T2>
void Counter<T>::check_cl_propagated_conflicted(T2& cl, uint32_t off) const {
  Lit unk = NOT_A_LIT;
  uint32_t num_unknown = 0;
  bool satisfied = false;
  for(const auto& l: cl) {
    if (is_true(l)) {satisfied = true; break;}
    if (is_unknown(l)) {num_unknown++; unk = l;}
    if (num_unknown > 1) break;
  }

  if (!satisfied && num_unknown == 1) {
    cout << "ERROR! Clause offs: " << off << " should have propagated: " << unk << endl;
    print_cl(cl);
    assert(false);
  }
  if (!satisfied && num_unknown == 0) {
    cout << "ERROR! Clause offs: " << off << " should have conflicted" << endl;
    print_cl(cl);
    assert(false);
  }
}

template<typename T>
void Counter<T>::check_all_propagated_conflicted() const {
  // Everything that should have propagated, propagated
  for(const auto& t: unit_cls) {
    if (val(t) != T_TRI) {
      cout << "Unit clause: " << t << " is set/falsified on trail." << endl;
      assert(false);
    }
  }
  for(const auto& off: long_irred_cls) {
    const Clause& cl = *alloc->ptr(off);
    check_cl_propagated_conflicted(cl, off);
  }
  for(const auto& off: long_red_cls) {
    const Clause& cl = *alloc->ptr(off);
    check_cl_propagated_conflicted(cl, off);
  }

  all_lits(i) {
    Lit lit(i/2, i%2);
    if (val(lit) == T_TRI) continue;
    for(const auto& ws: watches[lit].binaries) {
      if (val(ws.lit()) == T_TRI) continue;
      if (val(lit) == F_TRI) {
        if (val(ws.lit()) == X_TRI) {
          cout << "Should have propagated lit: " << ws.lit() << " due to binary clause " << lit << " " << ws.lit() << endl;
          assert(false);
        }
        if (val(ws.lit()) == F_TRI) {
          cout << "Falsified binary clause: " << lit << " " << ws.lit() << endl;
          assert(false);
        }
      }
      if (val(ws.lit()) == F_TRI) {
        if (val(lit) == X_TRI) {
          cout << "Should have propagated lit: " << lit << " due to binary clause " << lit << " " << ws.lit() << endl;
          assert(false);
        }
        if (val(lit) == F_TRI) {
          cout << "Falsified binary clause: " << lit << " " << ws.lit() << endl;
          assert(false);
        }
      }
    }
  }

#ifdef SLOW_DEBUG
  for(const auto& cl: debug_irred_cls) check_cl_propagated_conflicted(cl);
#endif
}

template<typename T>
void Counter<T>::v_backup() {
  for(const auto& off: long_irred_cls) {
    const Clause& cl = *alloc->ptr(off);
    vector<Lit> lits(cl.begin(), cl.end());
    v_backup_cls.push_back(lits);
  }
  for(const auto& off: long_red_cls) {
    const Clause& cl = *alloc->ptr(off);
    vector<Lit> lits(cl.begin(), cl.end());
    v_backup_cls.push_back(lits);
  }
  for(const auto& ws: watches) {
    vector<ClOffsBlckL> tmp(ws.watch_list_.begin(), ws.watch_list_.end());
    v_backup_watches.push_back(tmp);
  }
}

template<typename T>
void Counter<T>::v_restore() {
  uint32_t at = 0;
  for(const auto& off: long_irred_cls) {
    Clause& cl = *alloc->ptr(off);
    auto& lits = v_backup_cls[at];
    for(uint32_t i = 0; i < cl.size(); i++) {
      cl[i]=lits[i];
    }
    at++;
  }
  for(const auto& off: long_red_cls) {
    Clause& cl = *alloc->ptr(off);
    auto& lits = v_backup_cls[at];
    for(uint32_t i = 0; i < cl.size(); i++) {
      cl[i]=lits[i];
    }
    at++;
  }

  at = 0;
  for(auto& ws: watches) {
    ws.watch_list_.clear();
    ws.watch_list_ = v_backup_watches[at];
    at++;
  }
  v_backup_watches.clear();
  v_backup_watches.shrink_to_fit();
  v_backup_cls.clear();
  v_backup_cls.shrink_to_fit();
}

template<typename T>
void Counter<T>::set_lit(const Lit lit, int32_t dec_lev, Antecedent ant) {
  assert(val(lit) == X_TRI);
  if (ant.isNull())
    debug_print("set_lit called with a decision. Lit: " << lit << " lev: " << dec_lev << " cur dec lev: " << dec_level());
  else debug_print("-> lit propagated: " << lit << " trail pos will be: " << trail.size() << " cur dec lev: " << dec_level());

  debug_print("setting lit: " << lit << " to lev: " << dec_lev << " cur val: " << lit_val_str(lit) << " ante: " << ant << " sublev: " << trail.size());
  var(lit).decision_level = dec_lev;
  var(lit).ante = ant;
  if (!ant.isNull()) {
    var(lit).last_polarity = lit.sign();
  }
  var(lit).sublevel = trail.size();
  qhead = std::min<uint32_t>(qhead, trail.size());
  trail.push_back(lit);
  __builtin_prefetch(watches[lit.neg()].binaries.data());
  __builtin_prefetch(watches[lit.neg()].watch_list_.data());
  if constexpr (weighted) if (dec_lev <= dec_level() && get_weight(lit) != 1) {
    int32_t until = decisions.size();
    if (sat_mode()) until = std::min((int)decisions.size(), sat_start_dec_level);
    for(int32_t i = dec_lev; i < until; i++) {
      debug_print("set_lit, compensating weight. i: " << i << " dec_lev: " << dec_lev);
      bool found = false;
      if (vars_act_dec.size() <= i*(nVars()+1)) break;
      uint64_t* at = vars_act_dec.data()+i*(nVars()+1);
      found = (at[0] == at[lit.var()]);
      /* debug_print("dec val compare: " << at[0]); */
      // Not found in parent, so not in any children for sure
      if (!found) {
        debug_print("Var not found in parent, so not in children for sure. Exiting");
        break;
      } else debug_print("Var found in parent.");
      if (i > dec_lev && found) decisions[i].include_solution_left_side(1/get_weight(lit));

      bool found_in_children = false;
      const auto& s = decisions.at(i);
      debug_print("s.get_unprocessed_comps_end(): " << s.get_unprocessed_comps_end()
          << " s.remaining_comps_ofs(): " << s.remaining_comps_ofs()
          << " comp_manager->size: " << comp_manager->get_comp_stack().size());
      for(int comp_at = s.get_unprocessed_comps_end()-1; comp_at >= (int)s.remaining_comps_ofs() &&
          comp_at < (int)comp_manager->get_comp_stack().size(); comp_at--) {
        const auto& c2 = comp_manager->at(comp_at);
        VERBOSE_DEBUG_DO(cout << "vars in side comp: ";
          all_vars_in_comp(*c2, v) VERBOSE_DEBUG_DO(cout << *v << " ");
          cout << endl;);
        all_vars_in_comp(*c2, v) if (*v == lit.var()) {found_in_children = true;break;}
      }
      debug_print("found in children: " << found_in_children);
      if (!found_in_children) {
        // Not found in children, so it must have been already processed and multiplied in. Compensate.
        assert((int)decisions.size() > i);
        if (decisions[i].get_branch_sols() != 0) decisions[i].include_solution(1/get_weight(lit));
      }
    }
  }
  values[lit] = T_TRI;
  values[lit.neg()] = F_TRI;
}

template<typename T>
template<class T2>
bool Counter<T>::clause_falsified(const T2& cl) const {
  for(const auto&l: cl) {
    if (val(l) != F_TRI) return false;
  }
  return true;
}

template<typename T>
void Counter<T>::end_irred_cls() {
  seen.clear();
  seen.resize(2*(nVars()+2), 0);
  stats.max_cache_size_bytes = conf.maximum_cache_size_MB*1024*1024;

  delete comp_manager;
  comp_manager = new CompManager(conf, stats, values, indep_support_end, this);
  comp_manager->getrandomseedforclhash();

  init_decision_stack();
  simple_preprocess();
  ended_irred_cls = true;

  // This below will initialize the disjoint component analyzer (ana)
  comp_manager->initialize(watches, alloc, long_irred_cls);
}

#ifdef BUDDY_ENABLED
void my_gbchandler(int pre, bddGbcStat *) {
   if (!pre) {
      /* printf("Garbage collection #%d: %d nodes / %d free", s->num, s->nodes, s->freenodes); */
      /* printf(" / %.1fs / %.1fs total\n", */
      /* (double)s->time/(double)(CLOCKS_PER_SEC), */
      /* (double)s->sumtime/(double)CLOCKS_PER_SEC); */
   }
}
#endif

template<typename T>
Counter<T>::Counter(const CounterConfiguration& _conf) :
    conf(_conf)
    , stats(_conf)
    , mtrand(_conf.seed)
    , order_heap(VarOrderLt(Counter<T>::watches)) {
  sat_solver = new CMSat::SATSolver;
  sat_solver->set_prefix("c o ");
  alloc = new ClauseAllocator<T>(_conf);
  lbd_cutoff = conf.base_lbd_cutoff;
  BUDDY_DO(if (conf.do_buddy) {
    bdd_init(10000, 100000);
    bdd_gbc_hook(my_gbchandler);
    bdd_setvarnum(63);
    bdd_autoreorder(BDD_REORDER_NONE);
  });
}

template<typename T>
Counter<T>::~Counter() {
  delete comp_manager;
  BUDDY_DO(if (conf.do_buddy) bdd_done());
  delete alloc;
  delete sat_solver;
}

template<typename T>
void Counter<T>::simple_preprocess() {
  verb_print(2, "[simple-preproc] Running.");
  for (const auto& lit : unit_cls) {
    assert(!exists_unit_cl_of(lit.neg()) && "Formula is not UNSAT, we ran CMS before");
    if (val(lit) == X_TRI) {
      set_lit(lit, 0);
      verb_print(2, "[simple-preproc] set: " << lit);
    }
    assert(val(lit) == T_TRI);
  }

  verb_print(2, "[simple-preproc] propagating.");
  bool succeeded = propagate();
  release_assert(succeeded && "We ran CMS before, so it cannot be UNSAT");
  for(const auto& t: trail) if (!exists_unit_cl_of(t)) unit_cls.push_back(t);
  init_decision_stack();
  qhead = 0;

  // Remove for reasons for 0-level clauses, these may interfere with
  // deletion of clauses during subsumption
  for(auto& v: var_data) {v.ante = Antecedent();}
  verb_print(2, "[simple-preproc] finished.");
}

// TODO Yash we should do Jeroslow-Wang heuristic, i.e. 1/2 for binary, 1/3 for tertiary, etc.
template<typename T>
void Counter<T>::init_activity_scores() {
  if (!conf.do_init_activity_scores) return;
  all_lits(x) {
    Lit l(x/2, x%2);
    for (const auto& ws: watches[l].binaries) {
      if (!ws.red()) watches[l].activity++;
    }
  }
  for(const auto& off: long_irred_cls) {
    const auto& cl = *alloc->ptr(off);
    for(const auto& l: cl) watches[l].activity++;
  }
}

template<typename T>
void Counter<T>::check_all_cl_in_watchlists() const {
  auto red_cls2 = long_red_cls;
  // check for duplicates
  std::sort(red_cls2.begin(), red_cls2.end());
  for(uint32_t i = 1; i < red_cls2.size(); i++) {
    assert(red_cls2[i-1] != red_cls2[i]);
  }

  for(const auto& offs: long_red_cls) {
    const auto& cl = *alloc->ptr(offs);
    if (!find_offs_in_watch(watches[cl[0]].watch_list_, offs)) {
      cout << "ERROR: Did not find watch cl[0]!!" << endl;
      assert(false);
      exit(-1);
    }
    if (!find_offs_in_watch(watches[cl[1]].watch_list_, offs)) {
      cout << "ERROR: Did not find watch cl[1]!!" << endl;
      assert(false);
      exit(-1);
    }
  }
}

template<typename T>
bool Counter<T>::find_offs_in_watch(const vector<ClOffsBlckL>& ws, ClauseOfs off)  const
{
  for (auto& w: ws) if (w.ofs == off) { return true; }
  return false;
}

template<typename T>
struct ClSorter {
  ClSorter(ClauseAllocator<T>* _alloc, uint32_t _lbd_cutoff) :
    alloc(_alloc), lbd_cutoff(_lbd_cutoff) {}

  bool operator()(ClauseOfs& a, ClauseOfs& b) const {
    const auto& ah = *alloc->ptr(a);
    const auto& bh = *alloc->ptr(b);
    assert(ah.red);
    assert(bh.red);
    if (ah.lbd <= lbd_cutoff || bh.lbd <= lbd_cutoff) return ah.lbd < bh.lbd;
    if (ah.used != bh.used) return ah.used > bh.used;
    return ah.total_used > bh.total_used;
  }
  ClauseAllocator<T>* alloc;
  const uint32_t lbd_cutoff;
};

template<typename T>
void Counter<T>::reduce_db() {
  stats.reduce_db++;
  const auto cls_before = long_red_cls.size();

  vector<ClauseOfs> tmp_red_cls = long_red_cls;
  long_red_cls.clear();
  num_low_lbd_cls = 0;
  num_used_cls = 0;
  uint32_t cannot_be_del = 0;
  sort(tmp_red_cls.begin(), tmp_red_cls.end(), ClSorter(alloc, lbd_cutoff));
  int64_t new_decs = stats.decisions - last_reducedb_dec;
  int64_t new_confls = stats.conflicts - last_reducedb_confl;
  uint32_t target = conf.rdb_cls_target;
  if (new_confls*4 > new_decs) target *= 2;
  else if (new_confls*8 > new_decs) target *= 1.5;
  else if (new_confls*16 > new_decs) target *= 1;
  else if (new_confls*32 > new_decs) target *= 0.8;
  else target *= 0.4;

  for(uint32_t i = 0; i < tmp_red_cls.size(); i++){
    const ClauseOfs& off = tmp_red_cls[i];
    auto& h = *alloc->ptr(off);
    if (h.lbd <= lbd_cutoff) num_low_lbd_cls++;
    else if (h.used) num_used_cls++;

    bool can_be_del = red_cl_can_be_deleted(off);
    cannot_be_del += !can_be_del;
    if (can_be_del && h.lbd > lbd_cutoff
        && (!conf.rdb_keep_used || !h.used)
        && i > target + num_low_lbd_cls + (conf.rdb_keep_used ? num_used_cls : 0)) {
      delete_cl(off);
      stats.cls_deleted_since_compaction++;
      stats.cls_removed++;
    } else {
      long_red_cls.push_back(off);
      h.used = 0;
    }
  }

  // Update LBD cutoff
  if (stats.conflicts > (100ULL*1000ULL) && lbd_cutoff == conf.base_lbd_cutoff
      && num_low_lbd_cls < 50 && conf.update_lbd_cutoff) {
    verb_print(1, " [rdb] bumping rdb cutoff to 3");
    lbd_cutoff++;
  }

  if (conf.verb >= 2 || stats.reduce_db % 3 == 1) {
    verb_print(1, "[rdb] cls before: " << cls_before << " after: " << long_red_cls.size()
      << " low lbd: " << num_low_lbd_cls
      << " lbd cutoff: " << lbd_cutoff
      << " target computed: " << target
      << " cannot be del : " << cannot_be_del
      << " used: " << num_used_cls << " rdb: " << stats.reduce_db);
    verb_print(2, "Time until now: " << cpu_time());}
}

template<typename T>
bool Counter<T>::red_cl_can_be_deleted(ClauseOfs off){
  // only first literal may possibly have cl_ofs as antecedent
  Clause& cl = *alloc->ptr(off);
  if (is_antec_of(off, cl[0])) return false;
  return true;
}

template<typename T>
void Counter<T>::delete_cl(const ClauseOfs off){
  Clause& cl = *alloc->ptr(off);
  watches[cl[0]].del_c(off);
  watches[cl[1]].del_c(off);
  alloc->clause_free(off);
}

template<typename T>
void Counter<T>::new_vars(const uint32_t n) {
  if (num_vars_set) {
    cout << "ERROR: you can only call new_vars() once!" << endl;
    exit(-1);
  }
  sat_solver->new_vars(n);

  assert(var_data.empty());
  assert(values.empty());
  assert(watches.empty());
  assert(unit_cls.empty());
  assert(long_red_cls.empty());
  assert(weights.empty());
  assert(sat_solution.empty());

  var_data.resize(n + 1);
  values.resize(n + 1, X_TRI);
  watches.resize(n + 1);
  lbd_helper.resize(n+1, 0);
  if constexpr (weighted) sat_solution.resize(n+1);
  if constexpr (weighted) weights.resize(2*(n + 1), 1);
  num_vars_set = true;
}

template<typename T>
Clause* Counter<T>::add_cl(const vector<Lit> &lits, bool red) {
  if (lits.size() == 1) {
    assert(!exists_unit_cl_of(lits[0].neg()) && "UNSAT is not dealt with");
    if (!exists_unit_cl_of(lits[0])) unit_cls.push_back(lits[0]);
    return nullptr;
  }

  if (lits.size() == 2) {
    add_bin_cl(lits[0], lits[1], red);
    return nullptr;
  }

  Clause* cl = alloc->new_cl(red, lits.size());
  for(uint32_t i = 0; i < lits.size(); i ++) (*cl)[i] = lits[i];
  attach_cl(alloc->get_offset(cl), lits);
  return cl;
}

template<typename T>
bool Counter<T>::add_irred_cl(const vector<Lit>& lits_orig) {
  if (!ok) return ok;
  if (!sat_solver->add_clause(ganak_to_cms_cl(lits_orig))) { ok = false; return ok; }

  vector<Lit> lits;
  for(const auto& l: lits_orig) {
    if (val(l) == T_TRI) return ok;
    if (val(l) == X_TRI) lits.push_back(l);
  }
  if (lits.empty()) {
    cout << "ERROR: UNSAT should have been caught by external SAT solver" << endl;
    exit(-1);
  }
  for(const auto& l: lits) assert(l.var() <= nVars());
  if (!remove_duplicates(lits)) return ok;

  stats.incorporateIrredClauseData(lits);
  Clause* cl = add_cl(lits, false);
  auto off = alloc->get_offset(cl);
  if (cl) long_irred_cls.push_back(off);
  SLOW_DEBUG_DO(debug_irred_cls.push_back(lits));
  return ok;
}

template<typename T>
bool Counter<T>::add_red_cl(const vector<Lit>& lits_orig, int lbd) {
  if (!sat_solver->add_clause(ganak_to_cms_cl(lits_orig))) { ok = false; return ok; }

  vector<Lit> lits;
  for(const auto& l: lits_orig) {
    if (val(l) == T_TRI) return ok;
    if (val(l) == X_TRI) lits.push_back(l);
  }
  if (lits.empty()) {
    cout << "ERROR: UNSAT should have been caught by external SAT solver" << endl;
    exit(-1);
  }
  for(const auto& l: lits) assert(l.var() <= nVars());
  if (!remove_duplicates(lits)) return ok;
  Clause* cl = add_cl(lits, true);
  if (cl) {
    auto off = alloc->get_offset(cl);
    long_red_cls.push_back(off);
    if (lbd == -1) lbd = lits.size();
    cl->lbd = lbd;
    assert(cl->red);
  }
  return ok;
}

template<typename T>
void Counter<T>::reactivate_comps_and_backtrack_trail([[maybe_unused]] bool check_ws) {
  debug_print("->reactivate and backtrack. Dec lev: " << dec_level() << " top declevel sublev: " << var(decisions.top().var).sublevel <<  "...");
  auto jt = top_declevel_trail_begin();
  auto it = jt;
  int32_t off_by = 0;
  for (; it != trail.end(); it++) {
    debug_print("Backing up, checking lit: " << std::setw(5) << *it << " at: " << it-trail.begin());
    SLOW_DEBUG_DO(assert(it->var() != 0));
    const int32_t dl = var(*it).decision_level;
    assert(dl != -1);
    if (dl < dec_level()) {
      off_by++;
      var(*it).sublevel = jt - trail.begin();
      *jt++ = *it;
      debug_print("Backing up, setting: " << std::setw(5) << *it << " at lev: " << std::setw(4) << dl
          << " to sublev: " << var(*it).sublevel);
    } else {
      debug_print("Backing up, unsetting: " << std::right << std::setw(8) << *it
          << " lev: " << std::setw(4) << var(*it).decision_level
          << " ante was: " << var(*it).ante);
      if (sat_mode() && !order_heap.in_heap(it->var())) order_heap.insert(it->var());
      unset_lit(*it);
    }
  }
  VERY_SLOW_DEBUG_DO(if (check_ws && !check_watchlists()) {
      print_trail(false, false);assert(false);});
  if (!sat_mode()) comp_manager->clean_remain_comps_of(decisions.top());
  trail.resize(jt - trail.begin());
  if (dec_level() == 0) qhead = 0;
  else qhead = std::min<int32_t>(trail.size()-off_by, qhead);
  if (!sat_mode()) {
    decisions.top().reset_remain_comps();
  }
}

template<typename T>
void Counter<T>::set_lit_weight(Lit l, const T& w) {
  if (l.var() >= opt_indep_support_end) {
    cerr << "ERROR: Trying to set weight of a variable that is not in the (optional) independent support."
      " Var: " << l << " opt_indep_support_end: " << opt_indep_support_end << endl;
    exit(-1);
  }
  verb_print(2, "Setting weight of " << l << " to " << w);
  weights[l.raw()] = w;
  if (w == 0) add_irred_cl({l.neg()});
}

template<typename T>
void Counter<T>::init_decision_stack() {
    decisions.clear();
    trail.clear();
    // initialize the stack to contain at least level zero
    decisions.push_back(StackLevel<T>(
          1, // super comp
          2)); //comp stack offset

    // I guess this is needed so the system later knows it's fully counted
    // since this is only a dummy.
    decisions.back().change_to_right_branch();
  }

template<typename T>
string Counter<T>::lit_val_str(Lit lit) const {
    if (values[lit] == F_TRI) return "FALSE";
    else if (values[lit] == T_TRI) return "TRUE";
    else return "UNKN";
  }

template<typename T>
string Counter<T>::val_to_str(const TriValue& tri) const {
    if (tri == F_TRI) return "FALSE";
    else if (tri == T_TRI) return "TRUE";
    else return "UNKN";
  }

template class Counter<mpz_class>;
template class Counter<mpfr::mpreal>;
template class Counter<mpq_class>;
