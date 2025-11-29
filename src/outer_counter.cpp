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

#include "outer_counter.hpp"
#include <iostream>
#include <algorithm>
#include <future>
#include <thread>
#include "TreeDecomposition.hpp"
#include "IFlowCutter.hpp"
#include "arjun/arjun.h"
#include "counter.hpp"
#include "time_mem.hpp"
#include "common.hpp"

namespace GanakInt {

CMSat::Lit ganak_to_cms_lit(const GanakInt::Lit& l) {
  return CMSat::Lit(l.var()-1, !l.sign());
}

GanakInt::Lit cms_to_ganak_lit(const CMSat::Lit& l) {
  return GanakInt::Lit(l.var()+1, !l.sign());
}

template<typename T>
vector<uint32_t> ganak_to_cms_vars(const T& vars) {
  vector<uint32_t> cms_vars; cms_vars.reserve(vars.size());
  for(const auto& v: vars) cms_vars.push_back(v-1);
  return cms_vars;
}

inline vector<GanakInt::Lit> cms_to_ganak_cl(const vector<CMSat::Lit>& cl) {
  vector<GanakInt::Lit> ganak_cl; ganak_cl.reserve(cl.size());
  for(const auto& l: cl) ganak_cl.push_back(cms_to_ganak_lit(l));
  return ganak_cl;
}

inline vector<CMSat::Lit> ganak_to_cms_cl(const vector<GanakInt::Lit>& cl) {
  vector<CMSat::Lit> cms_cl; cms_cl.reserve(cl.size());
  for(const auto& l: cl) cms_cl.push_back(ganak_to_cms_lit(l));
  return cms_cl;
}

FF OuterCounter::count(uint8_t bits_jobs, int num_threads) {
  verb_print(2, "[par] bits_jobs: " << (uint32_t)bits_jobs);
  verb_print(2, "[par] num_threads: " << (uint32_t)num_threads);
  verb_print(2, "[par] TD: " << (conf.do_td ? "enabled" : "disabled"));
  verb_print(2, "[par] Number of variables: " << nvars);
  verb_print(2, "[par] Independent support size: " << indep_support.size());
  verb_print(2, "[par] TD variable limit: " << conf.td_varlim);
  verb_print(2, "[par] AppMC timeout: " << conf.appmc_timeout);
  if ((num_threads > 1 || num_threads == -1) &&
      bits_jobs > 0 && conf.do_td && nvars > 30 && indep_support.size() > 10 &&  nvars <= conf.td_varlim && conf.appmc_timeout < 0) {
    verb_print(1, "[par] Attempting parallel counting with bits_jobs: " << bits_jobs
        << ", num_threads: " << (uint32_t)num_threads);
    return count_with_parallel(bits_jobs, num_threads);
  } else {
    verb_print(1, "[par] Using non-parallel counting");
    return count_regular();
  }
}

FF OuterCounter::count_regular() {
  auto counter = std::make_unique<Counter>(conf, fg);
  counter->new_vars(nvars);
  counter->set_indep_support(indep_support);
  counter->set_optional_indep_support(opt_indep_support);
  for (const auto& [lit, weight] : lit_weights)
    counter->set_lit_weight(lit, weight->dup());
  counter->set_generators(generators);

  for(const auto& cl : irred_cls) counter->add_irred_cl(cl);
  counter->end_irred_cls();
  for(const auto& p : red_cls) counter->add_red_cl(p.first, (int)p.second);
  auto ret = counter->outer_count();

  assert(num_cache_lookups == 0);
  assert(max_cache_elems == 0);
  num_cache_lookups = counter->get_stats().num_cache_look_ups;
  max_cache_elems = counter->get_cache()->get_max_num_entries();
  count_is_approximate |= counter->get_is_approximate();
  return ret;
}

void run_arjun(ArjunNS::SimplifiedCNF& cnf) {
  /* double my_time = cpu_time(); */
  ArjunNS::Arjun arjun;
  arjun.set_verb(0);
  /* arjun.set_or_gate_based(arjun_gates); */
  /* arjun.set_xor_gates_based(arjun_gates); */
  /* arjun.set_ite_gate_based(arjun_gates); */
  /* arjun.set_irreg_gate_based(arjun_gates); */
  /* arjun.set_extend_max_confl(arjun_extend_max_confl); */
  /* arjun.set_probe_based(do_probe_based); */
  arjun.set_simp(0);
  arjun.set_backw_max_confl(100);
  arjun.set_oracle_find_bins(0);
  arjun.set_cms_glob_mult(0.0001);

  ArjunNS::SimpConf simp_conf;
  /* simp_conf.iter1 = 1; */
  simp_conf.iter1 = 0;
  ArjunNS::Arjun::ElimToFileConf etof_conf;
  arjun.standalone_minimize_indep(cnf, false);
  if (cnf.get_sampl_vars().size() >= 10) {
    arjun.standalone_elim_to_file(cnf, etof_conf, simp_conf);
  } else cnf.renumber_sampling_vars_for_ganak();
  /* verb_print(1, "Arjun T: " << (cpu_time()-my_time)); */
}

void setup_ganak(const ArjunNS::SimplifiedCNF& cnf, OuterCounter& counter) {
  cnf.check_sanity();
  counter.new_vars(cnf.nVars());

  set<uint32_t> tmp;
  for(auto const& s: cnf.sampl_vars) tmp.insert(s+1);
  counter.set_indep_support(tmp);
  if (cnf.get_opt_sampl_vars_set()) {
    tmp.clear();
    for(auto const& s: cnf.opt_sampl_vars) tmp.insert(s+1);
  }
  counter.set_optional_indep_support(tmp);

  if (cnf.weighted) {
    for(const auto& t: cnf.weights) {
      counter.set_lit_weight(Lit(t.first+1, true), t.second.pos);
      counter.set_lit_weight(Lit(t.first+1, false), t.second.neg);
    }
  }

  for(const auto& cl: cnf.clauses) counter.add_irred_cl(cms_to_ganak_cl(cl));
  for(const auto& cl: cnf.red_clauses) counter.add_red_cl(cms_to_ganak_cl(cl));
}

FF OuterCounter::count_with_parallel(uint8_t bits_jobs, int num_threads) {
  assert(conf.appmc_timeout < 0 && "TD-parallel not compatible with AppMC");
  assert(bits_jobs > 0);
  assert(bits_jobs <= 20);
  double td_start_time = cpu_time();

  // Build primal graph from clauses
  TWD::Graph primal(nvars);
  for (const auto& cl : irred_cls) {
    for (size_t i = 0; i < cl.size(); i++) {
      for (size_t j = i + 1; j < cl.size(); j++) {
        uint32_t v1 = cl[i].var() - 1;
        uint32_t v2 = cl[j].var() - 1;
        if (v1 != v2) primal.addEdge(v1, v2);
      }
    }
  }

  auto nodes = opt_indep_support.size();
  for(uint32_t i = nodes; i < nvars; i++) {
    primal.contract(i, conf.td_max_edges*100);
    if (primal.numEdges() > conf.td_max_edges*100 ) break;
  }
  verb_print(1, "[par] nodes: " << nodes << " nvars: " << nvars << " edges: " << primal.numEdges());

  // Run FlowCutter to get tree decomposition
  TWD::IFlowCutter fc(primal.numNodes(), primal.numEdges(), conf.verb);
  fc.importGraph(primal);

  // Compute TD with reduced steps to avoid timeout
  auto tdec  = fc.constructTD(conf.td_steps / 3, conf.td_iters / 3);

  // Find centroid
  int centroid_id = tdec.centroid(primal.numNodes(), conf.verb);
  vector<int> centroid_bag = tdec.Bags()[centroid_id];

  verb_print(1, "[par] TD width: " << tdec.width()
          << ", centroid bag size: " << centroid_bag.size()
          << ", TD time: " << td_start_time-cpu_time() << "s");

  // If centroid bag is empty or very small, just use regular counting
  if (centroid_bag.size() < bits_jobs) {
    bits_jobs = (int)std::log2(centroid_bag.size());
    verb_print(2, "[par] Centroid bag smaller than 2**bits_jobs, using bits_jobs: " << bits_jobs);
  }
  vector<uint32_t> var_freq(nvars+1, 0);
  for(const auto& cl: irred_cls) {
    for(const auto& l: cl) var_freq[l.var()]++;
  }
  std::sort(centroid_bag.begin(), centroid_bag.end(),
            [&](uint32_t a, uint32_t b) {
              return var_freq[a+1] > var_freq[b+1];
            });
  for(uint32_t i = 0; i < centroid_bag.size(); i++) {
    verb_print(2, "[par] var " << centroid_bag[i]+1
                << " with frequency " << var_freq[centroid_bag[i]+1]);
  }

  uint64_t num_jobs = 1ULL << bits_jobs;
  assert(1ULL<<bits_jobs == num_jobs);
  if (num_threads == -1) {
    num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 2; // Fallback
    verb_print(1, "[par] Auto-detected num_threads: " << num_threads);
  }

  if (conf.verb >= 1) {
    std::cout << "c o [par] Using " << num_jobs << " jobs for "
      << num_threads << " concurrent threads"
      << " centroid variables: ";
    for (uint32_t i = 0; i < bits_jobs; i++) std::cout << (centroid_bag[i] + 1) << " ";
    std::cout << std::endl;
  }

  auto worker = [&](uint64_t num) -> FF {
    ArjunNS::SimplifiedCNF cnf(fg);
    cnf.new_vars(nvars);
    cnf.set_sampl_vars(ganak_to_cms_vars(indep_support));
    cnf.set_opt_sampl_vars(ganak_to_cms_vars(opt_indep_support));
    for (const auto& [lit, weight] : lit_weights)
      cnf.set_lit_weight(ganak_to_cms_lit(lit), weight->dup());
    for (const auto& cl : irred_cls) cnf.add_clause(ganak_to_cms_cl(cl));

    // Add unit clauses for centroid variable assignment
    for (size_t i = 0; i < bits_jobs; i++) {
      uint32_t var = centroid_bag[i];
      bool sign = (num >> i) & 1;
      CMSat::Lit unit_lit(var, sign);
      cnf.add_clause({unit_lit});
      verb_print(2, "[par] Thread " << num << " fixing var " << var
                  << " to " << (sign ? "true" : "false"));
    }
    for (const auto& [cl, lbd] : red_cls)
      cnf.add_red_clause(ganak_to_cms_cl(cl));
    if (true) run_arjun(cnf);
    if (cnf.multiplier_weight != fg->zero()) {
      auto local_conf = conf;
      local_conf.verb = 0; // disable verb for threads
      auto counter = std::make_unique<OuterCounter>(local_conf, fg);
      setup_ganak(cnf, *counter);
      auto ret = counter->count();
      num_cache_lookups += counter->get_num_cache_lookups();
      max_cache_elems = std::max(max_cache_elems, counter->get_max_cache_elems());
      count_is_approximate |= counter->get_is_approximate();
      *ret *= *cnf.multiplier_weight;
      return ret;
    } else {
      return fg->zero();
    }
  };

  // Thread pool: maintain at most max_concurrent active threads
  std::vector<std::pair<uint64_t, std::future<FF>>> active_futures;
  uint64_t next_task = 0;
  FF total_count = fg->zero();

  // Launch initial batch of threads
  assert(num_threads > 0);
  for (uint32_t i = 0; i < (uint32_t)num_threads && next_task < num_jobs; i++) {
    active_futures.push_back({next_task, std::async(std::launch::async, worker, next_task)});
    next_task++;
  }

  // Process results and launch new tasks as threads complete
  while (!active_futures.empty()) {
    for (size_t i = 0; i < active_futures.size(); ) {
      auto& [task_id, future] = active_futures[i];
      // Check if this future is ready
      if (future.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
        // Collect result
        FF partial = future.get();
        *total_count += *partial;
        verb_print(2, "[par] Task " << task_id << " completed");

        // Remove from active list
        active_futures.erase(active_futures.begin() + i);

        // Launch next task if available
        if (next_task < num_jobs) {
          active_futures.push_back({next_task, std::async(std::launch::async, worker, next_task)});
          verb_print(2, "[par] Launched task " << next_task);
          next_task++;
        }
      } else {
        i++;
      }
    }

    // Small sleep to avoid busy waiting
    if (!active_futures.empty()) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  }

  return total_count;
}

void OuterCounter::print_indep_distrib() const {
  cout << "c o indep/optional/none distribution: ";
  std::set<uint32_t> all_opt_indep = opt_indep_support;
  std::set<uint32_t> all_indep(indep_support);
  for(uint32_t i = 0; i <= nvars; i++) {
    if (all_opt_indep.count(i)) {
      if (all_indep.count(i)) cout << "I";
      else cout << "O";
    } else cout << "N";
  }
  cout << endl;
}

}
