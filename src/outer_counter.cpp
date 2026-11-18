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
#include <set>
#include <future>
#include <numeric>
#include <tuple>
#include <map>
#include <sstream>
#include <chrono>
#include <thread>
#include <treedecomp/TreeDecomposition.hpp>
#include <treedecomp/IFlowCutter.hpp>
#include "arjun/arjun.h"
#include "counter.hpp"
#include "time_mem.hpp"
#include "common.hpp"
#include "ganak.hpp"

using std::set;

namespace GanakInt {


GanakInt::Lit cms_to_ganak_lit(const CMSat::Lit& l) {
  return GanakInt::Lit(l.var()+1, !l.sign());
}

template<std::ranges::sized_range T>
vector<uint32_t> ganak_to_cms_vars(const T& vars) {
  vector<uint32_t> cms_vars;
  cms_vars.reserve(vars.size());
  std::ranges::transform(vars, std::back_inserter(cms_vars), [](uint32_t v) { return v - 1; });
  return cms_vars;
}

inline vector<GanakInt::Lit> cms_to_ganak_cl(const vector<CMSat::Lit>& cl) {
  vector<GanakInt::Lit> ganak_cl;
  ganak_cl.reserve(cl.size());
  std::ranges::transform(cl, std::back_inserter(ganak_cl), cms_to_ganak_lit);
  return ganak_cl;
}
// ganak_to_cms_cl for vector<Lit> is handled by the template in counter.hpp

FF OuterCounter::count(uint8_t bits_jobs, int num_threads, const bool debug_threads) {
  verb_print(2, "[par] bits_jobs: " << (uint32_t)bits_jobs);
  verb_print(2, "[par] num_threads: " << (uint32_t)num_threads);
  verb_print(2, "[par] TD: " << (conf.do_td ? "enabled" : "disabled"));
  verb_print(2, "[par] Number of variables: " << nvars);
  verb_print(2, "[par] Independent support size: " << indep_support.size());
  verb_print(2, "[par] TD variable limit: " << conf.td_varlim);
  verb_print(2, "[par] AppMC timeout: " << conf.appmc_timeout);
  if ((debug_threads || num_threads > 1 || num_threads == -1) &&
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
  for(const auto& p : red_cls) counter->add_red_cl(p.first, (int)p.second);
  auto ret = counter->outer_count();

  assert(num_cache_lookups == 0);
  assert(max_cache_elems == 0);
  num_cache_lookups = counter->get_stats().num_cache_look_ups;
  // if CNF is UNSAT, we early exit, and cache is not initialized, so we need this check
  if (counter->has_comp_manager()) {
    max_cache_elems = counter->get_cache()->get_max_num_entries();
  }
  count_is_approximate |= counter->get_is_approximate();
  return ret;
}

void run_arjun(ArjunNS::SimplifiedCNF& cnf, uint32_t verb) {
  /* double my_time = cpu_time(); */
  ArjunNS::Arjun arjun;
  arjun.set_verb(verb >= 3 ? verb : 0);
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
  ArjunNS::Arjun::ElimToFileConf const etof_conf;
  ArjunNS::Arjun::InterpConf iconf;
  arjun.standalone_minimize_indep(cnf, iconf, false);
  if (cnf.get_sampl_vars().size() >= td_at_or_above_indep) {
    arjun.standalone_elim_to_file(cnf, etof_conf, simp_conf, iconf);
  } else {
    if (verb >= 5)
      cout<< "c o skipping strong Arjun simp because size sampl_vars < " << td_at_or_above_indep << endl;
    cnf.renumber_sampling_vars_for_ganak();
  }
  /* verb_print(1, "Arjun T: " << (cpu_time()-my_time)); */
}

// Cheap residual-work proxy: after propagation, measure the largest connected
// variable set and the sum of squared component sizes. This is only a ranking
// heuristic; probes never modify the CNF or contribute to the count.
static std::pair<unsigned, uint64_t> probe_job(const vector<vector<Lit>> &cs, unsigned n,
                                               const vector<Lit> &prefix, const vector<Lit> &predicate,
                                               bool truth, uint64_t &budget) {
  vector<int8_t> values(n + 1, -1);
  auto assign = [&](Lit l) { values[l.var()] = !l.sign(); };
  for (auto l : prefix)
    assign(l);
  if (!truth || predicate.size() == 1)
    for (auto l : predicate)
      assign(truth ? l : l.neg());
  vector<unsigned> parent(n + 1), size(n + 1, 1);
  std::iota(parent.begin(), parent.end(), 0);
  auto root = [&](unsigned v) {
    while (parent[v] != v) {
      parent[v] = parent[parent[v]];
      v = parent[v];
    }
    return v;
  };
  vector<bool> active(n + 1);
  bool changed;
  do {
    changed = false;
    for (size_t i = 0; i < cs.size() + unsigned(truth && predicate.size() > 1); ++i) {
      const auto &c = i == cs.size() ? predicate : cs[i];
      if (budget < c.size() + 1) {
        budget = 0;
        return {};
      }
      budget -= c.size() + 1;
      unsigned unknown = 0;
      Lit last;
      bool satisfied = false, has_a = false, has_b = false;
      for (Lit l : c) {
        if (values[l.var()] < 0) {
          ++unknown;
          last = l;
        } else if (bool(values[l.var()]) != l.sign()) {
          satisfied = true;
          break;
        }
        if (truth && predicate.size() == 2 && i != cs.size()) {
          has_a |= l == predicate[0];
          has_b |= l == predicate[1];
        }
      }
      if (satisfied || (has_a && has_b)) continue;
      if (!unknown) return {};
      if (unknown == 1) {
        assign(last);
        changed = true;
      }
    }
  } while (changed);
  for (size_t i = 0; i < cs.size() + unsigned(truth && predicate.size() > 1); ++i) {
    const auto &c = i == cs.size() ? predicate : cs[i];
    if (budget < c.size() + 1) {
      budget = 0;
      return {};
    }
    budget -= c.size() + 1;
    if (std::any_of(c.begin(), c.end(),
                    [&](Lit l) { return values[l.var()] >= 0 && bool(values[l.var()]) != l.sign(); }))
      continue;
    if (truth && predicate.size() == 2 && i != cs.size() &&
        std::find(c.begin(), c.end(), predicate[0]) != c.end() &&
        std::find(c.begin(), c.end(), predicate[1]) != c.end())
      continue;
    unsigned first = 0;
    for (auto l : c)
      if (values[l.var()] < 0) {
        unsigned v = l.var();
        active[v] = true;
        if (!first)
          first = v;
        else {
          unsigned a = root(first), b = root(v);
          if (a != b) {
            parent[b] = a;
            size[a] += size[b];
          }
        }
      }
  }
  unsigned largest = 0;
  uint64_t total = 0;
  for (unsigned v = 1; v <= n; ++v)
    if (active[v] && parent[v] == v) {
      largest = std::max(largest, size[v]);
      total += uint64_t(size[v]) * size[v];
    }
  return {largest, total};
}

FF OuterCounter::count_with_parallel(uint8_t bits_jobs, int num_threads) {
  assert(conf.appmc_timeout < 0 && "TD-parallel not compatible with AppMC");
  assert(bits_jobs > 0);
  assert(bits_jobs <= 20);
  double const td_start_time = cpu_time();

  // Build primal graph from clauses
  TWD::Graph primal(nvars);
  for (const auto& cl : irred_cls) {
    for (auto it1 = cl.begin(); it1 != cl.end(); ++it1) {
      for (auto it2 = std::next(it1); it2 != cl.end(); ++it2) {
        uint32_t const v1 = it1->var() - 1;
        uint32_t const v2 = it2->var() - 1;
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
  auto tdec  = fc.constructTD(conf.td_steps / 3, conf.td_iters / 3,
      conf.td_band_pct, conf.td_dense_pct);

  // Find centroid
  int const centroid_id = tdec.centroid(conf.verb);
  vector<int> centroid_bag = tdec.Bags()[centroid_id];

  verb_print(1, "[par] TD width: " << tdec.width()
          << ", centroid bag size: " << centroid_bag.size()
          << ", TD time: " << cpu_time()-td_start_time << "s");

  // If centroid bag is empty or of size 1, just use regular counting
  if (centroid_bag.size() <= 1) {
    verb_print(1, "[par] Centroid bag is empty, using regular counting");
    return count_regular();
  }
  if (centroid_bag.size() < bits_jobs) {
    assert(centroid_bag.size() >= 2);
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
  for(const auto v: centroid_bag) {
    verb_print(2, "[par] var " << v+1 << " with frequency " << var_freq[v+1]);
  }

  const uint64_t num_jobs = 1ULL << bits_jobs;
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

  // Since optional indep support cannot be provided, vars with weights must be in indep support
  for (const auto& [lit, weight] : lit_weights) {
    if (weight->is_one()) continue;
    indep_support.insert(lit.var());
  }

  // Keep the first k-1 centroid bits. Pair the last one with a literal
  // from an incident clause; maximize the number of clauses discharged
  // when this OR is true. A counted-variable predicate preserves disjointness
  // even when input assignments are projected onto the sampling set.
  vector<Lit> gate;
  if (conf.mixed_jobs && !fg->weighted() && bits_jobs > 0) {
    const uint32_t pivot = centroid_bag[bits_jobs-1]+1;
    const set<int> fixed(centroid_bag.begin(), centroid_bag.begin()+bits_jobs-1);
    std::map<std::pair<Lit,Lit>, unsigned> freq;
    if (indep_support.count(pivot)) for (const auto& cl : irred_cls) {
      if (cl.size() < 3) continue;
      auto p = std::find_if(cl.begin(),cl.end(),[&](Lit l){return l.var()==pivot;});
      if (p==cl.end()) continue;
      for (Lit l : cl) {
        if (l.var()==pivot || fixed.count(l.var()-1) || !indep_support.count(l.var())) continue;
        ++freq[{*p,l}];
      }
    }
    if (conf.mixed_jobs == 1) {
      unsigned best = 0;
      for (const auto &[pair, score] : freq)
        if (score > best) {
          best = score;
          gate = {pair.first, pair.second};
        }
    } else if (!freq.empty()) {
      // Use the same prefix cubes for every candidate. With more than two
      // prefix bits use four fixed samples spanning the job range.
      uint64_t budget = 8000000;
      const unsigned prefix_jobs = 1u << (bits_jobs - 1), samples = std::min(prefix_jobs, 4u);
      auto score = [&](const vector<Lit> &pred) {
        unsigned worst = 0;
        uint64_t squares = 0;
        for (unsigned j = 0; j < samples && budget; ++j) {
          unsigned mask = samples == 1 ? 0 : j * (prefix_jobs - 1) / (samples - 1);
          vector<Lit> prefix;
          for (unsigned i = 0; i + 1 < bits_jobs; ++i)
            prefix.emplace_back(centroid_bag[i] + 1, (mask >> i) & 1);
          for (bool truth : {false, true}) {
            auto [width, work] = probe_job(irred_cls, nvars, prefix, pred, truth, budget);
            worst = std::max(worst, width);
            squares += work;
          }
        }
        return std::pair{worst, squares};
      };
      auto best = score({Lit(pivot, false)});
      const unsigned variable_width = best.first;
      if (conf.profile_jobs >= 3)
        cout << "c gate_probe variable " << pivot << " " << best.first << " " << best.second << endl;
      vector<std::pair<unsigned, std::pair<Lit, Lit>>> candidates;
      for (const auto &[pair, hits] : freq)
        candidates.emplace_back(hits, pair);
      std::stable_sort(candidates.begin(), candidates.end(),
                       [](auto &a, auto &b) { return a.first > b.first; });
      for (unsigned i = 0; i < std::min<size_t>(8, candidates.size()) && budget; ++i) {
        auto [a, b] = candidates[i].second;
        vector<Lit> candidate{a, b};
        auto rank = score(candidate);
        if (conf.profile_jobs >= 3 && budget)
          cout << "c gate_probe or " << a << " " << b << " " << rank.first << " " << rank.second << endl;
        // A true OR can leave both inputs live, where a variable decision
        // necessarily removes one. Allow this single-variable difference in
        // mode 3, then compare the aggregate post-propagation cost.
        const bool improves = conf.mixed_jobs == 3 ? rank.first <= variable_width + 1 &&
                                                         std::tie(rank.second, rank.first) <
                                                             std::tie(best.second, best.first)
                                                   : rank < best;
        if (budget && improves) {
          best = rank;
          gate = candidate;
        }
      }
    }
  }
  using JobClock = std::chrono::steady_clock;
  vector<double> job_seconds(conf.profile_jobs ? num_jobs : 0);
  vector<std::string> job_counts(conf.profile_jobs > 1 ? num_jobs : 0);
  auto worker = [&](const uint64_t num) -> FF {
    auto start = conf.profile_jobs ? JobClock::now() : JobClock::time_point{};
    const bool gate_true = ((num >> (bits_jobs-1)) & 1) == 0;
    // Each thread needs its own FieldGen copy: fg may have mutable internal state
    FG thread_fg = fg->dup();
    ArjunNS::SimplifiedCNF cnf(thread_fg);
    cnf.new_vars(nvars);
    cnf.set_sampl_vars(ganak_to_cms_vars(indep_support));
    if (thread_fg->weighted()) {
      cnf.set_weighted(true);
      for (const auto& [lit, weight] : lit_weights)
        cnf.set_lit_weight(ganak_to_cms_lit(lit), weight->dup());
    }
    for (const auto& cl : irred_cls) {
      // Under gate=true these clauses are already satisfied. When false,
      // the two unit assumptions below let propagation shorten them.
      if (!gate.empty() && gate_true &&
          std::find(cl.begin(),cl.end(),gate[0]) != cl.end() &&
          std::find(cl.begin(),cl.end(),gate[1]) != cl.end()) continue;
      cnf.add_clause(ganak_to_cms_cl(cl));
    }

    // Add unit clauses that differentiates this thread
    for (size_t i = 0; i < bits_jobs; i++) {
      if (!gate.empty() && i+1 == bits_jobs) {
        if (gate_true) cnf.add_clause(ganak_to_cms_cl(gate));
        else for (Lit l : gate) cnf.add_clause({ganak_to_cms_lit(l.neg())});
        continue;
      }
      const uint32_t var = centroid_bag[i];
      const bool sign = (num >> i) & 1;
      CMSat::Lit const unit_lit(var, sign);
      cnf.add_clause({unit_lit});
      verb_print(2, "[par] Thread " << num << " fixing var " << var
                  << " to " << (sign ? "true" : "false"));
    }
    for (const auto& [cl, lbd] : red_cls)
      cnf.add_red_clause(ganak_to_cms_cl(cl));
    run_arjun(cnf, conf.verb >= 5 ? conf.verb : 0);
    FF ret = cnf.get_multiplier_weight()->dup();
    if (!ret->is_zero()) {
      auto local_conf = conf;
      local_conf.verb = conf.verb >= 3 ? conf.verb : 0;
      auto counter = std::make_unique<Ganak>(local_conf, thread_fg);
      if (conf.verb >= 3) cnf.write_simpcnf("input_to_thread_" + std::to_string(num) + ".cnf");
      setup_ganak(cnf, *counter);
      *ret *= *counter->count();
      num_cache_lookups += counter->get_num_cache_lookups();
      {
        std::lock_guard<std::mutex> const lock(stats_mutex);
        max_cache_elems = std::max(max_cache_elems, counter->get_max_cache_elems());
        count_is_approximate |= counter->get_is_approximate();
      }
    }
    if (conf.profile_jobs > 1) {
      std::ostringstream text;
      text << *ret;
      job_counts[num] = text.str();
    }
    if (conf.profile_jobs)
      job_seconds[num] = std::chrono::duration<double>(JobClock::now()-start).count();
    return ret;
  };

  // Thread pool: maintain at most max_concurrent active threads
  std::vector<std::pair<uint64_t, std::future<FF>>> active_futures;
  uint64_t next_task = 0;
  FF total_count = fg->zero();

  // Launch initial batch of threads
  assert(num_threads > 0);
  for (uint32_t i = 0; i < (uint32_t)num_threads && next_task < num_jobs; i++) {
    active_futures.emplace_back(next_task, std::async(std::launch::async, worker, next_task));
    next_task++;
  }

  // Process results and launch new tasks as threads complete
  while (!active_futures.empty()) {
    for (size_t i = 0; i < active_futures.size(); ) {
      auto& [task_id, future] = active_futures[i];
      // Check if this future is ready
      if (future.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
        // Collect result
        FF const partial = future.get();
        *total_count += *partial;
        verb_print(2, "[par] Task " << task_id << " completed");

        // Remove from active list
        active_futures.erase(active_futures.begin() + i);

        // Launch next task if available
        if (next_task < num_jobs) {
          active_futures.emplace_back(next_task, std::async(std::launch::async, worker, next_task));
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

  if (conf.profile_jobs > 1) {
    if (!gate.empty()) cout << "c job_gate " << gate[0] << " " << gate[1] << endl;
    for (uint64_t i=0; i<num_jobs; ++i)
      cout << "c job_detail " << i << " " << job_seconds[i] << " " << job_counts[i] << endl;
  }
  if (conf.profile_jobs) {
    cout << "c job_profile " << num_jobs << " " << !gate.empty() << " "
         << *std::max_element(job_seconds.begin(),job_seconds.end()) << endl;
  }
  return total_count;
}

void OuterCounter::print_indep_distrib() const {
  cout << "c o indep/optional/none distribution: ";
  std::set<uint32_t> const all_opt_indep = opt_indep_support;
  std::set<uint32_t> const all_indep(indep_support);
  for(uint32_t i = 0; i <= nvars; i++) {
    if (all_opt_indep.count(i)) {
      if (all_indep.count(i)) cout << "I";
      else cout << "O";
    } else cout << "N";
  }
  cout << endl;
}

}
