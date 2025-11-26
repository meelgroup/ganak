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
#include "TreeDecomposition.hpp"
#include "IFlowCutter.hpp"
#include "time_mem.hpp"
#include "common.hpp"

namespace GanakInt {

FF OuterCounter::count() {
  if (conf.do_td && nvars > 100 && indep_support.size() > 10 &&  nvars <= conf.td_varlim && conf.appmc_timeout == 0)
    return count_with_td_parallel();
  else
    return count_regular();
}

FF OuterCounter::count_regular() {
  auto counter = std::make_unique<Counter>(conf, fg);
  counter->new_vars(nvars);
  counter->set_indep_support(indep_support);
  counter->set_optional_indep_support(opt_indep_support);
  counter->set_generators(generators);

  for(const auto& cl : irred_cls) counter->add_irred_cl(cl);
  counter->end_irred_cls();
  for(const auto& p : red_cls) counter->add_red_cl(p.first, (int)p.second);
  auto ret = counter->outer_count();

  assert(num_cache_lookups == 0);
  assert(max_cache_elems == 0);
  num_cache_lookups = counter->get_stats().num_cache_look_ups;
  max_cache_elems = counter->get_cache()->get_max_num_entries();
  return ret;
}

FF OuterCounter::count_with_td_parallel() {
  assert(conf.appmc_timeout == 0 && "TD-parallel not compatible with AppMC");
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
  verb_print(1, "[td-par] nodes: " << nodes << " nvars: " << nvars << " edges: " << primal.numEdges());

  // Run FlowCutter to get tree decomposition
  TWD::IFlowCutter fc(primal.numNodes(), primal.numEdges(), conf.verb);
  fc.importGraph(primal);

  // Compute TD with reduced steps to avoid timeout
  auto tdec  = fc.constructTD(conf.td_steps / 10, conf.td_iters / 10);

  // Find centroid
  //
  /* const auto& bags = tdec.Bags(); */
  /* const auto& adj = tdec.get_adj_list(); */
  int centroid_id = tdec.centroid(primal.numNodes(), conf.verb);
  auto& centroid_bag = tdec.Bags()[centroid_id];

  verb_print(1, "[td-par] TD width: " << tdec.width()
          << ", centroid bag size: " << centroid_bag.size()
          << ", TD time: " << td_start_time-cpu_time() << "s");

  // If centroid bag is empty or very small, just use regular counting
  if (centroid_bag.size() <= 2) {
    if (conf.verb >= 1) {
      std::cout << "c o [td-par] Centroid bag too small, using regular counting" << std::endl;
    }
    return count_regular();
  }

  // Launch parallel threads
  uint64_t nthreads = 4;//std::thread::hardware_concurrency();
  uint8_t num_bits = 2;
  assert(1<<num_bits == nthreads);

  if (conf.verb >= 1) {
    std::cout << "c o [td-par] Launching " << nthreads
              << " parallel threads for centroid variables: ";
    for (uint32_t i = 0; i < 4; i++) std::cout << (centroid_bag[i] + 1) << " ";
    std::cout << std::endl;
  }

  std::vector<std::future<FF>> futures;
  auto worker = [&](uint64_t num) -> FF {
    // Create a new Counter for this configuration
    Counter local_counter(conf, fg->dup());
    local_counter.new_vars(nvars);
    local_counter.set_indep_support(indep_support);
    local_counter.set_optional_indep_support(opt_indep_support);

    for (const auto& [lit, weight] : lit_weights)
      local_counter.set_lit_weight(lit, weight->dup());

    for (const auto& cl : irred_cls) local_counter.add_irred_cl(cl);

    // Add unit clauses for centroid variable assignment
    for (size_t i = 0; i < num_bits; i++) {
      uint32_t var = centroid_bag[i] + 1; // Convert from 0-indexed to 1-indexed
      bool sign = (num >> i) & 1;
      Lit unit_lit(var, sign);
      local_counter.add_irred_cl({unit_lit});
    }
    local_counter.end_irred_cls();

    for (const auto& [cl, lbd] : red_cls) local_counter.add_red_cl(cl, lbd);

    if (!generators.empty()) local_counter.set_generators(generators);

    auto ret = local_counter.outer_count();
    num_cache_lookups += local_counter.get_stats().num_cache_look_ups;
    max_cache_elems = std::max(max_cache_elems, local_counter.get_cache()->get_max_num_entries());
    return ret;
  };

  for (uint32_t t = 0; t < nthreads; t++) {
    futures.push_back(std::async(std::launch::async, worker, t));
  }

  // Collect results
  FF total_count = fg->zero();
  for (auto& future : futures) {
    FF partial = future.get();
    *total_count += *partial;
  }

  if (conf.verb >= 1) {
    std::cout << "c o [td-par] Parallel counting completed, total time: "
              << (cpu_time() - td_start_time) << "s" << std::endl;
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
