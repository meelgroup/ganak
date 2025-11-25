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
#include <thread>
#include <mutex>
#include <future>
#include "TreeDecomposition.hpp"
#include "IFlowCutter.hpp"
#include "time_mem.hpp"

namespace GanakInt {

FF OuterCounter::count() {
  // Try TD-based parallel counting if enabled and feasible
  if (conf.do_td && nvars > 5 && nvars <= conf.td_varlim) {
    return count_with_td_parallel();
  }

  // Fall back to regular counting
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
  return counter->outer_count();
}

FF OuterCounter::count_with_td_parallel() {
  double td_start_time = cpu_time();

  // Build primal graph from clauses
  TWD::Graph primal(nvars);

  for (const auto& cl : irred_cls) {
    for (size_t i = 0; i < cl.size(); i++) {
      for (size_t j = i + 1; j < cl.size(); j++) {
        uint32_t v1 = cl[i].var() - 1;
        uint32_t v2 = cl[j].var() - 1;
        if (v1 != v2) {
          primal.addEdge(v1, v2);
        }
      }
    }
  }

  // Also add edges from red clauses
  for (const auto& cl_pair : red_cls) {
    const auto& cl = cl_pair.first;
    for (size_t i = 0; i < cl.size(); i++) {
      for (size_t j = i + 1; j < cl.size(); j++) {
        uint32_t v1 = cl[i].var() - 1;
        uint32_t v2 = cl[j].var() - 1;
        if (v1 != v2) {
          primal.addEdge(v1, v2);
        }
      }
    }
  }

  if (conf.verb >= 1) {
    std::cout << "c o [td-par] Built primal graph with " << primal.numNodes()
              << " nodes and " << primal.numEdges() << " edges" << std::endl;
  }

  // Run FlowCutter to get tree decomposition
  TWD::IFlowCutter fc(primal.numNodes(), primal.numEdges(), conf.verb >= 2 ? 1 : 0);
  fc.importGraph(primal);

  // Compute TD with reduced steps to avoid timeout
  auto td = fc.constructTD(conf.td_steps / 10, conf.td_iters / 10);

  double td_time = cpu_time() - td_start_time;

  // Check if TD timed out or failed
  if (td.width() <= 1000) {
    if (conf.verb >= 1) {
      std::cout << "c o [td-par] TD computation timed out or failed (time: "
                << td_time << "s), falling back to regular counting" << std::endl;
    }
    return count_regular();
  }

  // Find centroid
  int centroid_id = td.centroid(primal.numNodes(), conf.verb >= 2 ? 1 : 0);
  auto& centroid_bag = td.Bags()[centroid_id];

  if (conf.verb >= 1) {
    std::cout << "c o [td-par] TD width: " << td.width()
              << ", centroid bag size: " << centroid_bag.size()
              << ", TD time: " << td_time << "s" << std::endl;
  }

  // Check if centroid bag is too large for parallel splitting
  const uint32_t max_centroid_size = 20;
  if (centroid_bag.size() > max_centroid_size) {
    if (conf.verb >= 1) {
      std::cout << "c o [td-par] Centroid bag too large (" << centroid_bag.size()
                << " > " << max_centroid_size
                << "), falling back to regular counting" << std::endl;
    }
    return count_regular();
  }

  // If centroid bag is empty or very small, just use regular counting
  if (centroid_bag.size() <= 2) {
    if (conf.verb >= 1) {
      std::cout << "c o [td-par] Centroid bag too small, using regular counting" << std::endl;
    }
    return count_regular();
  }

  // Launch parallel threads
  uint64_t num_configs = 1ULL << centroid_bag.size();

  if (conf.verb >= 1) {
    std::cout << "c o [td-par] Launching " << num_configs
              << " parallel threads for centroid variables: ";
    for (auto v : centroid_bag) std::cout << (v + 1) << " ";
    std::cout << std::endl;
  }

  // Determine number of hardware threads
  uint32_t num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) num_threads = 4; // Default fallback

  // We'll process configurations in batches
  std::vector<std::future<FF>> futures;
  std::atomic<uint64_t> next_config(0);
  std::mutex result_mutex;

  auto worker = [&](uint64_t start_config, uint64_t end_config) -> FF {
    FF local_sum = fg->zero();

    for (uint64_t config = start_config; config < end_config; config++) {
      // Create a new Counter for this configuration
      Counter local_counter(conf, fg->dup());
      local_counter.new_vars(nvars);
      local_counter.set_indep_support(indep_support);
      local_counter.set_optional_indep_support(opt_indep_support);

      // Set literal weights
      for (const auto& [lit, weight] : lit_weights) {
        local_counter.set_lit_weight(lit, weight->dup());
      }

      // Add irredundant clauses
      for (const auto& cl : irred_cls) {
        local_counter.add_irred_cl(cl);
      }

      // Add unit clauses for centroid variable assignment
      for (size_t i = 0; i < centroid_bag.size(); i++) {
        uint32_t var = centroid_bag[i] + 1; // Convert from 0-indexed to 1-indexed
        bool sign = (config >> i) & 1;
        Lit unit_lit(var, sign);
        local_counter.add_irred_cl({unit_lit});
      }

      local_counter.end_irred_cls();

      // Add redundant clauses
      for (const auto& [cl, lbd] : red_cls) {
        local_counter.add_red_cl(cl, lbd);
      }

      // Set generators if any
      if (!generators.empty()) {
        local_counter.set_generators(generators);
      }

      // Count
      FF thread_result = local_counter.outer_count();
      *local_sum += *thread_result;
    }

    return local_sum;
  };

  // Divide work among threads
  uint64_t configs_per_thread = (num_configs + num_threads - 1) / num_threads;

  for (uint32_t t = 0; t < num_threads; t++) {
    uint64_t start = t * configs_per_thread;
    uint64_t end = std::min(start + configs_per_thread, num_configs);

    if (start >= num_configs) break;

    futures.push_back(std::async(std::launch::async, worker, start, end));
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

}
