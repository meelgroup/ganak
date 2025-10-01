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

#pragma once

#include <cstdint>
#include <vector>
#include <gmpxx.h>

#include "counter_config.hpp"
#include "comp_cache_if.hpp"
#include "common.hpp"
#include "lit.hpp"

using std::vector;
using std::cout;
using std::endl;

namespace GanakInt {

class Counter;
class CompCacheIF;

class DataAndStatistics {
public:
  DataAndStatistics (const CounterConfiguration& _conf, const FG& _fg): conf(_conf), fg(_fg->dup()) {}
  const CounterConfiguration& conf;
  FG fg;
  uint64_t num_bin_irred_cls = 0;
  uint64_t num_bin_red_cls = 0;

  // Cubes
  uint64_t num_cubes_orig = 0;
  uint64_t num_cubes_final = 0;
  uint64_t num_cubes_symm = 0;
  uint64_t cube_lit_extend = 0;
  uint64_t cube_lit_rem = 0;

  // Clause db management
  uint64_t reduce_db = 0;
  uint32_t cls_deleted_since_compaction = 0;
  uint32_t cls_removed = 0;

  /// number of all decisions made
  uint64_t decisions = 0;

  // number of all conflicts occurred
  uint64_t conflicts = 0;
  uint64_t learnt_cls_added = 0;

  // number of clauses overall learned
  uint64_t uip_cls = 0;
  uint64_t final_cl_sz = 0;
  uint64_t uip_lits_ccmin = 0;
  uint64_t rem_lits_with_bins = 0;
  uint32_t rem_lits_tried = 0;

  uint64_t orig_uip_lits = 0;
  uint64_t last_restart_num_conflicts = 0;
  uint64_t last_restart_num_decisions = 0;
  uint32_t num_restarts = 0;

  // Vivification
  uint64_t vivif_tried = 0;
  uint64_t vivif_tried_cl = 0;
  uint64_t vivif_cl_minim = 0;
  uint64_t vivif_lit_rem = 0;

  // Toplevel probing
  uint64_t toplevel_probe_runs = 0;
  uint64_t toplevel_probe_fail = 0;
  uint64_t toplevel_bothprop_fail = 0;

  // Subsumption
  uint64_t subsume_runs = 0;
  uint64_t subsumed_bin_irred_cls = 0;
  uint64_t subsumed_bin_red_cls = 0;
  uint64_t subsumed_long_irred_cls = 0;
  uint64_t subsumed_long_red_cls = 0;

  // SAT solver stats
  uint64_t sat_called = 0;
  uint64_t sat_found_sat = 0;
  uint64_t sat_found_unsat = 0;
  uint64_t sat_conflicts = 0;
  uint64_t sat_rst = 0;

  // buddy stats
  uint64_t buddy_called = 0;
  uint64_t buddy_unsat = 0;
  uint64_t buddy_num_vars = 0;
  uint64_t buddy_num_bin_cls = 0;
  uint64_t buddy_num_long_cls = 0;
  uint64_t buddy_max_bin_cls = 0;
  uint64_t buddy_max_long_cls = 0;
  uint64_t buddy_max_num_vars = 0;

  /* cache statistics */
  uint64_t max_cache_size_bytes = 0;
  uint64_t num_cache_hits = 0;
  uint64_t num_cache_dels = 0;
  uint64_t num_cache_look_ups = 0;
  uint64_t last_restart_num_cache_look_ups = 0;
  uint64_t sum_cache_hit_sizes = 0;
  uint64_t sum_cache_store_sizes = 0;

  // Components
  uint64_t comp_sorts = 0;
  uint64_t comp_sizes = 0;
  uint64_t comps_reset = 0;
  uint64_t comps_non_reset = 0;
  uint32_t comps_recorded = 0;

  uint64_t num_cached_comps = 0;
  uint64_t total_num_cached_comps = 0;
  uint64_t cache_pollutions_removed = 0;
  uint64_t cache_pollutions_called = 0;

  // Lookahead
  uint64_t lookaheads = 0;
  uint64_t lookahead_computes = 0;

  // the number of bytes occupied by all comps
  uint64_t sum_extra_bytes = 0;
  uint64_t cache_infra_bytes_mem_usage = 0;

  bool cache_full(const uint64_t empty_size, uint64_t extra_will_be_added) {
    return (cache_bytes_memory_usage() - empty_size + extra_will_be_added)
      >= max_cache_size_bytes;
  }

  uint64_t cache_bytes_memory_usage() const {
    return cache_infra_bytes_mem_usage
           + sum_extra_bytes;
  }

  void incorporate_cache_store(const uint64_t& extra_bytes, const uint32_t comp_nvars) {
    sum_extra_bytes += extra_bytes;
    sum_cache_store_sizes += comp_nvars;
    num_cached_comps++;
    total_num_cached_comps++;
  }

  void incorporate_cache_erase(const uint64_t extra_bytes){
    sum_extra_bytes -= extra_bytes;
    num_cached_comps--;
  }

  void incorporate_cache_hit(const uint32_t comp_nvars){
    num_cache_hits++;
    sum_cache_hit_sizes += comp_nvars;
  }

  void incorporateIrredClauseData(const vector<Lit>& clause) {
    if (clause.size() == 1) return;
    if (clause.size() == 2) num_bin_irred_cls++;
  }

  void print_short(const Counter* counter, const std::unique_ptr<CompCacheIF>& cache) const;
  void print_short_formula_info(const Counter* counter) const;

  double get_avg_comp_hit_size() const {
    if (num_cache_hits == 0) return 0.0L;
    return (double)sum_cache_hit_sizes / (double) num_cache_hits;
  }

  double cache_miss_rate() const {
    if(num_cache_look_ups == 0) return 0.0;
    return (double)(num_cache_look_ups - num_cache_hits)
        / (double) num_cache_look_ups;
  }

  long double get_avg_cache_store_sz() const {
    if(num_cache_hits == 0) return 0.0L;
    return sum_cache_hit_sizes / (long double) num_cache_hits;
  }

  long double get_avg_cache_store_size() const {
    if(total_num_cached_comps == 0) return 0.0L;
    return sum_cache_store_sizes / (long double) total_num_cached_comps;
  }
};

}
