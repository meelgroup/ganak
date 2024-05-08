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

#include <string>
#include <cstdint>
#include <vector>
#include <gmpxx.h>

#include "structures.hpp"
#include "comp_types/cacheable_comp.hpp"
#include "primitive_types.hpp"
#include "boundedqueue.hpp"
#include "counter_config.hpp"

using std::vector;
using std::cout;
using std::endl;

class Instance;
class Counter;
class CompCache;

class DataAndStatistics {
public:
  DataAndStatistics (const Instance* _inst, CounterConfiguration& _conf): conf(_conf) {
    inst = _inst;
    comp_size_times_depth_q.clearAndResize(10000);
  }
  CounterConfiguration& conf;
  uint64_t maximum_cache_size_bytes_ = 0;

  uint64_t num_binary_irred_clauses_ = 0;

  uint64_t num_binary_red_clauses_ = 0;

  // Cubes
  uint64_t num_cubes_orig = 0;
  uint64_t num_cubes_final = 0;
  uint64_t num_cubes_symm = 0;
  uint64_t cube_lit_extend = 0;
  uint64_t cube_lit_rem = 0;

  // Clause db management
  uint64_t reduce_db = 0;
  uint32_t cls_deleted_since_compaction = 0;
  uint32_t compactions = 0;
  uint32_t cls_removed = 0;

  /// number of all decisions made
  uint64_t decisions = 0;

  // number of all conflicts occurred
  uint64_t conflicts = 0;

  // number of clauses overall learned
  uint32_t num_clauses_learned_ = 0;
  uint64_t uip_cls = 0;
  uint64_t final_cl_sz = 0;
  uint64_t uip_lits_ccmin = 0;
  uint64_t rem_lits_with_bins = 0;
  uint32_t rem_lits_tried = 0;

  uint64_t  orig_uip_lits = 0;
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

  // buddy stats
  uint64_t buddy_called = 0;
  uint64_t buddy_num_bin_cls = 0;
  uint64_t buddy_num_long_cls = 0;

  /* cache statistics */
  uint64_t num_cache_hits_ = 0;
  uint64_t num_cache_dels_ = 0;
  uint64_t num_cache_look_ups_ = 0;
  uint64_t last_restart_num_cache_look_ups = 0;
  uint64_t sum_cache_hit_sizes_ = 0;
  uint64_t sum_cache_store_sizes_ = 0;

  // Components
  uint64_t comp_sorts = 0;
  uint64_t comp_sizes = 0;

  uint64_t num_cached_comps_ = 0;
  uint64_t total_num_cached_comps_ = 0;
  uint64_t cache_pollutions_removed = 0;
  uint64_t cache_pollutions_called = 0;

  bqueue<double, double> comp_size_times_depth_q;

  // Lookahead
  uint64_t lookaheads = 0;
  uint64_t lookahead_computes = 0;

  // the number of bytes occupied by all comps
  uint64_t sum_bignum_bytes = 0;
  uint64_t cache_infrastructure_bytes_memory_usage_ = 0;

  const Instance* inst;

  bool cache_full(const uint64_t empty_size, uint64_t extra_will_be_added) {
    return (cache_bytes_memory_usage() - empty_size + extra_will_be_added)
      >= maximum_cache_size_bytes_;
  }

  uint64_t cache_bytes_memory_usage() const {
    return cache_infrastructure_bytes_memory_usage_
           + sum_bignum_bytes;
  }

  void incorporate_cache_store(const CacheableComp &ccomp, const uint32_t comp_nvars) {
    sum_bignum_bytes += ccomp.bignum_bytes();
    sum_cache_store_sizes_ += comp_nvars;
    num_cached_comps_++;
    total_num_cached_comps_++;
  }

  void incorporate_cache_erase(const CacheableComp &ccomp){
    sum_bignum_bytes -= ccomp.bignum_bytes();
    num_cached_comps_--;
  }

  void incorporate_cache_hit(const uint32_t comp_nvars){
    num_cache_hits_++;
    sum_cache_hit_sizes_ += comp_nvars;
  }

  void incorporateConflictClauseData(const vector<Lit> &clause) {
    if (clause.size() == 2) num_binary_red_clauses_++;
  }

  void incorporateIrredClauseData(const vector<Lit> &clause) {
    if (clause.size() == 1) return;
    if (clause.size() == 2) num_binary_irred_clauses_++;
  }

  void print_short(const Counter* counter, const CompCache* cache) const;
  void print_short_formula_info(const Counter* counter) const;

  double get_avg_comp_hit_size() const {
    if (num_cache_hits_ == 0) return 0.0L;
    return (double)sum_cache_hit_sizes_ / (double) num_cache_hits_;
  }

  uint64_t cached_comp_count() const { return num_cached_comps_; }
  uint64_t cache_hits() const { return num_cache_hits_; }

  double cache_miss_rate() const {
    if(num_cache_look_ups_ == 0) return 0.0;
    return (num_cache_look_ups_ - num_cache_hits_)
        / (double) num_cache_look_ups_;
  }

  long double get_avg_cache_store_sz() const {
    if(num_cache_hits_ == 0) return 0.0L;
    return sum_cache_hit_sizes_ / (long double) num_cache_hits_;
  }

  long double get_avg_cache_store_size() const {
    if(total_num_cached_comps_ == 0) return 0.0L;
    return sum_cache_store_sizes_ / (long double) total_num_cached_comps_;
  }
};
