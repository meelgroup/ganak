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

#include "structures.h"
#include "comp_types/cacheable_comp.h"
#include "primitive_types.h"
#include "boundedqueue.h"
#include "counter_config.h"

using std::vector;
using std::cout;
using std::endl;

class Instance;
class Counter;
class ComponentCache;

class DataAndStatistics {
public:
  DataAndStatistics (const Instance* _inst, CounterConfiguration& _conf): conf(_conf) {
    inst = _inst;
    cache_hits_misses_q.clearAndResize(10000);
    comp_size_times_depth_q.clearAndResize(10000);
  }
  CounterConfiguration& conf;
  uint64_t maximum_cache_size_bytes_ = 0;
  uint64_t numcachedec_ = 0;

  uint64_t num_unit_irred_clauses_ = 0;
  uint64_t num_long_irred_clauses_ = 0;
  uint64_t num_binary_irred_clauses_ = 0;

  uint64_t num_unit_red_clauses_ = 0;
  uint64_t num_binary_red_clauses_ = 0;

  // Clause db management
  uint64_t reduceDBs = 0;
  uint32_t cls_deleted_since_compaction = 0;
  uint32_t compactions = 0;
  uint32_t cls_removed = 0;

  /// number of all decisions made
  uint64_t decisions = 0;
  /// number of all implications derived
  uint64_t num_implications_ = 0;
  // number of all failed literal detections
  uint64_t num_failed_literals_detected_ = 0;
  uint64_t num_failed_bprop_literals_failed = 0;
  uint64_t num_failed_lit_tests_ = 0;

  // number of all conflicts occurred
  uint64_t conflicts = 0;

  // number of clauses overall learned
  uint32_t num_clauses_learned_ = 0;
  uint64_t uip_not_added = 0;
  uint64_t uip_cls = 0;
  uint64_t final_cl_sz = 0;
  uint64_t uip_lits_ccmin = 0;
  uint64_t rem_lits_with_bins = 0;
  uint32_t rem_lits_tried = 0;

  uint64_t  orig_uip_lits = 0;
  uint64_t last_restart_num_conflicts = 0;
  uint64_t last_restart_num_decisions = 0;

  uint64_t saved_uip = 0;
  uint64_t saved_uip_used = 0;
  uint64_t saved_uip_thrown = 0;

  uint64_t  saved_uip_used_falsified = 0;
  uint64_t  saved_uip_used_asserting = 0;
  uint64_t  saved_uip_used_sat_or_unk = 0;

  // Vivification
  uint64_t vivif_tried = 0;
  uint64_t vivif_tried_cl = 0;
  uint64_t vivif_cl_minim = 0;
  uint64_t vivif_lit_rem = 0;

  /* cache statistics */
  uint64_t num_cache_hits_ = 0;
  uint64_t num_cache_look_ups_ = 0;
  uint64_t last_restart_num_cache_look_ups = 0;
  uint64_t sum_cache_hit_sizes_ = 0;
  uint64_t sum_cache_store_sizes_ = 0;
  bqueue<uint32_t> cache_hits_misses_q;

  uint64_t num_cached_comps_ = 0;
  uint64_t total_num_cached_comps_ = 0;
  uint64_t cache_pollutions_removed = 0;
  uint64_t cache_pollutions_called = 0;

  bqueue<uint64_t, double> comp_size_times_depth_q;

  // Lookahead
  uint64_t lookaheads = 0;
  uint64_t lookahead_computes = 0;

  // the number of bytes occupied by all comps
  uint64_t sum_bytes_cached_comps_ = 0;

  uint64_t cache_infrastructure_bytes_memory_usage_ = 0;

  const Instance* inst;

  bool cache_full(uint64_t empty_size) {
    return (cache_bytes_memory_usage() - empty_size) >= maximum_cache_size_bytes_;
  }

  uint64_t cache_bytes_memory_usage() const {
    return cache_infrastructure_bytes_memory_usage_
           + sum_bytes_cached_comps_;
  }

  void incorporate_cache_store(const CacheableComponent &ccomp, const uint32_t comp_nvars) {
    sum_bytes_cached_comps_ += ccomp.SizeInBytes();
    sum_cache_store_sizes_ += comp_nvars;
    num_cached_comps_++;
    total_num_cached_comps_++;
  }

  void incorporate_cache_erase(const CacheableComponent &ccomp){
    sum_bytes_cached_comps_ -= ccomp.SizeInBytes();
    num_cached_comps_--;
  }

  void incorporate_cache_hit(const uint32_t comp_nvars){
    num_cache_hits_++;
    sum_cache_hit_sizes_ += comp_nvars;
  }

  double implicitBCP_miss_rate() const {
      if(num_failed_lit_tests_ == 0) return 0.0;
      return (num_failed_lit_tests_ - num_failed_literals_detected_) / (double) num_failed_lit_tests_;
  }
  uint64_t num_irred_clauses() const {
    return num_long_irred_clauses_ + num_binary_irred_clauses_ + num_unit_irred_clauses_;
  }

  void incorporateConflictClauseData(const vector<Lit> &clause) {
    if (clause.size() == 1) num_unit_red_clauses_++;
    else if (clause.size() == 2) num_binary_red_clauses_++;
  }

  void incorporateIrredClauseData(const vector<Lit> &clause) {
    if (clause.size() == 1) num_unit_irred_clauses_++;
    else if (clause.size() == 2) num_binary_irred_clauses_++;
    else num_long_irred_clauses_++;
  }

  void printShort(const Counter* counter, const ComponentCache* cache_) const;
  void printShortFormulaInfo() const {
    verb_print(1, "irred cls (all/long/bin/unit): "
      << num_irred_clauses() << "/" << num_long_irred_clauses_
      << "/" << num_binary_irred_clauses_ << "/" << num_unit_irred_clauses_);
  }

  double getAvgComponentHitSize() const {
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

  long double getAvgCacheHitSize() const {
    if(num_cache_hits_ == 0) return 0.0L;
    return sum_cache_hit_sizes_ / (long double) num_cache_hits_;
  }

  long double getAvgCacheStoreSize() const {
    if(total_num_cached_comps_ == 0) return 0.0L;
    return sum_cache_store_sizes_ / (long double) total_num_cached_comps_;
  }
};
