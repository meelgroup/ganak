/*
 * statistics.h
 *
 *  Created on: Feb 13, 2013
 *      Author: mthurley
 */

#pragma once

#include <string>
#include <cstdint>
#include <vector>
#include <gmpxx.h>

#include "structures.h"
#include "comp_types/cacheable_comp.h"
#include "primitive_types.h"
#include "boundedqueue.h"

using std::vector;
using std::cout;
using std::endl;

class Instance;
class Counter;
class ComponentCache;

class DataAndStatistics {
public:
  DataAndStatistics (const Instance* _inst) {
    inst = _inst;
    cache_hits_misses_q.clearAndResize(10000);
    comp_size_per_depth_q.clearAndResize(10000);
  }
  uint64_t maximum_cache_size_bytes_ = 0;
  uint64_t numcachedec_ = 0;

  uint64_t num_unit_irred_clauses_ = 0;
  uint64_t num_long_irred_clauses_ = 0;
  uint64_t num_binary_irred_clauses_ = 0;

  uint64_t num_unit_red_clauses_ = 0;
  uint64_t num_long_red_clauses_ = 0;
  uint64_t num_binary_red_clauses_ = 0;

  uint64_t times_conflict_clauses_cleaned_ = 0;

  /// number of all decisions made
  uint64_t num_decisions_ = 0;
  /// number of all implications derived
  uint64_t num_implications_ = 0;
  // number of all failed literal detections
  uint64_t num_failed_literals_detected_ = 0;
  uint64_t num_failed_bprop_literals_failed = 0;
  uint64_t num_failed_lit_tests_ = 0;

  // number of all conflicts occurred
  uint64_t num_conflicts_ = 0;

  // number of clauses overall learned
  uint32_t num_clauses_learned_ = 0;

  // restarts
  uint64_t last_restart_num_decisions = 0;

  /* cache statistics */
  uint64_t num_cache_hits_ = 0;
  uint64_t num_cache_look_ups_ = 0;
  uint64_t last_restart_num_cache_look_ups = 0;
  uint64_t sum_cache_hit_sizes_ = 0;
  bqueue<uint32_t> cache_hits_misses_q;

  uint64_t num_cached_comps_ = 0;
  uint64_t total_num_cached_comps_ = 0;
  uint64_t sum_size_cached_comps_ = 0;
  uint64_t cache_pollutions_removed = 0;
  uint64_t cache_pollutions_called = 0;

  bqueue<double, double> comp_size_per_depth_q;


  // Lookahead
  uint64_t lookaheads = 0;
  uint64_t lookahead_computes = 0;

  // the number of bytes occupied by all comps
  uint64_t sum_bytes_cached_comps_ = 0;
  // the same number, summing over all comps ever stored
  uint64_t overall_bytes_comps_stored_ = 0;

  uint64_t sys_overhead_sum_bytes_cached_comps_ = 0;
    // the same number, summing over all comps ever stored
  uint64_t sys_overhead_overall_bytes_comps_stored_ = 0;

  uint64_t cache_infrastructure_bytes_memory_usage_ = 0;


  uint64_t overall_num_cache_stores_ = 0;
  const Instance* inst;

  bool cache_full(){
    return cache_bytes_memory_usage() >= maximum_cache_size_bytes_;
  }

  uint64_t cache_bytes_memory_usage() const {
    return cache_infrastructure_bytes_memory_usage_
           + sum_bytes_cached_comps_;
  }

  uint64_t overall_cache_bytes_memory_stored() const {
      return cache_infrastructure_bytes_memory_usage_
             + overall_bytes_comps_stored_;
    }

  void incorporate_cache_store(CacheableComponent &ccomp, const BPCSizes& sz){
    sum_bytes_cached_comps_ += ccomp.SizeInBytes(sz);
    sum_size_cached_comps_ += ccomp.nVars(sz);
    num_cached_comps_++;
    total_num_cached_comps_++;
    overall_bytes_comps_stored_ += ccomp.SizeInBytes(sz);
    overall_num_cache_stores_ += ccomp.nVars(sz);
    sys_overhead_sum_bytes_cached_comps_ += ccomp.sys_overhead_SizeInBytes(sz);
    sys_overhead_overall_bytes_comps_stored_ += ccomp.sys_overhead_SizeInBytes(sz);
  }

  void incorporate_cache_erase(CacheableComponent &ccomp, const BPCSizes& sz){
    sum_bytes_cached_comps_ -= ccomp.SizeInBytes(sz);
    sum_size_cached_comps_ -= ccomp.nVars(sz);
    num_cached_comps_--;
    sys_overhead_sum_bytes_cached_comps_ -= ccomp.sys_overhead_SizeInBytes(sz);
  }

  void incorporate_cache_hit(CacheableComponent &ccomp, const BPCSizes& sz){
      num_cache_hits_++;
      sum_cache_hit_sizes_ += ccomp.nVars(sz);
  }
  uint64_t cache_MB_memory_usage() const {
      return cache_bytes_memory_usage() / 1000000;
  }

  double implicitBCP_miss_rate() const {
      if(num_failed_lit_tests_ == 0) return 0.0;
      return (num_failed_lit_tests_ - num_failed_literals_detected_) / (double) num_failed_lit_tests_;
  }
  uint64_t num_irred_clauses() const {
    return num_long_irred_clauses_ + num_binary_irred_clauses_ + num_unit_irred_clauses_;
  }

  uint64_t clause_deletion_interval() const {
    return 10000 + 10 * times_conflict_clauses_cleaned_;
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

  void printShort(const Counter* solver, const ComponentCache* cache_) const;
  void printShortFormulaInfo() const {
    cout << "c irred cls (all/long/bin/unit): "
      << num_irred_clauses() << "/" << num_long_irred_clauses_
      << "/" << num_binary_irred_clauses_ << "/" << num_unit_irred_clauses_ << endl;
  }

  double getAvgComponentSize() const {
    if (num_cached_comps_ == 0) return 1.0L;
    return sum_size_cached_comps_ / (double) num_cached_comps_;
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
};
