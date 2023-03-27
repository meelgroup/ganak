/*
 * statistics.h
 *
 *  Created on: Feb 13, 2013
 *      Author: mthurley
 */

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <string>
#include <cstdint>
#include <vector>
#include <gmpxx.h>

#include "structures.h"
#include "comp_types/cacheable_comp.h"
#include "primitive_types.h"

using std::vector;
using std::cout;
using std::endl;

class Instance;

class DataAndStatistics {
public:
  DataAndStatistics (const Instance* _inst) { inst = _inst; }
  double time_elapsed_ = 0.0;
  uint64_t maximum_cache_size_bytes_ = 0;

  SOLVER_StateT exit_state_ = NO_STATE;
  // different variable counts
  // number of variables  and clauses before preprocessing
  uint64_t num_original_clauses_ = 0;

  // number of variables remaining
  uint64_t nVars_ = 0;
  // number of variables that actually occurs in clauses

  /// different clause counts

  // number of clauses after preprocessing
  uint64_t num_long_clauses_ = 0;
  uint64_t num_binary_clauses_ = 0;

  uint64_t num_long_conflict_clauses_ = 0;
  uint64_t num_binary_conflict_clauses_ = 0;

  uint64_t times_conflict_clauses_cleaned_ = 0;

  uint64_t num_unit_clauses_ = 0;
  /// number of all decisions made
  uint64_t num_decisions_ = 0;
  /// number of all implications derived
  uint64_t num_implications_ = 0;
  // number of all failed literal detections
  uint64_t num_failed_literals_detected_ = 0;
  uint64_t num_failed_lit_tests_ = 0;
  // number of all conflicts occurred
  uint64_t num_conflicts_ = 0;

  // number of clauses overall learned
  uint32_t num_clauses_learned_ = 0;

  uint32_t last_restart_decisions = 0;
  uint32_t num_restarts = 0;
  uint32_t next_restart = 100000000U;
  uint32_t next_restart_diff = 1000;

  /* cache statistics */
  uint64_t num_cache_hits_ = 0;
  uint64_t num_cache_look_ups_ = 0;
  uint64_t sum_cache_hit_sizes_ = 0;

  uint64_t num_cached_comps_ = 0;
  uint64_t total_num_cached_comps_ = 0;
  uint64_t sum_size_cached_comps_ = 0;

  // the number of bytes occupied by all
  // comps
  uint64_t sum_bytes_cached_comps_ = 0;
  // the same number, summing over all comps ever stored
  uint64_t overall_bytes_comps_stored_ = 0;

  // the above numbers, but without any overhead,
  // counting only the pure data size of the comps - without model counts
  uint64_t sum_bytes_pure_cached_comp_data_ = 0;
  // the same number, summing over all comps ever stored
  uint64_t overall_bytes_pure_stored_comp_data_ = 0;


  uint64_t sys_overhead_sum_bytes_cached_comps_ = 0;
    // the same number, summing over all comps ever stored
  uint64_t sys_overhead_overall_bytes_comps_stored_ = 0;

  uint64_t cache_infrastructure_bytes_memory_usage_ = 0;


  uint64_t overall_num_cache_stores_ = 0;
  const Instance* inst;

  void print_cache_state() {
    cout << "c printing cache state " << endl;
    cout << "c " <<  cache_bytes_memory_usage() <<" "<< maximum_cache_size_bytes_ << endl;
  }

  bool cache_full(){
    return cache_bytes_memory_usage() >= maximum_cache_size_bytes_;
  }

  uint64_t cache_bytes_memory_usage(){
    return cache_infrastructure_bytes_memory_usage_
           + sum_bytes_cached_comps_;
  }

  uint64_t overall_cache_bytes_memory_stored(){
      return cache_infrastructure_bytes_memory_usage_
             + overall_bytes_comps_stored_;
    }

  void incorporate_cache_store(CacheableComponent &ccomp, bool pccflag){
    if (pccflag){
      sum_bytes_cached_comps_ += ccomp.SizeInBytes_CLHASH();
    }
    else{
      sum_bytes_cached_comps_ += ccomp.SizeInBytes();
    }
    sum_size_cached_comps_ += ccomp.nVars();
    num_cached_comps_++;
    total_num_cached_comps_++;
    overall_bytes_comps_stored_ += ccomp.SizeInBytes();
    overall_num_cache_stores_ += ccomp.nVars();
    sys_overhead_sum_bytes_cached_comps_ += ccomp.sys_overhead_SizeInBytes();
    sys_overhead_overall_bytes_comps_stored_ += ccomp.sys_overhead_SizeInBytes();


    sum_bytes_pure_cached_comp_data_ += ccomp.data_only_byte_size();
    overall_bytes_pure_stored_comp_data_ += ccomp.data_only_byte_size();
  }
  void incorporate_cache_erase(CacheableComponent &ccomp, bool pccflag){
    if (pccflag){
      sum_bytes_cached_comps_ -= ccomp.SizeInBytes_CLHASH();
    }
    else{
      sum_bytes_cached_comps_ -= ccomp.SizeInBytes();
    }
    sum_size_cached_comps_ -= ccomp.nVars();
    num_cached_comps_--;
    sum_bytes_pure_cached_comp_data_ -= ccomp.data_only_byte_size();

    sys_overhead_sum_bytes_cached_comps_ -= ccomp.sys_overhead_SizeInBytes();
  }

  void incorporate_cache_hit(CacheableComponent &ccomp){
      num_cache_hits_++;
      sum_cache_hit_sizes_ += ccomp.nVars();
  }
  uint64_t cache_MB_memory_usage() {
      return cache_bytes_memory_usage() / 1000000;
  }
  mpz_class final_solution_count_ = 0;

  double implicitBCP_miss_rate() {
      if(num_failed_lit_tests_ == 0) return 0.0;
      return (num_failed_lit_tests_ - num_failed_literals_detected_) / (double) num_failed_lit_tests_;
  }
  uint64_t num_clauses() {
    return num_long_clauses_ + num_binary_clauses_ + num_unit_clauses_;
  }
  uint64_t num_conflict_clauses() {
    return num_long_conflict_clauses_ + num_binary_conflict_clauses_;
  }

  uint64_t clause_deletion_interval() {
    return 10000 + 10 * times_conflict_clauses_cleaned_;
  }

  void set_final_solution_count_projected(const mpz_class &count);
  void set_final_solution_count(const mpz_class &count);

  const mpz_class &final_solution_count() const {
    return final_solution_count_;
  }

  void incorporateConflictClauseData(const vector<Lit> &clause) {
    if (clause.size() == 1)
      num_unit_clauses_++;
    else if (clause.size() == 2)
      num_binary_conflict_clauses_++;
    num_long_conflict_clauses_++;
  }
  void incorporateClauseData(const vector<Lit> &clause) {
    if (clause.size() == 1)
      num_unit_clauses_++;
    else if (clause.size() == 2)
      num_binary_clauses_++;
    else
      num_long_clauses_++;
  }

  std::string getFinalSolutionCountStr();
  void printShort();
  void printShortFormulaInfo() {
    cout << "c variables (all/used/free): \t" << endl;

    cout << "c clauses (all/long/binary/unit): ";
    cout << num_clauses() << "/" << num_long_clauses_;
    cout << "/" << num_binary_clauses_ << "/" << num_unit_clauses_ << endl;
  }
  uint32_t getNumDecisions() {
    return num_decisions_;
  }

  double avgCachedSize() {
    if (num_cache_hits_ == 0)
      return 0.0;
    return (double) sum_size_cached_comps_
        / (double) num_cached_comps_;
  }

  double avgCacheHitSize() {
    if (num_cache_hits_ == 0)
      return 0.0;
    return (double) sum_cache_hit_sizes_ / (double) num_cache_hits_;
  }

  long double getAvgComponentSize() {
    return sum_size_cached_comps_ / (long double) num_cached_comps_;
  }

  uint64_t cached_comp_count() {
    return num_cached_comps_;
  }

  uint64_t cache_hits() {
    return num_cache_hits_;
  }

  double cache_miss_rate() {
    if(num_cache_look_ups_ == 0) return 0.0;
    return (num_cache_look_ups_ - num_cache_hits_)
        / (double) num_cache_look_ups_;
  }

  long double getAvgCacheHitSize() {
    if(num_cache_hits_ == 0) return 0.0;
    return sum_cache_hit_sizes_ / (long double) num_cache_hits_;
  }
};

#endif /* STATISTICS_H_ */
