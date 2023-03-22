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
  unsigned long num_original_clauses_ = 0;

  //This variable is used to periodically increase cache scores
  unsigned long numcachedec_ = 0;

  // number of variables remaining
  unsigned long num_variables_ = 0;
  // number of variables that actually occurs in clauses
  unsigned long num_used_variables_ = 0;
  unsigned long num_free_variables_ = 0;
  unsigned long num_free_projected_variables_ = 0;

  /// different clause counts

  // number of clauses after preprocessing
  unsigned long num_long_clauses_ = 0;
  unsigned long num_binary_clauses_ = 0;

  unsigned long num_long_conflict_clauses_ = 0;
  unsigned long num_binary_conflict_clauses_ = 0;

  unsigned long times_conflict_clauses_cleaned_ = 0;

  unsigned long num_unit_clauses_ = 0;
  /// number of all decisions made
  unsigned long num_decisions_ = 0;
  /// number of all implications derived
  unsigned long num_implications_ = 0;
  // number of all failed literal detections
  unsigned long num_failed_literals_detected_ = 0;
  unsigned long num_failed_literal_tests_ = 0;
  // number of all conflicts occurred
  unsigned long num_conflicts_ = 0;

  // number of clauses overall learned
  unsigned num_clauses_learned_ = 0;

  unsigned last_restart_decisions = 0;
  unsigned num_restarts = 0;
  unsigned next_restart = 1000;
  unsigned next_restart_diff = 1000;

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
    sum_size_cached_comps_ += ccomp.num_variables();
    num_cached_comps_++;
    total_num_cached_comps_++;
    overall_bytes_comps_stored_ += ccomp.SizeInBytes();
    overall_num_cache_stores_ += ccomp.num_variables();
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
    sum_size_cached_comps_ -= ccomp.num_variables();
    num_cached_comps_--;
    sum_bytes_pure_cached_comp_data_ -= ccomp.data_only_byte_size();

    sys_overhead_sum_bytes_cached_comps_ -= ccomp.sys_overhead_SizeInBytes();
  }

  void incorporate_cache_hit(CacheableComponent &ccomp){
      num_cache_hits_++;
      sum_cache_hit_sizes_ += ccomp.num_variables();
  }
  unsigned long cache_MB_memory_usage() {
      return cache_bytes_memory_usage() / 1000000;
  }
  mpz_class final_solution_count_ = 0;

  double implicitBCP_miss_rate() {
      if(num_failed_literal_tests_ == 0) return 0.0;
      return (num_failed_literal_tests_ - num_failed_literals_detected_) / (double) num_failed_literal_tests_;
  }
  unsigned long num_clauses() {
    return num_long_clauses_ + num_binary_clauses_ + num_unit_clauses_;
  }
  unsigned long num_conflict_clauses() {
    return num_long_conflict_clauses_ + num_binary_conflict_clauses_;
  }

  unsigned long clause_deletion_interval() {
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
    cout << "c variables (all/used/free): \t";
    cout << num_variables_ << "/" << num_used_variables_ << "/";
    cout << num_variables_ - num_used_variables_ << endl;

    cout << "c clauses (all/long/binary/unit): ";
    cout << num_clauses() << "/" << num_long_clauses_;
    cout << "/" << num_binary_clauses_ << "/" << num_unit_clauses_ << endl;
  }
  unsigned getNumDecisions() {
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

  unsigned long cached_comp_count() {
    return num_cached_comps_;
  }

  unsigned long cache_hits() {
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
