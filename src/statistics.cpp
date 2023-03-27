/*
 * statistics.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: mthurley
 */

#include "statistics.h"
#include "solver.h"

#include <iostream>
#include <iomanip>
#include <fstream>

std::string DataAndStatistics::getFinalSolutionCountStr() const {
    return final_solution_count_.get_str();
}

void DataAndStatistics::set_final_solution_count_projected(const mpz_class &count) {
  mpz_mul_2exp(
    final_solution_count_.get_mpz_t (),
    count.get_mpz_t (),
    inst->get_must_mult_exp2());
}

void DataAndStatistics::set_final_solution_count(const mpz_class &count) {
  // set final_solution_count_ = count * 2^(nVars_ - num_used_variables_)
  mpz_mul_2exp(
    final_solution_count_.get_mpz_t (),
    count.get_mpz_t (),
    nVars_ - inst->get_must_mult_exp2());
}

static double in_MB(uint64_t bytes) {
  return (double)bytes/(double)(1024*1024);
}

void DataAndStatistics::printShort() const {
  cout << "c " << endl;
  cout << "c  -- FINISHED ---" << endl;
  cout << "c vars (total / active / free)   " << nVars_ << "/" << endl;
  cout << "c cls (removed)                  "
    << num_original_clauses_ << " (" << num_original_clauses_ - num_clauses() << ")" << endl;
  cout << "c decisions                      " << num_decisions_ << endl;
  cout << "c conflicts                      " << num_conflicts_ << endl;
  cout << "c conflict cls (all/bin/unit)    " << num_conflict_clauses();
  cout << "/" << num_binary_conflict_clauses_ << "/" << num_unit_clauses_ << endl;
  cout << "c failed lits found by iBCP      " << num_failed_literals_detected_ << endl;

  cout << "c implicit BCP miss rate         " << implicitBCP_miss_rate() * 100 << "%";
  cout << endl;
  cout << "c MB cache size                  "
    << std::setprecision(3) << in_MB(cache_bytes_memory_usage()) << "\t" << endl;
  cout << "c MB cache (overall)             "
    << std::setprecision(3) << in_MB(overall_cache_bytes_memory_stored()) << "" << endl;
  cout << "c MB cache comps (ovrall)        "
    << std::setprecision(3) << in_MB(overall_bytes_comps_stored_) << endl;
  cout << "c MB cache (infra / comps)       "
    << std::setprecision(3) << in_MB(cache_infrastructure_bytes_memory_usage_) << "/"
    << std::setprecision(3) << in_MB(sum_bytes_cached_comps_) << endl;
  cout << "c MB pure comp dat (curr)        "
    << std::setprecision(3) << in_MB(sum_bytes_pure_cached_comp_data_) << endl;
  cout << "c MB pure comp dat (ovrall)      "
    << std::setprecision(3) << in_MB(overall_bytes_pure_stored_comp_data_) << endl;
  cout << "c MB cache w/ sysoverh (curr)    "
    << std::setprecision(3) << in_MB(sys_overhead_sum_bytes_cached_comps_) << endl;
  cout << "c MB cache w/ sysoverh (ovrall)  "
    << std::setprecision(3) << in_MB(sys_overhead_overall_bytes_comps_stored_) << endl;
  cout << "c cache (lookup / stores / hits) "
    << num_cache_look_ups_ << "/"
    << total_num_cached_comps_ << "/"
    << num_cache_hits_ << endl;
  cout << "c cache miss rate                "
    << std::setprecision(3) << cache_miss_rate() * 100 << "%" << endl;
  cout << "c avg. variable count (stores / hits) \t"
    << getAvgComponentSize()
    << "/" << getAvgCacheHitSize() << endl;
  cout << "c time: " << time_elapsed_ << "s" << endl;
  cout << "c " << endl;
}

void DataAndStatistics::print_solution() const {
  if (final_solution_count_ == 0) {
    cout << "s UNSATISFIABLE " << endl;
  } else {
    cout << "s SATISFIABLE " << endl;
  }
  if (inst->get_indep_support_given()) {
    cout << "s pmc ";
  } else {
    cout << "s mc ";
  }
  cout << getFinalSolutionCountStr() << endl;
}
