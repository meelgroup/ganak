/*
 * statistics.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: mthurley
 */

#include "statistics.h"
#include "solver.h"

#include <iostream>
#include <fstream>

std::string DataAndStatistics::getFinalSolutionCountStr() {
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

void DataAndStatistics::printShort() {
  cout << "c " << endl;
  cout << "c variables (total / active / free)\t" << nVars_ << "/" << endl;
  cout << "c clauses (removed) \t\t\t" << num_original_clauses_ << " ("
      << num_original_clauses_ - num_clauses() << ")" << endl;
  cout << "c decisions \t\t\t\t" << num_decisions_ << endl;
  cout << "c conflicts \t\t\t\t" << num_conflicts_ << endl;
  cout << "c conflict clauses (all/bin/unit) \t";
  cout << num_conflict_clauses();
  cout << "/" << num_binary_conflict_clauses_ << "/" << num_unit_clauses_
      << endl;
  cout << "c " << endl;
  cout << "c failed literals found by implicit BCP \t "
      << num_failed_literals_detected_ << endl;


  cout << "c implicit BCP miss rate \t " << implicitBCP_miss_rate() * 100 << "%";
  cout << endl;
  cout << "c bytes cache size     \t" << cache_bytes_memory_usage()  << "\t"
      << endl;
  cout << "c bytes cache (overall) \t" << overall_cache_bytes_memory_stored()
      << "" << endl;
  cout << "c bytes cache comps (overall) \t"  << overall_bytes_comps_stored_ << ""<< endl;
  cout << "c bytes cache (infra / comps) "
      << (cache_infrastructure_bytes_memory_usage_) << "/"
      << sum_bytes_cached_comps_  << "\t" << endl;

  cout << "c bytes pure comp data (curr)    " << sum_bytes_pure_cached_comp_data_  << "" << endl;
  cout << "c bytes pure comp data (overall) " <<overall_bytes_pure_stored_comp_data_ << "" << endl;
  cout << "c bytes cache with sysoverh (curr)    " << sys_overhead_sum_bytes_cached_comps_  << "" << endl;
  cout << "c bytes cache with sysoverh (overall) " << sys_overhead_overall_bytes_comps_stored_ << "" << endl;
  cout << "c cache (lookup / stores / hits) \t\t\t" << num_cache_look_ups_ << "/"
      << total_num_cached_comps_ << "/"
      << num_cache_hits_ << endl;
  cout << "c cache miss rate " << cache_miss_rate() * 100 << "%" << endl;
  cout << "c avg. variable count (stores / hits) \t" << getAvgComponentSize()
      << "/" << getAvgCacheHitSize() << endl;
  cout << "c " << endl;
  cout << "c " << endl;
  if (final_solution_count_ == 0) {
    cout << "s UNSATISFIABLE " << endl;
  } else {
    cout << "s SATISFIABLE " << endl;
  }
  cout << "c # solutions " << endl;
  if (inst->get_indep_support_given()) {
    cout << "s pmc ";
  } else {
    cout << "s mc ";
  }
  cout << getFinalSolutionCountStr() << endl;

  cout << endl;
  cout << "c # END" << endl;
  cout << "c " << endl;
  cout << "c time: " << time_elapsed_ << "s" << endl;
}
