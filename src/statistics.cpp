/*
 * statistics.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: mthurley
 */

#include "statistics.h"

#include <iostream>
#include <fstream>

using namespace std;

void DataAndStatistics::print_final_solution_count() {
  cout << final_solution_count_.get_str();
}

void DataAndStatistics::writeToFile(const string & file_name, const bool pmc) {
  ofstream out(file_name, ios_base::app);
  if (exit_state_ == TIMEOUT) {
    out << endl << "c TIMEOUT !" << endl;
    return;
  }
  out << "c " << endl;
  out << "c " << endl;
  out << "c variables (total / active / free)\t" << num_variables_ << "/"
      << num_used_variables_ << "/" << num_variables_ - num_used_variables_
      << endl;
  out << "c clauses (removed) \t\t\t" << num_original_clauses_ << " ("
      << num_original_clauses_ - num_clauses() << ")" << endl;
  out << "c decisions \t\t\t\t" << num_decisions_ << endl;
  out << "c Max. decision level \t\t\t" << max_decision_level_ << endl;
  out << "c conflicts \t\t\t\t" << num_conflicts_ << endl;
  out << "c conflict clauses (all/bin/unit) \t";
  out << num_conflict_clauses();
  out << "/" << num_binary_conflict_clauses_ << "/" << num_unit_clauses_
      << endl;
  out << "c " << endl;
  out << "c failed literals found by implicit BCP \t "
      << num_failed_literals_detected_ << endl;


  out << "c implicit BCP miss rate \t " << implicitBCP_miss_rate() * 100 << "%";
  out << endl;
  out << "c bytes cache size     \t" << cache_bytes_memory_usage()  << "\t"
      << endl;
  out << "c bytes cache (overall) \t" << overall_cache_bytes_memory_stored()
      << "" << endl;
  out << "c bytes cache components (overall) \t"  << overall_bytes_components_stored_ << ""<< endl;
  out << "c bytes cache (infra / comps) "
      << (cache_infrastructure_bytes_memory_usage_) << "/"
      << sum_bytes_cached_components_  << "\t" << endl;

  out << "c bytes pure comp data (curr)    " << sum_bytes_pure_cached_component_data_  << "" << endl;
  out << "c bytes pure comp data (overall) " <<overall_bytes_pure_stored_component_data_ << "" << endl;
  out << "c bytes cache with sysoverh (curr)    " << sys_overhead_sum_bytes_cached_components_  << "" << endl;
  out << "c bytes cache with sysoverh (overall) " << sys_overhead_overall_bytes_components_stored_ << "" << endl;
  out << "c cache (lookup / stores / hits) \t\t\t" << num_cache_look_ups_ << "/"
      << total_num_cached_components_ << "/"
      << num_cache_hits_ << endl;
  out << "c cache miss rate " << cache_miss_rate() * 100 << "%" << endl;
  out << "c avg. variable count (stores / hits) \t" << getAvgComponentSize()
      << "/" << getAvgCacheHitSize() << endl;
  out << "c " << endl;
  out << "c " << endl;
  out << "c # solutions " << endl;
  if (pmc) {
    out << "s pmc " << std::flush;
  } else {
    out << "s mc " << std::flush;
  }
  out << final_solution_count_.get_str();
  out << endl;
  out << "c # END" << endl;
  out << "c " << endl;
  out << "c time: " << time_elapsed_ << "s" << endl;
  
}

void DataAndStatistics::printShort(const bool pmc) {
  if (exit_state_ == TIMEOUT) {
    cout << endl << "c TIMEOUT !" << endl;
    return;
  }
  cout << "c " << endl;
  cout << "c " << endl;
  cout << "c variables (total / active / free)\t" << num_variables_ << "/"
      << num_used_variables_ << "/" << num_variables_ - num_used_variables_
      << endl;
  cout << "c clauses (removed) \t\t\t" << num_original_clauses_ << " ("
      << num_original_clauses_ - num_clauses() << ")" << endl;
  cout << "c decisions \t\t\t\t" << num_decisions_ << endl;
  cout << "c Max. decision level \t\t\t" << max_decision_level_ << endl;
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
  cout << "c bytes cache components (overall) \t"  << overall_bytes_components_stored_ << ""<< endl;
  cout << "c bytes cache (infra / comps) "
      << (cache_infrastructure_bytes_memory_usage_) << "/"
      << sum_bytes_cached_components_  << "\t" << endl;

  cout << "c bytes pure comp data (curr)    " << sum_bytes_pure_cached_component_data_  << "" << endl;
  cout << "c bytes pure comp data (overall) " <<overall_bytes_pure_stored_component_data_ << "" << endl;
  cout << "c bytes cache with sysoverh (curr)    " << sys_overhead_sum_bytes_cached_components_  << "" << endl;
  cout << "c bytes cache with sysoverh (overall) " << sys_overhead_overall_bytes_components_stored_ << "" << endl;
  cout << "c cache (lookup / stores / hits) \t\t\t" << num_cache_look_ups_ << "/"
      << total_num_cached_components_ << "/"
      << num_cache_hits_ << endl;
  cout << "c cache miss rate " << cache_miss_rate() * 100 << "%" << endl;
  cout << "c avg. variable count (stores / hits) \t" << getAvgComponentSize()
      << "/" << getAvgCacheHitSize() << endl;
  cout << "c " << endl;
  cout << "c " << endl;
  cout << "c # solutions " << endl;
  if (pmc) {
    cout << "s pmc " << std::flush;
  } else {
    cout << "s mc " << std::flush;
  }
  print_final_solution_count();
  cout << endl;
  cout << "c # END" << endl;
  cout << "c " << endl;
  cout << "c time: " << time_elapsed_ << "s" << endl;
}
