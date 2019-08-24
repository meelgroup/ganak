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
  mp_exp_t exp;
  std::string significand = final_solution_count_.get_str(exp);
  cout <<  significand << "x10^" << exp-(int)significand.size();
}

void DataAndStatistics::writeToFile(const string & file_name) {
  ofstream out(file_name, ios_base::app);
  unsigned pos = input_file_.find_last_of("/\\");
  out << "<tr>" << endl;
  out << "<td>" << input_file_.substr(pos + 1) << "</td>" << endl;
  out << "<td>" << num_original_variables_ << "</td>" << endl;
  out << "<td>" << num_original_clauses_ << "</td>" << endl;
  out << "<td>" << num_decisions_ << "</td>" << endl;
  out << "<td>" << time_elapsed_ << "</td>" << endl;

  mp_exp_t exp;
  string s = final_solution_count_.get_str(exp);
  if (final_solution_count_ == 0)
    s = "UNSAT";
  out << "<td>" << s << "</td>" << endl;
  out << "</tr>" << endl;
}

void DataAndStatistics::printShort() {
  if (exit_state_ == TIMEOUT) {
    cout << endl << " TIMEOUT !" << endl;
    return;
  }
  cout << endl << endl;
  cout << "variables (total / active / free)\t" << num_variables_ << "/"
      << num_used_variables_ << "/" << num_variables_ - num_used_variables_
      << endl;
  cout << "clauses (removed) \t\t\t" << num_original_clauses_ << " ("
      << num_original_clauses_ - num_clauses() << ")" << endl;
  cout << "decisions \t\t\t\t" << num_decisions_ << endl;
  cout << "Max. decision level \t\t\t\t" << max_decision_level_ << endl;
  cout << "conflicts \t\t\t\t" << num_conflicts_ << endl;
  cout << "conflict clauses (all/bin/unit) \t";
  cout << num_conflict_clauses();
  cout << "/" << num_binary_conflict_clauses_ << "/" << num_unit_clauses_
      << endl << endl;
  cout << "failed literals found by implicit BCP \t "
      << num_failed_literals_detected_ << endl;


  cout << "implicit BCP miss rate \t " << implicitBCP_miss_rate() * 100 << "%";
  cout << endl;
  cout << "bytes cache size     \t" << cache_bytes_memory_usage()  << "\t"
      << endl;
  cout << "bytes cache (overall) \t" << overall_cache_bytes_memory_stored()
      << "" << endl;
  cout << "bytes cache components (overall) \t"  << overall_bytes_components_stored_ << ""<< endl;
  cout << "bytes cache (infra / comps) "
      << (cache_infrastructure_bytes_memory_usage_) << "/"
      << sum_bytes_cached_components_  << "\t" << endl;

  cout << "bytes pure comp data (curr)    " << sum_bytes_pure_cached_component_data_  << "" << endl;
  cout << "bytes pure comp data (overall) " <<overall_bytes_pure_stored_component_data_ << "" << endl;

  cout << "bytes cache with sysoverh (curr)    " << sys_overhead_sum_bytes_cached_components_  << "" << endl;
  cout << "bytes cache with sysoverh (overall) " << sys_overhead_overall_bytes_components_stored_ << "" << endl;


  cout << "cache (lookup / stores / hits) \t\t\t" << num_cache_look_ups_ << "/"
      << total_num_cached_components_ << "/"
      << num_cache_hits_ << endl;
  cout << "cache miss rate " << cache_miss_rate() * 100 << "%" << endl;
  cout << " avg. variable count (stores / hits) \t" << getAvgComponentSize()
      << "/" << getAvgCacheHitSize() << endl << endl;
  cout << "\n# solutions " << endl;

  print_final_solution_count();

  cout << "\n# END" << endl << endl;
  cout << "time: " << time_elapsed_ << "s\n\n";
}
