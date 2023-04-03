/*
 * statistics.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: mthurley
 */

#include "statistics.h"
#include "comp_cache.h"
#include "solver.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <solver.h>

static double in_MB(uint64_t bytes) {
  return (double)bytes/(double)(1024*1024);
}

double safe_div(double a, double b) {
  if (b == 0) return 0;
  else return a/b;
}

void DataAndStatistics::printShort(const Solver* solver, const ComponentCache* cache_) const {
  cout << "c " << endl;
  cout << "c      --- FINISHED ---" << endl;
  cout << "c cls irred                      " << num_irred_clauses() << endl;
  cout << "c decisions                      " << num_decisions_ << endl;
  cout << "c conflicts                      " << num_conflicts_ << endl;
  cout << "c conflict cls (long/bin/u)      "
    << solver->get_num_irred_long_cls() << "/"
    << num_binary_red_clauses_ << "/" << num_unit_red_clauses_ << endl;
  cout << "c conflict cls compacted         " << solver->num_conflict_clauses_compacted() << endl;
  cout << "c failed lits by iBCP/tests      " << num_failed_literals_detected_ << "/" << num_failed_lit_tests_ << " -- " << safe_div((double)num_failed_literals_detected_, (double)num_failed_lit_tests_)  << endl;

  cout << "c implicit BCP miss rate         "
    << std::setprecision(2) << implicitBCP_miss_rate() * 100 << "%";
  cout << endl;
  cout << "c cache entries                  " << cache_->get_num_entries_used() << endl;
  cout << "c MB cache size                  "
    << std::setprecision(3) << in_MB(cache_bytes_memory_usage()) << "\t" << endl;
  cout << "c MB cache (overall)             "
    << std::setprecision(3) << in_MB(overall_cache_bytes_memory_stored()) << "" << endl;
  cout << "c MB cache comps (ovrall)        "
    << std::setprecision(3) << in_MB(overall_bytes_comps_stored_) << endl;
  cout << "c MB cache (infra / comps)       "
    << std::setprecision(3) << in_MB(cache_infrastructure_bytes_memory_usage_) << "/"
    << std::setprecision(3) << in_MB(sum_bytes_cached_comps_) << endl;
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
  cout << "c " << endl;
}

