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

void DataAndStatistics::printShort(const Counter* solver, const ComponentCache* cache_) const {
  solver->print_restart_data();
  cout << "c cls irred                      " << num_irred_clauses() << endl;
  cout << "c decisions                      " << num_decisions_ << endl;
  cout << "c conflicts                      " << num_conflicts_ << endl;
  cout << "c conflict cls (long/bin/u)      " << std::fixed
    << solver->get_num_irred_long_cls() << "/"
    << num_binary_red_clauses_ << "/" << num_unit_red_clauses_ << endl;
  cout << "c conflict cls compacted         " << solver->num_conflict_clauses_compacted() << endl;
  cout << "c looks/look-computes            "
    << lookaheads << "/" << lookahead_computes << endl;
  cout << "c flits/bplits/tests K           "
    << std::left
    << std::setw(6) << (num_failed_literals_detected_/1000ULL)
    << std::setw(6) << (num_failed_bprop_literals_failed/1000ULL)
    << std::setw(6) << (num_failed_lit_tests_/1000ULL)
    << " -- " << std::setprecision(2) << safe_div( num_failed_literals_detected_, num_failed_lit_tests_)
    << " -- " << std::setprecision(2) << safe_div( num_failed_literals_detected_+num_failed_bprop_literals_failed, num_failed_lit_tests_)  << endl;

  cout << "c implicit BCP miss rate         "
    << std::setprecision(2) << implicitBCP_miss_rate() * 100 << "%";
  cout << endl;
  cout << "c cache entries K                " << (cache_->get_num_entries_used()/1000ULL) << endl;
  /* cout << "c MB cache size                  " */
  /*   << std::setprecision(3) << in_MB(cache_bytes_memory_usage()) << "\t" << endl; */
  cout << "c MB cache (overall)             "
    << std::setprecision(3) << in_MB(overall_cache_bytes_memory_stored()) << "" << endl;
  cout << "c cache K (lookup/ stores/ hits) "
    << std::left
    << std::setw(6) << (num_cache_look_ups_/(1000ULL)) << " "
    << std::setw(6) << (total_num_cached_comps_ /(1000ULL)) << " "
    << std::setw(6) << (num_cache_hits_ /(1000ULL)) << " " << endl;

  cout << "c cache pollutions call/removed  "
    << cache_pollutions_called << "/"
    << cache_pollutions_removed << endl;
  cout << "c cache miss rate                "
    << std::setprecision(3) << cache_miss_rate() << endl;
  cout << "c avg. var count (stores / hits) "
    << getAvgComponentSize()
    << "/" << getAvgCacheHitSize() << endl;
  cout << "c " << endl;
}

