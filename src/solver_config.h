/*
 * basic_types.h
 *
 *  Created on: Jun 24, 2012
 *      Author: Marc Thurley
 */

#ifndef SOLVER_CONFIG_H_
#define SOLVER_CONFIG_H_

struct SolverConfiguration {
  // TODO comp caching cannot be deactivated for now!
  bool perform_non_chron_back_track = true;
  bool perform_comp_caching = true;
  bool perform_failed_lit_test = true;
  bool perform_pre_processing = true;
  bool perform_pcc = true; // probabilistic comp caching
  bool use_csvsads = true; // polarity heuristic
  int verb = 1;
  int restart = 1;

  uint64_t time_bound_seconds = 100000;
  uint64_t randomseed = 1000;
  unsigned hashrange = 1;
  unsigned maxdec = 5000000;
  unsigned minconflicts_ = 500;
  float delta = 0.05;
};

#endif /* SOLVER_CONFIG_H_ */
