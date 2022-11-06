/*
 * basic_types.h
 *
 *  Created on: Jun 24, 2012
 *      Author: Marc Thurley
 */

#ifndef SOLVER_CONFIG_H_
#define SOLVER_CONFIG_H_

enum polarity_type {
    polar_default, polaritycache
};

struct SolverConfiguration {
  // TODO component caching cannot be deactivated for now!
  bool perform_non_chron_back_track = true;
  bool perform_component_caching = true;
  bool perform_failed_lit_test = true;
  bool perform_pre_processing = true;
  bool perform_pcc = true; // probabilistic component caching
  bool use_csvsads = true; // polarity heuristic
  bool verbose = false;
  // quiet = true will override verbose;
  bool quiet = false;

  unsigned long time_bound_seconds = 100000;
  unsigned long randomseed = 1000;

  unsigned hashrange = 1;
  unsigned maxdec = 5000000;
  unsigned minconflicts_ = 500;

  float delta = 0.05;
  float csvsads_param = 0.9;

  polarity_type polarity_config = polaritycache;
};

#endif /* SOLVER_CONFIG_H_ */
