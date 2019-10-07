/*
 * basic_types.h
 *
 *  Created on: Jun 24, 2012
 *      Author: Marc Thurley
 */

#ifndef SOLVER_CONFIG_H_
#define SOLVER_CONFIG_H_

enum polarity_type {
    polar_false, polar_true, polar_default, polar_cached, reverse_cached, random_pol, new_random_pol
};

struct SolverConfiguration {

  bool perform_non_chron_back_track = true;

  // TODO component caching cannot be deactivated for now!
  bool perform_component_caching = true;
  bool perform_failed_lit_test = true;
  bool perform_pre_processing = true;
  bool allowrestart = false;
  bool allowactivitydecrease = false;
  bool tworandom = false;
  bool allowrandomness = false;
  bool usegraphdecom = false;
  bool usepcc = false;
  // float activityparam = 1;
  unsigned long randomseed = 1000;
  unsigned long restartafterdecisions = 1000;
  unsigned long decc = 0;
  unsigned long maxdec = 5000000;
  unsigned long minconflicts_ = 500;
  unsigned long cleardbafter = 10000;
  polarity_type polarity_config = polar_default;

  unsigned long time_bound_seconds = 100000;

  bool verbose = false;

  bool maxdecterminate = false;

  // quiet = true will override verbose;
  bool quiet = false;

  bool qbf = false;

  bool usecachetencoding = false;

  bool useposet = false;
  //independent support file 
  bool useindependentsupport = false;
  bool projected = false;
  // file for independent support
  string independent_support_file;

  // string posetfile;

  bool useIsomorphicComponentCaching = false;

  unsigned int isomorphicComponentCachingn = 5000; 
};

#endif /* SOLVER_CONFIG_H_ */
