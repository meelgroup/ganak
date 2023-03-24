/*
 * basic_types.h
 *
 *  Created on: Jun 24, 2012
 *      Author: Marc Thurley
 */

#pragma once

#include <cstdint>

struct SolverConfiguration {
  // TODO comp caching cannot be deactivated for now!
  bool perform_non_chron_back_track = true;
  bool perform_comp_caching = true;
  bool perform_failed_lit_probe = true;
  bool perform_pre_processing = true;
  bool perform_pcc = true; // probabilistic comp caching
  bool use_csvsads = true; // polarity heuristic
  int verb = 1;
  int restart = 1;

  uint64_t time_bound_seconds = 100000;
  uint64_t randomseed = 1000;
  uint32_t hashrange = 1;
  uint32_t maxdec = 5000000;
  uint32_t minconflicts_ = 500;
  float delta = 0.05;
};
