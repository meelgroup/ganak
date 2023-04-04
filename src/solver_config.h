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
  bool do_non_chron_back_track = true;
  bool do_comp_caching = true;
  int failed_lit_probe_type = 2;
  double exp = 0.5;
  bool do_pre_processing = true;
  bool do_pcc = true; // probabilistic comp caching
  int verb = 1;
  int do_restart = 1;
  uint64_t first_restart = 100000000U;
  uint64_t next_restart = 100000000U;
  uint64_t maximum_cache_size_bytes_ = 0;

  uint64_t seed = 0;
  uint32_t hashrange = 1;
  uint32_t maxdec = 5000000;
  uint32_t minconflicts_ = 500;
  double delta = 0.05;
};
