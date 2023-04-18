/*
 * basic_types.h
 *
 *  Created on: Jun 24, 2012
 *      Author: Marc Thurley
 */

#pragma once

#include <cstdint>

struct CounterConfiguration {
  // TODO comp caching cannot be deactivated for now!
  bool do_non_chron_back_track = true;
  bool do_comp_caching = true;
  int failed_lit_probe_type = 2;
  double exp = 1.0;
  bool do_pre_processing = true;
  int verb = 1;
  int do_restart = 0;
  uint64_t first_restart = 1000U;
  uint64_t next_restart = 100000000U;
  uint64_t maximum_cache_size_bytes_ = 0;
  int restart_type = 5;
  double lookahead_depth = 1.2;
  uint32_t lookahead_num = 2;
  double probe_only_after_ratio = 0.25;
  double num_probe_multi = 1.0;
  int alluip_inc_act = 1;
  int bprop = 1;
  int do_on_path_print = 0;
  int do_lookahead = 0;
  int do_cache_score = 1;
  int do_single_bump = 1; // non-single bump is OLD ganak

  uint64_t seed = 0;
  uint32_t maxdec = 5000000;
  uint32_t minconflicts_ = 500;
  double delta = 0.05;
};
