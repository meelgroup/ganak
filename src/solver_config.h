/*
 * basic_types.h
 *
 *  Created on: Jun 24, 2012
 *      Author: Marc Thurley
 */

#pragma once

#include <cstdint>
#include <string>
#include <cassert>
#include <iostream>
#include "common.h"
using std::string;

enum class branch_t {old_ganak, sharptd, gpmc};

struct CounterConfiguration {
  // TODO comp caching cannot be deactivated for now!
  bool do_non_chron_back_track = true;
  double act_exp = 0.99;
  bool do_pre_processing = true;
  int verb = 1;
  int do_restart = 0;
  uint64_t first_restart = 1000U;
  uint64_t next_restart = 100000000U;
  double restart_cutoff_mult = 0.8;
  uint64_t maximum_cache_size_MB = 0;
  int restart_type = 5;
  double lookahead_depth = 1.2;
  uint32_t lookahead_num = 2;
  int failed_lit_probe_type = 2;
  double probe_only_after_ratio = 0.25;
  double num_probe_multi = 1.0; //TODO should try 2.0
  int alluip_inc_act = 1;
  int bprop = 1;
  int do_on_path_print = 0;
  int do_lookahead = 0;
  int do_cache_score = 1;
  int do_single_bump = 1; // non-single bump is OLD ganak
  uint32_t rdb_cls_target = 10000;
  int rdb_keep_used = 1;
  int polar_type = 1; // 3 = TRUE, 2 = FALSE, 1 = standard, 0 = cache

  uint32_t td_varlim = 150000;
  double td_denselim = 0.10;
  double td_ratiolim = 30.0;
  double tw_vare_lim = 1.1;
  double tw_coef_tdscore = 100.0;
  branch_t branch_type = branch_t::sharptd;
  branch_t branch_fallback_type = branch_t::old_ganak;

  uint64_t seed = 0;
  double delta = 0.05;

  string get_branch_type_str() const {
    if (branch_type == branch_t::gpmc) return "gpmc";
    else if (branch_type == branch_t::sharptd) return "sharptd";
    else if (branch_type == branch_t::old_ganak) return "ganak";
    else release_assert(false);
  }
};

inline std::string branch_type_to_str(const branch_t& t) {
  if (t == branch_t::gpmc) return "gpmc";
  else if (t == branch_t::sharptd) return "sharptd";
  else if (t == branch_t::old_ganak) return "ganak";
  else {
    std::cout << "ERROR: Can't translate branch type" << std::endl;
    exit(-1);
  }
}
inline branch_t parse_branch_type(const std::string& name) {
  if (name == "gpmc") return branch_t::gpmc;
  else if (name == "sharptd") return branch_t::sharptd;
  else if (name == "ganak") return branch_t::old_ganak;
  else {
    std::cout << "ERROR: Wrong branch type: '" << name << "'" << std::endl;
    exit(-1);
  }
}
