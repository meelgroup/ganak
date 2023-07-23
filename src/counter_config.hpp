/******************************************
Copyright (C) 2023 Authors of GANAK, see AUTHORS file

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

#pragma once

#include <cstdint>
#include <string>
#include <cassert>
#include <iostream>
#include "common.hpp"
using std::string;

enum class branch_t {old_ganak, sharptd, gpmc};

struct CounterConfiguration {
  double act_exp = 0.99; // 1.0 is VERY slow. Seems good, also tried 0.985
  bool do_pre_processing = true;
  int verb = 1;
  int do_restart = 0;
  uint64_t first_restart = 1000U;
  uint64_t next_restart = 10000U;
  double restart_cutoff_mult = 0.8;
  uint64_t maximum_cache_size_MB = 0;

  // 7 == conflict/luby-based, taking into account first_restart
  // 6 == conflict/static, taking into account next_restart
  int restart_type = 7;

  double probe_only_after_ratio = 0.25;
  int alluip_inc_act = 1;
  int do_cache_score = 1;
  int do_single_bump = 1; // non-single bump is OLD ganak
  uint32_t rdb_cls_target = 10000;
  int rdb_keep_used = 1; // quite a bit faster on lower time cut-off
                         // but loses the edge after ~2000s

  // 3 = TRUE, 2 = FALSE, 1 = standard, 0 = cache
  int polar_type = 1;

  int do_vivify = 1;
  uint32_t vivif_every = 60000;
  double vivif_mult = 1.0;
  int do_extra_cl_bump = 1;
  int do_buddy = 0;
  uint32_t buddy_max_cls = 10;

  uint32_t td_varlim = 150000;
  double td_denselim = 0.10;
  double td_ratiolim = 30.0;
  int td_with_red_bins = 1; // a MUST or we are VERY slow
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
