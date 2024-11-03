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
using std::string;

struct CounterConfiguration {
  int verb = 1;
  int do_restart = 0;
  int do_chronobt = 1;
  uint64_t first_restart = 20000U;
  double restart_cutoff_mult = 0.8;
  uint64_t maximum_cache_size_MB = 0;
  int restart_type = 8;
  int do_readjust_for_restart = 1;
  int max_num_rst = -1;
  uint32_t lbd_cutoff_always_keep_cube = 3;
  uint32_t max_num_cubes_per_restart = 6;

  int cache_time_update = 2;

  int decide = 0; // 0 = TD, 1 = ignore TD
  uint32_t rdb_cls_target = 10000;
  int rdb_keep_used = 0; // quite a bit faster on lower time cut-off
                         // but loses the edge after ~2000s

  uint32_t reduce_db_everyN = 10000;
  uint32_t consolidate_every_n = 30000;
  int polar_type = 0;
  uint32_t base_lbd_cutoff = 2;
  bool update_lbd_cutoff = 0;

  int do_vivify = 1;
  uint32_t vivif_every = 60000;
  double vivif_mult = 1.0;
  uint32_t vivif_outer_every_n = 1;
  int do_buddy = 0;
  uint32_t buddy_max_cls = 6;
#ifdef CHECK_COUNT
  int do_cube_check_count = 1;
#else
  int do_cube_check_count = 0;
#endif
  int do_use_cache = 1;
  int do_init_activity_scores = 1;
  int vsads_readjust_every = 256;
  uint32_t tot_used_cutoff_vivif = 50;


  bool do_td = 1;
  uint32_t td_varlim = 150000;
  double td_ratiolim = 100.0;
  double td_maxweight = 60;
  double td_minweight = 7;
  double td_divider = 1e3;
  double do_td_weight = 1;
  double td_exp_mult = 1.1;
  int do_check_td_vs_ind = 1;
  int64_t td_steps = 1e5;
  int td_iters = 900;
  int td_lookahead = -1;
  int td_lookahead_tw_cutoff = 26;
  int td_lookahead_iters = 10;
  int td_look_only_weight = false;
  int td_max_edges = 70000;

  int do_use_sat_solver = 1;
  int sat_restart_mult = 300;
  int do_sat_restart = 1;
  int do_sat_polar_cache = 1;
  int do_sat_vsids = 1;

  uint64_t seed = 0;
  double delta = 0.2;

  double appmc_timeout = -1;
  double appmc_epsilon = 0.8;
};
