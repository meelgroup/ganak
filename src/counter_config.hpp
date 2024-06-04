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
  double act_exp = 0.99;
  bool do_pre_processing = true;
  int verb = 1;
  int do_restart = 0;
  uint64_t first_restart = 20000U;
  double restart_cutoff_mult = 0.8;
  uint64_t maximum_cache_size_MB = 0;
  double var_freq_divider = 20.0; //10 is best for vsads_readjust_every = 0

  int restart_type = 8;
  int do_readjust_for_restart = 1;

  int do_comp_sort = 0; // they are very similar, see: out-ganak-6910211
  double probe_only_after_ratio = 0.25;
  int alluip_inc_act = 1;
  int cache_time_update = 2;
  int do_cache_reverse_sort = 1;

  int decide = 0; // 0 = sstd, 1 = gpmc
  uint32_t rdb_cls_target = 10000;
  int rdb_keep_used = 1; // quite a bit faster on lower time cut-off
                         // but loses the edge after ~2000s

  uint32_t reduce_db_everyN = 10000;
  uint32_t consolidate_every_n = 30000;
  int polar_type = 0;
  uint32_t base_lbd_cutoff = 1;

  int do_vivify = 1;
  uint32_t vivif_every = 60000;
  double vivif_mult = 1.0;
  uint32_t vivif_outer_every_n = 1;
  int do_buddy = 0;
  uint32_t buddy_max_cls = 12; // see out-ganak-2060064.pbs101-7
#ifdef CHECK_COUNT
  int do_cube_check_count = 1;
#else
  int do_cube_check_count = 0;
#endif
  int do_use_cache = 1;
  int do_init_activity_scores = 1;
  int do_check_unkn_bin = 1;
  int vsads_readjust_every = 256;


  bool do_td = 1;
  uint32_t td_varlim = 150000;
  double td_ratiolim = 100.0;
  double td_maxweight = 10.0; //4.0 is best for vsads_readjust_every = 0
  double td_minweight = 0.05; //0.1 is best for vsads_readjust_every = 0
  double td_divider = 1e3;
  double do_td_weight = 1;
  double td_exp_mult = 1.2;
  int do_check_td_vs_ind = 1;
  int64_t td_steps = 1e5;
  int td_iters = 900;
  int td_lookahead = -1;
  int td_lookahead_tw_cutoff = 20;
  int td_lookahead_iters = 10;
  int td_look_only_weight = false;
  int td_max_edges = 70000;


  uint64_t seed = 0;
  double delta = 0.05;

  double appmc_timeout = -1;
  double appmc_epsilon = 0.2;
};
