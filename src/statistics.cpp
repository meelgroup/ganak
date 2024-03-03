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

#include "statistics.hpp"
#include "comp_cache.hpp"
#include "counter.hpp"
#include "time_mem.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>

static double in_mb(uint64_t bytes) {
  return (double)bytes/(double)(1024*1024);
}

static double safe_div(double a, double b) {
  if (b == 0) return 0;
  else return a/b;
}

void DataAndStatistics::print_short(const Counter* counter, const CompCache* cache) const {
  counter->print_restart_data();
  verb_print(1, "cls long irred                 " << counter->get_num_irred_long_cls());
  verb_print(1, "decisions K                    "
    << std::left << std::setw(9) << decisions/1000
    << std::setw(16) << " -- Kdec/s: "
    << std::setprecision(2) << std::setw(9) << std::left << std::fixed
    << safe_div(decisions,(1000.0*(cpuTime()-counter->get_start_time())))
  );
  verb_print(1, "conflicts                      "
    << std::left << std::setw(9) << conflicts
    << "   " << std::left << std::setw(9)
    << std::setw(16) << " -- confl/s: "
    << std::setprecision(2) << std::setw(9) << std::left
    << safe_div(conflicts,((cpuTime()-counter->get_start_time())))
  );
  verb_print(1, "conflict cls (long/bin)      " << std::fixed
    << counter->get_num_long_reds() << "/" << num_binary_red_clauses_);

  /* verb_print(1, "lits /rem lits ccmin           " */
  /*   << std::setw(9) << orig_uip_lits/uip_cls << " " */
  /*   << std::setw(9) << ccmin_uip_lits/uip_cls << " " */
  /* ); */

  verb_print(1, "rem lits triedK/rem lits remK  "
    << std::setw(9) << rem_lits_tried/1000 << " "
    << std::setw(9) << rem_lits_with_bins/1000 << " "
  );

  verb_print(1, "avg clsz/rem lits avg/finalavg "
    << std::setw(9) << safe_div(orig_uip_lits,uip_cls) << " "
    << std::setw(9) << safe_div(uip_lits_ccmin,uip_cls) << " "
    << std::setw(9) << safe_div(final_cl_sz, uip_cls)
  );

  verb_print(1, "rdbs/low lbd/rem               "
    << std::setw(5) << reduceDBs << " "
    << std::setw(6) << counter->get_num_low_lbds() << " "
    << std::setw(6) << cls_removed);
  verb_print(1, "looks/look-computes            "
    << lookaheads << "/" << lookahead_computes);

  verb_print(1, "sat called/sat/unsat/conflK    "
    << std::setw(5) << sat_called << " "
    << std::setw(5) << sat_found_sat << " "
    << std::setw(5) << sat_found_unsat << " "
    << std::setw(5) << sat_conflicts/1000 << " ");
  verb_print(1, "buddy called                   "
    << std::setw(5) << buddy_called << " ");


  verb_print(1, "vivif: try/cls/clviv/litsravg  "
    << std::setw(9) << vivif_tried << " "
    << std::setw(9) << vivif_tried_cl << " "
    << std::setw(9) << vivif_cl_minim << " "
    << std::setw(9) << safe_div(vivif_cl_minim, vivif_lit_rem) << " "
  );

  verb_print(1, "toplev subs runs/bins/long-cls "
    << std::setw(9) << subsume_runs << " "
    << std::setw(9) << subsumed_bin_cls << " "
    << std::setw(9) << subsumed_cls << " "
  );
  verb_print(1, "toplevel probes/fail/bprop     "
    << std::setw(9) << toplevel_probe_runs << " "
    << std::setw(9) << toplevel_probe_fail << " "
    << std::setw(9) << toplevel_bothprop_fail << " "
  );

  verb_print(1, "probes/flits/bplits K          "
    << std::left
    << std::setw(6) << (num_failed_lit_tests_/1000ULL) << " "
    << std::setw(6) << (num_failed_literals_detected_/1000ULL) << " "
    << std::setw(6) << (num_failed_bprop_literals_failed/1000ULL) << " "
    << " -- " << std::setprecision(2) << safe_div( num_failed_literals_detected_, num_failed_lit_tests_)
    << " -- " << std::setprecision(2) << safe_div( num_failed_literals_detected_+num_failed_bprop_literals_failed, num_failed_lit_tests_)
    << std::setw(16) <<" -- Kprobe/s: "
    << std::setprecision(2) << std::setw(9) << std::left
    << safe_div(num_failed_lit_tests_,(1000.0*(cpuTime()-counter->get_start_time())))
  );
  verb_print(1, "implicit BCP miss rate         "
    << std::setprecision(2) << implicitBCP_miss_rate() * 100 << "%");
  verb_print(1, "cache entries K                " << (cache->get_num_entries_used()/1000ULL));
  verb_print(1, "MB cache                       "
    << std::setprecision(3) << in_mb(cache_bytes_memory_usage()) << " "
    << std::setprecision(3) << in_mb(cache->get_num_entries_used()*72) << " "
  );
  verb_print(1, "cache K (lookup/ stores/ hits) "
    << std::left
    << std::setw(6) << (num_cache_look_ups_/(1000ULL)) << " "
    << std::setw(6) << (total_num_cached_comps_ /(1000ULL)) << " "
    << std::setw(6) << (num_cache_hits_ /(1000ULL)) << " "
    << std::setw(16) << " -- Klookup/s: "
    << std::setprecision(2) << std::setw(9) << std::left
    << safe_div(num_cache_look_ups_,(1000.0*(cpuTime()-counter->get_start_time())))
  );
  verb_print(1, "cache pollutions call/removed  "
    << cache_pollutions_called << "/"
    << cache_pollutions_removed);
  verb_print(1, "cache miss rate                "
    << std::setprecision(3) << cache_miss_rate());
  verb_print(1, "avg hit/store num vars "
    << get_avg_cache_store_sz()
    << " / "
    << get_avg_cache_store_size());
  if (conf.verb >= 2) counter->get_cache().debug_mem_data();
  verb_print(1, "");
}


void DataAndStatistics::print_short_formula_info(const Counter* counter) const {
  verb_print(1, "irred cls (all/long/bin/unit): "
    << counter->get_num_irred_long_cls() << "/" << num_binary_irred_clauses_);
}
