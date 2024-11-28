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
#include "common.hpp"
#include "comp_cache.hpp"
#include "counter.hpp"
#include "time_mem.hpp"
#include "mpreal.h"

#include <gmpxx.h>
#include <iostream>
#include <iomanip>

using std::setw;

using namespace GanakInt;

static double in_mb(uint64_t bytes) {
  return (double)bytes/(double)(1024*1024);
}

static double safe_div(double a, double b) {
  if (b == 0) return 0;
  else return a/b;
}

template<typename T>
void DataAndStatistics<T>::print_short(const Counter<T>* counter, const CompCache<T>* cache) const {
  verb_print(1, "total time so far: " << cpu_time());
  verb_print(1, "cls long irred                 " << counter->get_num_irred_long_cls());
  verb_print(1, "decisions K                    "
    << std::left << setw(9) << decisions/1000
    << setw(16) << " -- Kdec/s: "
    << std::setprecision(2) << setw(9) << std::left << std::fixed
    << safe_div(decisions,(1000.0*(cpu_time()-counter->get_start_time())))
  );
  verb_print(1, "conflicts                      "
    << std::left << setw(9) << conflicts
    << "   " << std::left << setw(9)
    << setw(16) << " -- confl/s: "
    << std::setprecision(2) << setw(9) << std::left
    << safe_div(conflicts,((cpu_time()-counter->get_start_time())))
  );
  verb_print(1, "conflict cls (long/bin)        " << std::fixed
    << counter->get_num_long_red_cls() << "/" << num_bin_red_cls);

  /* verb_print(1, "lits /rem lits ccmin           " */
  /*   << setw(9) << orig_uip_lits/uip_cls << " " */
  /*   << setw(9) << ccmin_uip_lits/uip_cls << " " */
  /* ); */

  verb_print(1, "rem lits triedK/rem lits remK  "
    << setw(9) << rem_lits_tried/1000 << " "
    << setw(9) << rem_lits_with_bins/1000 << " "
  );

  verb_print(1, "avg clsz/rem lits avg/finalavg "
    << setw(9) << safe_div(orig_uip_lits,uip_cls) << " "
    << setw(9) << safe_div(uip_lits_ccmin,uip_cls) << " "
    << setw(9) << safe_div(final_cl_sz, uip_cls)
  );

  verb_print(1, "rdbs/low lbd/rem               "
    << setw(5) << reduce_db << " "
    << setw(6) << counter->get_num_low_lbds() << " "
    << setw(6) << cls_removed);
  verb_print(1, "looks/look-computes            "
    << lookaheads << "/" << lookahead_computes);

  verb_print(1, "sat call/sat/unsat/confl/rst  "
    << setw(5) << sat_called << " "
    << setw(5) << sat_found_sat << " "
    << setw(5) << sat_found_unsat << " "
    << setw(5) << sat_conflicts << " "
    << setw(5) << sat_rst << " ");
  verb_print(1, "buddy called /unsat ratio           "
    << setw(5) << buddy_called << " / "
    << setw(5) << safe_div(buddy_unsat, buddy_called) << " ");
  verb_print(1, "buddy avg num bin/long cls/vs  "
    << setw(5) << safe_div(buddy_num_bin_cls,buddy_called) << " / "
    << setw(5) << safe_div(buddy_num_long_cls,buddy_called) << " / "
    << setw(5) << safe_div(buddy_num_vars,buddy_called));
  verb_print(1, "buddy max num bin/long cls/vs  "
    << setw(5) << buddy_max_bin_cls << " / "
    << setw(5) << buddy_max_long_cls << " / "
    << setw(5) << buddy_max_num_vars);

  verb_print(1, "orig cubes/lit-r avg/exten avg "
    << setw(5) << num_cubes_orig << " / "
    << setw(5) << std::setprecision(2) << safe_div(cube_lit_rem, num_cubes_orig) << " / "
    << setw(5) << std::setprecision(2) << safe_div(cube_lit_extend, num_cubes_orig));
  verb_print(1, "cubes orig/symm/final          "
    << setw(5) << num_cubes_orig << " / "
    << setw(5) << num_cubes_symm << " / "
    << setw(5) << num_cubes_final);

  verb_print(1, "tot restarts                   " << setw(5) << num_restarts);

  verb_print(1, "comp sortsK/avg sz             "
    << setw(5) << comp_sorts/1000 << " / "
    << setw(5) << std::setprecision(8) << safe_div(comp_sizes, comp_sorts))
    << std::setprecision(2);

  verb_print(1, "vivif: try/cls/clviv/litsravg  "
    << setw(9) << vivif_tried << " "
    << setw(9) << vivif_tried_cl << " "
    << setw(9) << vivif_cl_minim << " "
    << setw(9) << safe_div(vivif_cl_minim, vivif_lit_rem) << " "
  );

  verb_print(1, "toplev subs runs/bins/long-cls "
    << setw(9) << subsume_runs << " "
    << setw(9) << subsumed_bin_irred_cls << " "
    << setw(9) << subsumed_bin_red_cls << " "
    << setw(9) << subsumed_long_irred_cls << " "
    << setw(9) << subsumed_long_red_cls << " "
  );
  verb_print(1, "toplevel probes/fail/bprop     "
    << setw(9) << toplevel_probe_runs << " "
    << setw(9) << toplevel_probe_fail << " "
    << setw(9) << toplevel_bothprop_fail << " "
  );

  double vm_usage = 0;
  verb_print(1, "Mem used                       "
    << std::setprecision(2) << (double)mem_used(vm_usage) / (1e9)  << " GB");
  verb_print(1, "cache pollutions call/removed  "
    << cache_pollutions_called << "/"
    << cache_pollutions_removed);
  verb_print(1, "cache entries K                " << (cache->get_num_entries_used()/1000ULL));
  verb_print(1, "MB cache                       "
    << std::setprecision(3) << in_mb(cache_bytes_memory_usage()) << " "
    << std::setprecision(3) << in_mb(cache->get_num_entries_used()*72) << " "
  );
  verb_print(1, "cache K (lookup/ stores/ hits/ dels) "
    << std::left
    << setw(6) << (num_cache_look_ups/(1000ULL)) << " "
    << setw(6) << (total_num_cached_comps /(1000ULL)) << " "
    << setw(6) << (num_cache_hits /(1000ULL)) << " "
    << setw(6) << (num_cache_dels /(1000ULL)) << " "
    << setw(16) << " -- Klookup/s: "
    << std::setprecision(2) << setw(9) << std::left
    << safe_div(num_cache_look_ups,(1000.0*(cpu_time()-counter->get_start_time())))
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


template<typename T>
void DataAndStatistics<T>::print_short_formula_info(const Counter<T>* counter) const {
  verb_print(1, "irred cls long/bin: "
    << counter->get_num_irred_long_cls() << "/" << num_bin_irred_cls);
}

template class DataAndStatistics<mpz_class>;
template class DataAndStatistics<mpfr::mpreal>;
template class DataAndStatistics<mpq_class>;
