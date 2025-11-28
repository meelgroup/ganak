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
#include <memory>
#include <map>
#include <set>
#include <vector>

#include "counter.hpp"
#include "counter_config.hpp"
#include "common.hpp"

using std::unique_ptr;
using std::vector;
using std::map;
using std::set;

namespace GanakInt {

class OuterCounter {
public:
  OuterCounter(const CounterConfiguration& _conf, const FG& _fg) :
    conf(_conf), fg(_fg->dup()), nvars(0) {
  }

  void new_vars(const uint32_t n) { nvars = n; }
  void set_generators(const vector<map<Lit, Lit>>& _gens) { generators = _gens; }
  void set_indep_support(const set<uint32_t>& indeps) { indep_support = indeps; }
  void set_optional_indep_support(const set<uint32_t>& indeps) { opt_indep_support = indeps; }
  void add_red_cl(const vector<Lit>& lits, int lbd = -1) { red_cls.push_back({lits, (uint32_t)lbd}); }
  void add_irred_cl(const vector<Lit>& lits) { irred_cls.push_back(lits); }
  void set_lit_weight(const Lit l, const FF& w) { lit_weights[l] = w->dup(); }
  FF count(uint8_t bits_jobs = 0, int num_threads = 1);

  void print_indep_distrib() const;
  uint64_t get_num_cache_lookups() const { return num_cache_lookups; }
  uint64_t get_max_cache_elems() const { return max_cache_elems; }
  bool get_is_approximate() const { return count_is_approximate; }

private:
  FF count_with_parallel(uint8_t bits_jobs = 0, int num_threads = 1);
  FF count_regular();

  CounterConfiguration conf;
  FG fg;
  uint32_t nvars;

  vector<vector<Lit>> irred_cls;
  vector<std::pair<vector<Lit>, uint32_t>> red_cls;
  set<uint32_t> indep_support;
  set<uint32_t> opt_indep_support;
  map<Lit, FF> lit_weights;
  vector<map<Lit, Lit>> generators;

  // stats
  bool count_is_approximate = false;
  uint64_t num_cache_lookups = 0;
  uint64_t max_cache_elems = 0;
};

}
