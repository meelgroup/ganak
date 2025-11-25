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
  OuterCounter(const CounterConfiguration& conf, const FG& fg) {
    counter = std::make_unique<Counter>(conf, fg);
  }

  void set_generators(const vector<map<Lit, Lit>>& _gens) { counter->set_generators(_gens); }
  void end_irred_cls() { counter->end_irred_cls(); }
  void set_indep_support(const set<uint32_t>& indeps) { counter->set_indep_support(indeps); }
  FF count() { return counter->outer_count();}
  bool add_red_cl(const vector<Lit>& lits, int lbd = -1) { return counter->add_red_cl(lits, lbd); }
  bool get_is_approximate() const { return counter->get_is_approximate();}
  bool add_irred_cl(const vector<Lit>& lits) { return counter->add_irred_cl(lits); }
  void set_optional_indep_support(const set<uint32_t>& indeps) {
    counter->set_optional_indep_support(indeps);
  }
  void set_lit_weight(const Lit l, const FF& w) {
    return counter->set_lit_weight(l, w);
  }
  void new_vars(const uint32_t n) { counter->new_vars(n); }
  void print_indep_distrib() const { counter->print_indep_distrib(); }
  uint64_t get_num_cache_lookups() const {
    return counter->get_stats().num_cache_look_ups;
  }

  uint64_t get_max_cache_elems() const {
    return counter->get_cache()->get_max_num_entries();
  }
private:
  unique_ptr<Counter> counter = nullptr;
};

}
