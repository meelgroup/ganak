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
#include "counter_config.hpp"
#include "lit.hpp"
#include <vector>
#include <map>
#include <set>
#include <gmpxx.h>
#include <complex>

class Ganak {
public:
  Ganak(GanakInt::CounterConfiguration& conf, bool weighted, bool cpx = false );
  ~Ganak();
  void* count();
  void new_vars(const uint32_t n);
  bool add_irred_cl(const std::vector<GanakInt::Lit>& lits);
  void end_irred_cls();

  void set_indep_support(const std::set<uint32_t>& indeps);
  void set_generators(const std::vector<std::map<GanakInt::Lit, GanakInt::Lit>>& _gens);

  std::complex<mpq_class> cpx_outer_count();
  mpz_class unw_outer_count();
  mpq_class wq_outer_count();

  bool add_red_cl(const std::vector<GanakInt::Lit>& lits, int lbd = -1);
  bool get_is_approximate() const;
  void set_optional_indep_support(const std::set<uint32_t>& indeps);
  void set_lit_weight(const GanakInt::Lit l, const std::complex<mpq_class>& w);
  void print_indep_distrib() const;
private:
  void* counter = nullptr;
};
