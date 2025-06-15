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

#include "ganak.hpp"
#include "counter.hpp"

using namespace GanakInt;

#if defined _WIN32
    #define DLL_PUBLIC __declspec(dllexport)
#else
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#endif

DLL_PUBLIC Ganak::Ganak(CounterConfiguration& conf, FG& fg) {
  counter = new OuterCounter(conf, fg);
}

DLL_PUBLIC Ganak::~Ganak() {
  OuterCounter* c = (OuterCounter*)counter;
  delete c;
  counter = nullptr;
}

DLL_PUBLIC FF Ganak::count() {
  OuterCounter* c = (OuterCounter*)counter;
  return c->outer_count();
}

DLL_PUBLIC void Ganak::set_generators(const std::vector<std::map<GanakInt::Lit, GanakInt::Lit>>& _gens) {
  OuterCounter* c = (OuterCounter*)counter;
  c->set_generators(_gens);
}
DLL_PUBLIC void Ganak::end_irred_cls() {
  OuterCounter* c = (OuterCounter*)counter;
  c->end_irred_cls();
}
DLL_PUBLIC void Ganak::set_indep_support(const std::set<uint32_t>& indeps) {
  OuterCounter* c = (OuterCounter*)counter;
  c->set_indep_support(indeps);
}
DLL_PUBLIC bool Ganak::add_red_cl(const std::vector<GanakInt::Lit>& lits, int lbd) {
  OuterCounter* c = (OuterCounter*)counter;
  return c->add_red_cl(lits, lbd);
}
DLL_PUBLIC bool Ganak::get_is_approximate() const {
  OuterCounter* c = (OuterCounter*)counter;
  return c->get_is_approximate();
}
DLL_PUBLIC bool Ganak::add_irred_cl(const std::vector<GanakInt::Lit>& lits) {
  OuterCounter* c = (OuterCounter*)counter;
  return c->add_irred_cl(lits);
}
DLL_PUBLIC void Ganak::set_optional_indep_support(const std::set<uint32_t>& indeps) {
  OuterCounter* c = (OuterCounter*)counter;
  c->set_optional_indep_support(indeps);
}
DLL_PUBLIC void Ganak::set_lit_weight(const GanakInt::Lit l, const FF& w) {
  OuterCounter* c = (OuterCounter*)counter;
  c->set_lit_weight(l, w);
}
DLL_PUBLIC void Ganak::new_vars(const uint32_t n) {
  OuterCounter* c = (OuterCounter*)counter;
  c->new_vars(n);
}
DLL_PUBLIC void Ganak::print_indep_distrib() const {
  OuterCounter* c = (OuterCounter*)counter;
  c->print_indep_distrib();
}

DLL_PUBLIC uint64_t Ganak::get_num_cache_lookups() const {
  OuterCounter* c = (OuterCounter*)counter;
  return c->get_num_cache_lookups();
}

DLL_PUBLIC uint64_t Ganak::get_max_cache_elems() const {
  OuterCounter* c = (OuterCounter*)counter;
  return c->get_max_cache_elems();
}
