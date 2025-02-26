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
#include <arjun/arjun.h>
using namespace std;
using namespace GanakInt;

CounterConfiguration conf;
int do_pre_backbone = 0;
ArjunNS::SimpConf simp_conf;
ArjunNS::Arjun::ElimToFileConf etof_conf;

vector<Lit> cms_to_ganak_cl(const vector<CMSat::Lit>& cl) {
  vector<Lit> ganak_cl; ganak_cl.reserve(cl.size());
  for(const auto& l: cl) ganak_cl.push_back(Lit(l.var()+1, !l.sign()));
  return ganak_cl;
}

void run_arjun(ArjunNS::SimplifiedCNF& cnf) {
  ArjunNS::Arjun arjun;
  arjun.set_verb(0);
  if (do_pre_backbone) arjun.standalone_backbone(cnf);
  arjun.standalone_minimize_indep(cnf, false);
  assert(!etof_conf.all_indep);
  arjun.standalone_elim_to_file(cnf, etof_conf, simp_conf);
}

void setup_ganak(const ArjunNS::SimplifiedCNF& cnf, Ganak& counter) {
  // Setup independent support
  counter.new_vars(cnf.nVars());
  set<uint32_t> tmp;
  for(auto const& s: cnf.sampl_vars) tmp.insert(s+1);
  counter.set_indep_support(tmp);
  if (cnf.get_opt_sampl_vars_set()) {
    tmp.clear();
    for(auto const& s: cnf.opt_sampl_vars) tmp.insert(s+1);
    counter.set_optional_indep_support(tmp);
  }
  assert(!cnf.weighted);

  // Add clauses
  for(const auto& cl: cnf.clauses) counter.add_irred_cl(cms_to_ganak_cl(cl));
  for(const auto& cl: cnf.red_clauses) counter.add_red_cl(cms_to_ganak_cl(cl));
  counter.end_irred_cls();
}

int main() {
  ArjunNS::SimplifiedCNF cnf;
  cnf.new_vars(10);
  vector<CMSat::Lit> cl;
  cl.push_back(CMSat::Lit(0, false));
  cl.push_back(CMSat::Lit(1, false));
  cnf.add_clause(cl);
  cnf.set_weighted(false);
  cnf.set_sampl_vars(vector<uint32_t>{0, 1, 3});

  run_arjun(cnf);
  conf.verb = 0;
  Ganak counter(conf, cnf.weighted);
  setup_ganak(cnf, counter);

  mpz_class cnt;
  if (cnf.multiplier_weight == std::complex<mpq_class>()) cnt = mpz_class();
  else cnt = counter.unw_outer_count();
  assert(!counter.get_is_approximate());
  assert(cnf.multiplier_weight.imag() == 0);
  cnt *= cnf.multiplier_weight.real();

  cout << "c s exact arb int " << std::fixed << cnt << endl;
  return 0;
}
