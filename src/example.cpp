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

#include "counter.hpp"
#include <arjun/arjun.h>
using namespace std;
CounterConfiguration conf;
int arjun_gates = 1;
int sbva_steps = 1000;
int sbva_cls_cutoff = 4;
int sbva_lits_cutoff = 5;
int sbva_tiebreak = 1;
int do_bce = 1;
int all_indep = 0;
int arjun_extend_max_confl = 1000;
int do_extend_indep = 1;
int do_backbone = 0;
int do_probe_based = 1;
int arjun_simp_level = 2;
int arjun_backw_maxc = 20000;
ArjunNS::SimpConf simp_conf;
string debug_arjun_cnf;
int do_precise = 1;

int arjun_oracle_find_bins = 0;

vector<Lit> cms_to_ganak_cl(const vector<CMSat::Lit>& cl) {
  vector<Lit> ganak_cl;
  for(const auto& l: cl) ganak_cl.push_back(Lit(l.var()+1, !l.sign()));
  return ganak_cl;
}


void run_arjun(ArjunNS::SimplifiedCNF& cnf) {
  ArjunNS::Arjun arjun;
  arjun.set_or_gate_based(arjun_gates);
  arjun.set_xor_gates_based(arjun_gates);
  arjun.set_ite_gate_based(arjun_gates);
  arjun.set_irreg_gate_based(arjun_gates);
  arjun.set_extend_max_confl(arjun_extend_max_confl);
  arjun.set_probe_based(do_probe_based);
  arjun.set_simp(arjun_simp_level);
  arjun.set_backw_max_confl(arjun_backw_maxc);
  arjun.set_oracle_find_bins(arjun_oracle_find_bins);
  arjun.set_verb(0);
  if (do_backbone) arjun.only_backbone(cnf);
  arjun.only_run_minimize_indep(cnf);
  bool do_unate = false;
  assert(!all_indep);
  arjun.elim_to_file(cnf, all_indep, do_extend_indep, do_bce, do_unate,
      simp_conf, sbva_steps, sbva_cls_cutoff, sbva_lits_cutoff, sbva_tiebreak);
}

void setup_ganak(const ArjunNS::SimplifiedCNF& cnf, OuterCounter& counter) {
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
  OuterCounter counter(conf, cnf.weighted, do_precise);
  setup_ganak(cnf, counter);

  mpz_class cnt;
  if (cnf.multiplier_weight == 0) cnt = 0;
  else cnt = counter.unw_outer_count();
  assert(!counter.get_is_approximate());
  cnt *= cnf.multiplier_weight;

  cout << "c s exact arb int " << std::fixed << cnt << endl;
  return 0;
}
