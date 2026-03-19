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

#include <cryptominisat5/solvertypesmini.h>
#include <ganak.hpp>
#include <arjun/arjun.h>
#include <memory>
using namespace std;
using namespace GanakInt;

CounterConfiguration conf;
ArjunNS::SimpConf simp_conf;
ArjunNS::Arjun::ElimToFileConf etof_conf;

void run_arjun(ArjunNS::SimplifiedCNF& cnf) {
  ArjunNS::Arjun arjun;
  arjun.set_verb(0);
  arjun.standalone_minimize_indep(cnf, false);
  assert(!etof_conf.all_indep);
  arjun.standalone_elim_to_file(cnf, etof_conf, simp_conf);
}

constexpr CMSat::Lit mklit(int lit) {
  assert (lit != 0);
  if (lit > 0) return CMSat::Lit(lit - 1, false);
  else return CMSat::Lit((-lit) - 1, true);
}

int main() {
  std::unique_ptr<CMSat::FieldGen> fg = std::make_unique<ArjunNS::FGenMpq>();
  ArjunNS::SimplifiedCNF cnf(fg);
  cnf.new_vars(10);
  vector<CMSat::Lit> cl;
  cl.push_back(mklit(1));
  cl.push_back(mklit(2));
  cnf.add_clause(cl);
  cnf.set_weighted(true);
  cnf.set_lit_weight(mklit(4), ArjunNS::FMpq(10));
  cnf.set_lit_weight(mklit(-4), ArjunNS::FMpq(1));
  cnf.set_sampl_vars(vector<uint32_t>{mklit(1).var(), mklit(2).var(), mklit(4).var()});

  run_arjun(cnf);
  conf.verb = 0;
  Ganak counter(conf, fg);
  setup_ganak(cnf, counter);

  auto cnt = cnf.get_multiplier_weight()->dup();
  if (!cnf.get_multiplier_weight()->is_zero()) *cnt *= *counter.count();
  cout << "count is: " << std::fixed << *cnt << endl;
  return 0;
}
