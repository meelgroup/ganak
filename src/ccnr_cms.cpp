/******************************************
Copyright (C) 2009-2020 Authors of CryptoMiniSat, see AUTHORS file <soos.mate@gmail.com>

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

#include <cstdint>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "ccnr_cms.hpp"
#include "ccnr.hpp"
#include "common.hpp"
#include <arjun/arjun.h>
#include <iomanip>
#include "time_mem.hpp"

using namespace CCNR;
using std::setprecision;
using std::fixed;

Ganak_ccnr::Ganak_ccnr(uint32_t _verb) {
    conf.verb = _verb;
    ls_s = new LSSolver();
    ls_s->set_verbosity(conf.verb);
}

Ganak_ccnr::~Ganak_ccnr() { delete ls_s; }

int Ganak_ccnr::main(const ArjunNS::SimplifiedCNF* _cnf) {
    cnf = _cnf;

    //It might not work well with few number of variables
    //rnovelty could also die/exit(-1), etc.
    if (cnf->nVars() == 0 || cnf->clauses.size() == 0) {
        verb_print(1, "[ccnr] too few variables & clauses");
        return 0;
    }
    double start_time = cpu_time();

    if (cnf->nVars() == cnf->get_sampl_vars().size()) {
        verb_print(1, "[ccnr] all vars are sampling vars, exiting");
        return 0;
    }
    init_problem();
    int res = ls_s->local_search(conf.yalsat_max_mems*50LL*1000LL, "c o");
    if (res == 1) {
      ls_s->print_solution(true);
    }

    double time_used = cpu_time()-start_time;
    verb_print(1, "[ccnr] T: " << setprecision(2) << fixed << time_used << " res: " << res);
    return res;
}

template<class T>
void Ganak_ccnr::add_this_clause(const T& cl) {
    uint32_t sz = 0;
    yals_lits.clear();
    for(auto& lit: cl) {
        int l = lit.var()+1;
        l *= lit.sign() ? -1 : 1;
        yals_lits.push_back(l);
        sz++;
    }
    assert(sz > 0);

    for(auto& lit: yals_lits) {
        ls_s->cls[cl_num].lits.push_back(CCNR::lit(lit, cl_num));
    }
    cl_num++;
}

bool Ganak_ccnr::init_problem() {
    ls_s->num_vars = cnf->nVars();
    ls_s->num_cls = cnf->clauses.size();
    ls_s->indep_map.clear();
    ls_s->indep_map.resize(cnf->nVars()+1, 0);
    for(const auto& v: cnf->get_sampl_vars()) ls_s->indep_map[v+1] = 1;
    ls_s->make_space();
    for(auto& cl: cnf->clauses) add_this_clause(cl);

    for (int c=0; c < ls_s->num_cls; c++) {
        for(auto& l: ls_s->cls[c].lits) {
            ls_s->vars[l.var_num].lits.push_back(l);
        }
    }
    ls_s->build_neighborhood();
    return true;
}
