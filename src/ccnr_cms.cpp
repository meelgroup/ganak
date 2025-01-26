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
#include "time_mem.hpp"

using namespace CCNR;

Ganak_ccnr::Ganak_ccnr(const ArjunNS::SimplifiedCNF* _cnf, uint32_t _verb): cnf(_cnf)  {
    conf.verb = _verb;
    ls_s = new LS_solver();
    ls_s->set_verbosity(conf.verb);
}

Ganak_ccnr::~Ganak_ccnr() { delete ls_s; }

int Ganak_ccnr::main()
{
    //It might not work well with few number of variables
    //rnovelty could also die/exit(-1), etc.
    if (cnf->nVars() < 50 || cnf->clauses.size() < 10) {
        verb_print(1, "[ccnr] too few variables & clauses");
        return 0;
    }
    double start_time = cpu_time();

    init_problem();

    vector<bool> phases(cnf->nVars()+1, false);
    int res = ls_s->local_search(&phases, conf.yalsat_max_mems*2*1000*1000, "c o");

    double time_used = cpu_time()-start_time;
    verb_print(1, "[ccnr] time: " << time_used);
    return res;
}

template<class T>
void Ganak_ccnr::add_this_clause(const T& cl) {
    uint32_t sz = 0;
    yals_lits.clear();
    for(size_t i3 = 0; i3 < cl.size(); i3++) {
      CMSat::Lit lit = cl[i3];
        int l = lit.var()+1;
        l *= lit.sign() ? -1 : 1;
        yals_lits.push_back(l);
        sz++;
    }
    assert(sz > 0);
    cl_num++;

    for(auto& lit: yals_lits) {
        ls_s->cls[cl_num].literals.push_back(CCNR::lit(lit, cl_num));
    }
}

bool Ganak_ccnr::init_problem() {
    cl_num = 0;
    ls_s->num_vars = cnf->nVars();
    ls_s->num_vars = cnf->clauses.size();
    ls_s->make_space();
    for(auto& cl: cnf->clauses) add_this_clause(cl);

    for (int c=0; c < ls_s->num_vars; c++) {
        for(CCNR::lit item: ls_s->cls[c].literals) {
            int v = item.var_num;
            ls_s->vars[v].literals.push_back(item);
        }
    }
    ls_s->build_neighborhood();
    return true;
}
