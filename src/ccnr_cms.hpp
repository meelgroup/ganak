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

#pragma once

#include <cstdint>
#include <cstdio>
#include <utility>
#include <arjun/arjun.h>

using std::pair;
using std::make_pair;
using std::vector;

namespace CCNR {

class ls_solver;

struct CCNRConf {
  uint32_t verb = 1;
  uint64_t yalsat_max_mems = 1000;
};

class Ganak_ccnr {
public:
    int main();
    Ganak_ccnr(ArjunNS::SimplifiedCNF& _cnf, uint32_t verbosity);
    ~Ganak_ccnr();

private:
    void flipvar(uint32_t toflip);
    void parse_parameters();
    void init_for_round();
    bool init_problem();
    ls_solver* ls_s = nullptr;
    uint32_t cl_num = 0;

    enum class add_cl_ret {added_cl, skipped_cl, unsat};
    template<class T>
    void add_this_clause(const T& cl);
    vector<int> yals_lits;

    CCNRConf conf;
    ArjunNS::SimplifiedCNF& cnf;
};

}
