/******************************************
Copyright (c) 2019, Shaowei Cai

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
#include <cstdlib>
#include <ostream>
#include <vector>
#include "ccnr_mersenne.hpp"

using std::vector;
using std::abs;

namespace CCNR {

struct lit {
    int var_num;             //variable num, begin with 1
    int sense : 1;           //is 1 for true literals, 0 for false literals.
    int cl_num : 31;         //clause ID it belongs to, begin with 0
    lit(int the_lit, int the_clause) {
        var_num = std::abs(the_lit);
        cl_num = the_clause;
        sense = the_lit > 0 ? 1 : 0;
    }

    struct lit &operator^=(const struct lit &l) {
        sense ^= l.sense;
        cl_num ^= l.cl_num;
        var_num ^= l.var_num;
        return *this;
    }

    void reset(void) {
        sense = 0;
        cl_num = 0;
        var_num = 0;
    }

    bool operator==(const struct lit &l) const {
        return sense == l.sense && cl_num == l.cl_num && var_num == l.var_num;
    }

    bool operator!=(const struct lit &l) const {
        return !(*this == l);
    }
};

struct variable {
    vector<lit> lits;
    vector<int> neighbor_vars;
    long long score;
    long long last_flip_step;
    int unsat_appear; //how many unsat clauses it appears in
};

struct clause {
    vector<lit> lits;
    int sat_count; //no. of satisfied literals
    int touched_cnt; // no. of literals that are touched but not satisfied
    int sat_var;
    long long weight;
};

inline std::ostream& operator<<(std::ostream &os, const lit &l) {
    os << (l.sense ? "" : "-") << l.var_num;
    return os;
}

inline std::ostream& operator<<(std::ostream &os, const clause &cl) {
  for(const auto &l: cl.lits) {
    os << l << " ";
  }
  os << "0";
  return os;
}

class LSSolver {
 public:
    LSSolver();
    bool local_search(
        long long int _mems_limit = 100*1000*1000
        , const char* prefix = "c "
    );
    void print_solution(bool need_verify = 0);
    void set_verbosity(uint32_t _verb) { verb = _verb; }

    //formula
    vector<variable> vars;
    vector<clause> cls;
    int num_vars;
    int num_cls;

    //data structure used
    vector<int> unsat_cls; // list of unsatisfied clauses
    vector<int> idx_in_unsat_cls; // idx_in_unsat_cls[var] tells where "var" is in unsat_vars
    vector<int> unsat_vars; // clauses are UNSAT due to these vars
    vector<int> idx_in_unsat_vars;
    vector<int> indep_map; // always num_vars+1 size, contains 0/1 if it's indep or not
    vector<int> touched_cls; // cls that have been touched but not satisfied
    uint32_t set_vars = 0; //vars that have been set
    vector<uint8_t> sol; //solution information. 0 = false, 1 = true, 3 = unset

    //functions for building data structure
    bool make_space();
    void build_neighborhood();
    int get_cost() { return unsat_cls.size(); }

  private:
    Mersenne random_gen;

    // Config
    long long max_steps;
    int max_tries;
    uint32_t verb = 0;
    float swt_p = 0.3;
    float swt_q = 0.7;
    int swt_thresh = 50;

    // internal stats
    long long step;
    long long mems = 0;
    int avg_cl_weight;
    long long delta_tot_cl_weight;

    //main functions
    void initialize();
    void initialize_variable_datas();
    int pick_var();
    void flip(int flipv);
    void update_clause_weights();
    void smooth_clause_weights();

    //funcitons for basic operations
    void sat_a_clause(int the_clause);
    void unsat_a_clause(int the_clause);
};

} // namespace CCNR
