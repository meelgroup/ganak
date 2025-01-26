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
    unsigned char sense : 1; //is 1 for true literals, 0 for false literals.
    int clause_num : 31;     //clause num, begin with 0
    int var_num;             //variable num, begin with 1
    lit(int the_lit, int the_clause) {
        var_num = std::abs(the_lit);
        clause_num = the_clause;
        sense = the_lit > 0 ? 1 : 0;
    }

    struct lit &operator^=(const struct lit &l) {
        sense ^= l.sense;
        clause_num ^= l.clause_num;
        var_num ^= l.var_num;
        return *this;
    }

    void reset(void) {
        sense = 0;
        clause_num = 0;
        var_num = 0;
    }

    bool operator==(const struct lit &l) const {
        return sense == l.sense && clause_num == l.clause_num && var_num == l.var_num;
    }

    bool operator!=(const struct lit &l) const {
        return !(*this == l);
    }
};

struct variable {
    vector<lit> literals;
    vector<int> neighbor_var_nums;
    long long score;
    long long last_flip_step;
    int unsat_appear; //how many unsat clauses it appears in
};

struct clause {
    vector<lit> literals;
    int sat_count; //no. of satisfied literals
    int sat_var;
    long long weight;
};

inline std::ostream& operator<<(std::ostream &os, const lit &l) {
    os << (l.sense ? "" : "-") << l.var_num;
    return os;
}

inline std::ostream& operator<<(std::ostream &os, const clause &cl) {
  for(const auto &l: cl.literals) {
    os << l << " ";
  }
  os << "0";
  return os;
}


class LS_solver {
 public:
    LS_solver();
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
    vector<int> conflict_ct;
    vector<int> unsat_cls; // list of unsatisfied clauses
    vector<int> idx_in_unsat_cls; // idx_in_unsat_cls[var] tells where "var" is in unsat_vars
    vector<int> unsat_vars; // clauses are UNSAT due to these vars
    vector<int> idx_in_unsat_vars;

    vector<uint8_t> sol; //solution information

    //functions for buiding data structure
    bool make_space();
    void build_neighborhood();
    int get_cost() { return unsat_cls.size(); }

  private:
    long long mems = 0;
    long long step;
    long long max_steps;
    int max_tries;

    //aiding data structure
    Mersenne random_gen; //random generator

    //clause weighting
    int swt_thresh;
    float swt_p; //w=w*p+ave_w*q
    float swt_q;
    int avg_cl_weight;
    //-------------------

    //=================
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

    //--------------------
    uint32_t verb = 0;
};

} // namespace CCNR
