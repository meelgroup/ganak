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

#include "ccnr.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cassert>
#include <cstdint>

using namespace CCNR;

using std::cout;
using std::endl;
using std::string;


//constructor with default setting.
LS_solver::LS_solver() {
    max_tries = 100;
    max_steps = 1*1000 * 1000;
    random_gen.seed(1337); // random seed
    swt_thresh = 50;
    swt_p = 0.3;
    swt_q = 0.7;
    verb = 0;
}

bool LS_solver::make_space() {
    if (num_vars == 0) return false;
    vars.resize(num_vars+1);
    clauses.resize(num_vars+1);
    solution.resize(num_vars+1);
    idx_in_unsat_cls.resize(num_vars+1);
    idx_in_unsat_vars.resize(num_vars+1);

    return true;
}

void LS_solver::build_neighborhood() {
    vector<bool> neighbor_flag(num_vars+1, 0);
    for (int v = 1; v <= num_vars; ++v) {
        variable *vp = &(vars[v]);
        for (lit lv: vp->literals) {
            int c = lv.clause_num;
            for (lit lc: clauses[c].literals) {
                if (!neighbor_flag[lc.var_num] && lc.var_num != v) {
                    neighbor_flag[lc.var_num] = 1;
                    vp->neighbor_var_nums.push_back(lc.var_num);
                }
            }
        }
        for (uint32_t j = 0; j < vp->neighbor_var_nums.size(); ++j) {
            neighbor_flag[vp->neighbor_var_nums[j]] = 0;
        }
    }
}

bool LS_solver::local_search(
    const vector<bool> *init_solution
    , long long int _mems_limit
    , const char* prefix
) {
    bool result = false;
    for (int t = 0; t < max_tries; t++) {
        initialize(init_solution);
        if (0 == unsat_cls.size()) {
            result = true;
            break;
        }

        for (step = 0; step < max_steps; step++) {
            int flipv = pick_var();
            flip(flipv);
            if (mems > _mems_limit) return result;

            int cost = unsat_cls.size();
            if (verb && (cost == 0 || (step & 0x3ffff) == 0x3ffff)) {
                cout << prefix << "[ccnr] tries: "
                << t << " steps: " << step
                << " best found: " << cost
                << endl;
            }

            if (cost == 0) {
                result = true;
                break;
            }
        }
        if (0 == unsat_cls.size()) {
            result = true;
            break;
        }
    }
    return result;
}

void LS_solver::clear_prev_data() {
    unsat_cls.clear();
    _ccd_vars.clear();
    unsat_vars.clear();
    for (int &item: idx_in_unsat_cls)
        item = 0;
    for (int &item: idx_in_unsat_vars)
        item = 0;
}

void LS_solver::initialize(const vector<bool> *init_solution) {
    clear_prev_data();
    if (!init_solution) {
        //default random generation
        for (int v = 1; v <= num_vars; v++) {
            solution[v] = (random_gen.next(2) == 0 ? 0 : 1);
        }
    } else {
        if ((int)init_solution->size() != num_vars+1) {
            cout
            << "ERROR: the init solution's size"
            " is not equal to the number of variables."
            << endl;
            exit(-1);
        }
        for (int v = 1; v <= num_vars; v++) {
            solution[v] = init_solution->at(v);
        }
    }

    //unsat_appears, will be updated when calling unsat_a_clause function.
    for (int v = 1; v <= num_vars; v++) {
        vars[v].unsat_appear = 0;
    }

    //initialize data structure of clauses according to init solution
    for (int c = 0; c < num_vars; c++) {
        clauses[c].sat_count = 0;
        clauses[c].sat_var = -1;
        clauses[c].weight = 1;

        for (lit l: clauses[c].literals) {
            if (solution[l.var_num] == l.sense) {
                clauses[c].sat_count++;
                clauses[c].sat_var = l.var_num;
            }
        }
        if (0 == clauses[c].sat_count) {
            unsat_a_clause(c);
        }
    }
    avg_cl_weight = 1;
    delta_tot_cl_weight = 0;
    initialize_variable_datas();
}
void LS_solver::initialize_variable_datas()
{
    variable *vp;
    //scores
    for (int v = 1; v <= num_vars; v++) {
        vp = &(vars[v]);
        vp->score = 0;
        for (lit l: vp->literals) {
            int c = l.clause_num;
            if (0 == clauses[c].sat_count) {
                vp->score += clauses[c].weight;
            } else if (1 == clauses[c].sat_count && l.sense == solution[l.var_num]) {
                vp->score -= clauses[c].weight;
            }
        }
    }
    //last flip step
    for (int v = 1; v <= num_vars; v++) {
        vars[v].last_flip_step = 0;
    }
    //cc datas
    for (int v = 1; v <= num_vars; v++) {
        vp = &(vars[v]);
        vp->cc_value = 1;
        if (vp->score > 0) //&&vars[v].cc_value==1
        {
            _ccd_vars.push_back(v);
            vp->is_in_ccd_vars = 1;
        } else {
            vp->is_in_ccd_vars = 0;
        }
    }
    //the virtual var 0
    vp = &(vars[0]);
    vp->score = 0;
    vp->cc_value = 0;
    vp->is_in_ccd_vars = 0;
    vp->last_flip_step = 0;
}


int LS_solver::pick_var() {
    //First, try to get the var with the highest score from _ccd_vars if any
    //----------------------------------------
    int best_var = 0;
    mems += _ccd_vars.size()/8;
    if (_ccd_vars.size() > 0) {
        best_var = _ccd_vars[0];
        for (int v: _ccd_vars) {
            if (vars[v].score > vars[best_var].score) {
                best_var = v;
            } else if (vars[v].score == vars[best_var].score &&
                       vars[v].last_flip_step < vars[best_var].last_flip_step) {
                best_var = v;
            }
        }
        return best_var;
    }

    /**Diversification Mode**/
    update_clause_weights();

    /*focused random walk*/
    int c = unsat_cls[random_gen.next(unsat_cls.size())];
    clause *cp = &(clauses[c]);
    best_var = cp->literals[0].var_num;
    for (size_t k = 1; k < cp->literals.size(); k++) {
        int v = cp->literals[k].var_num;
        if (vars[v].score > vars[best_var].score) {
            best_var = v;
        } else if (vars[v].score == vars[best_var].score &&
                   vars[v].last_flip_step < vars[best_var].last_flip_step) {
            best_var = v;
        }
    }
    return best_var;
}

void LS_solver::flip(int flipv) {
    solution[flipv] = 1 - solution[flipv];
    int org_flipv_score = vars[flipv].score;
    mems += vars[flipv].literals.size();

    // Go through each clause the literal is in and update status
    for (lit l: vars[flipv].literals) {
        clause *cp = &(clauses[l.clause_num]);
        if (solution[flipv] == l.sense) {
            cp->sat_count++;
            if (1 == cp->sat_count) {
                sat_a_clause(l.clause_num);
                cp->sat_var = flipv;
                for (lit lc: cp->literals) {
                    vars[lc.var_num].score -= cp->weight;
                }
            } else if (2 == cp->sat_count) {
                vars[cp->sat_var].score += cp->weight;
            }
        } else {
            cp->sat_count--;
            if (0 == cp->sat_count) {
                unsat_a_clause(l.clause_num);
                for (lit lc: cp->literals) {
                    vars[lc.var_num].score += cp->weight;
                }
            } else if (1 == cp->sat_count) {
                for (lit lc: cp->literals) {
                    if (solution[lc.var_num] == lc.sense) {
                        vars[lc.var_num].score -= cp->weight;
                        cp->sat_var = lc.var_num;
                        break;
                    }
                }
            }
        }
    }
    vars[flipv].score = -org_flipv_score;
    vars[flipv].last_flip_step = step;
    //update cc_values
    update_cc_after_flip(flipv);
}

void LS_solver::update_cc_after_flip(int flipv) {
    int last_item;
    variable *vp = &(vars[flipv]);
    vp->cc_value = 0;
    mems += _ccd_vars.size()/4;
    for (int index = _ccd_vars.size() - 1; index >= 0; index--) {
        int v = _ccd_vars[index];
        if (vars[v].score <= 0) {
            last_item = _ccd_vars.back();
            _ccd_vars.pop_back();
            if (index < (int)_ccd_vars.size()) {
                _ccd_vars[index] = last_item;
            }

            vars[v].is_in_ccd_vars = 0;
        }
    }

    //update all flipv's neighbor's cc to be 1
    mems += vp->neighbor_var_nums.size()/4;
    for (int v: vp->neighbor_var_nums) {
        vars[v].cc_value = 1;
        if (vars[v].score > 0 && !(vars[v].is_in_ccd_vars)) {
            _ccd_vars.push_back(v);
            vars[v].is_in_ccd_vars = 1;
        }
    }
}

void LS_solver::sat_a_clause(int the_clause) {
    //use the position of the clause to store the last unsat clause in stack
    int last_item = unsat_cls.back();
    unsat_cls.pop_back();
    int index = idx_in_unsat_cls[the_clause];
    if (index < (int)unsat_cls.size()) {
        unsat_cls[index] = last_item;
    }
    idx_in_unsat_cls[last_item] = index;
    //update unsat_appear and unsat_vars
    for (lit l: clauses[the_clause].literals) {
        vars[l.var_num].unsat_appear--;
        if (0 == vars[l.var_num].unsat_appear) {
            last_item = unsat_vars.back();
            unsat_vars.pop_back();
            index = idx_in_unsat_vars[l.var_num];
            if (index < (int)unsat_vars.size()) {
                unsat_vars[index] = last_item;
            }
            idx_in_unsat_vars[last_item] = index;
        }
    }
}

void LS_solver::unsat_a_clause(int the_clause) {
    idx_in_unsat_cls[the_clause] = unsat_cls.size();
    unsat_cls.push_back(the_clause);
    //update unsat_appear and unsat_vars
    for (lit l: clauses[the_clause].literals) {
        vars[l.var_num].unsat_appear++;
        if (1 == vars[l.var_num].unsat_appear) {
            idx_in_unsat_vars[l.var_num] = unsat_vars.size();
            unsat_vars.push_back(l.var_num);
        }
    }
}

void LS_solver::update_clause_weights() {
    for (int c: unsat_cls) {
        clauses[c].weight++;
    }
    for (int v: unsat_vars) {
        vars[v].score += vars[v].unsat_appear;
        if (vars[v].score > 0 && 1 == vars[v].cc_value && !(vars[v].is_in_ccd_vars)) {
            _ccd_vars.push_back(v);
            vars[v].is_in_ccd_vars = 1;
        }
    }
    delta_tot_cl_weight += unsat_cls.size();
    if (delta_tot_cl_weight >= num_vars) {
        avg_cl_weight += 1;
        delta_tot_cl_weight -= num_vars;
        if (avg_cl_weight > swt_thresh) {
            smooth_clause_weights();
        }
    }
}

void LS_solver::smooth_clause_weights() {
    for (int v = 1; v <= num_vars; v++) {
        vars[v].score = 0;
    }
    int scale_avg = avg_cl_weight * swt_q;
    avg_cl_weight = 0;
    delta_tot_cl_weight = 0;
    mems += num_vars;
    for (int c = 0; c < num_vars; ++c) {
        clause *cp = &(clauses[c]);
        cp->weight = cp->weight * swt_p + scale_avg;
        if (cp->weight < 1)
            cp->weight = 1;
        delta_tot_cl_weight += cp->weight;
        if (delta_tot_cl_weight >= num_vars) {
            avg_cl_weight += 1;
            delta_tot_cl_weight -= num_vars;
        }
        if (0 == cp->sat_count) {
            for (lit l: cp->literals) {
                vars[l.var_num].score += cp->weight;
            }
        } else if (1 == cp->sat_count) {
            vars[cp->sat_var].score -= cp->weight;
        }
    }

    //reset ccd_vars
    _ccd_vars.clear();
    for (int v = 1; v <= num_vars; v++) {
        variable* vp = &(vars[v]);
        if (vp->score > 0 && 1 == vp->cc_value) {
            _ccd_vars.push_back(v);
            vp->is_in_ccd_vars = 1;
        } else {
            vp->is_in_ccd_vars = 0;
        }
    }
}

void LS_solver::print_solution(bool need_verify) {
    if (0 == get_cost()) cout << "s SATISFIABLE" << endl;
    else cout << "s UNKNOWN" << endl;

    bool sat_flag = false;
    if (need_verify) {
        for (int c = 0; c < num_vars; c++) {
            sat_flag = false;
            for (lit l: clauses[c].literals) {
                if (solution[l.var_num] == l.sense) {
                    sat_flag = true;
                    break;
                }
            }
            if (!sat_flag) {
                cout << "c Error: verify error in clause " << c << endl;
                return;
            }
        }
        cout << "c Verified." << endl;
    }

    if (verb > 0) {
        cout << "v";
        for (int v = 1; v <= num_vars; v++) {
            cout << ' ';
            if (solution[v] == 0)
                cout << '-';
            cout << v;
        }
        cout << endl;
    }
}

