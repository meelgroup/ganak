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
#include <cassert>

using namespace CCNR;

using std::cout;
using std::endl;
using std::string;


LSSolver::LSSolver() {
    max_tries = 100;
    max_steps = 1*1000 * 1000;
    random_gen.seed(1337);
    verb = 0;
}

bool LSSolver::make_space() {
    if (num_vars == 0) return false;
    vars.resize(num_vars+1);
    sol.resize(num_vars+1, 3);
    idx_in_unsat_vars.resize(num_vars+1);

    cls.resize(num_cls);
    idx_in_unsat_cls.resize(num_cls);

    return true;
}

void LSSolver::build_neighborhood() {
    vector<uint8_t> flag(num_vars+1, 0);
    for (int v = 1; v <= num_vars; ++v) {
        variable& vp = vars[v];
        for (lit lv: vp.literals) {
            int c = lv.clause_num;
            for (lit lc: cls[c].lits) {
                if (!flag[lc.var_num] && lc.var_num != v) {
                    flag[lc.var_num] = 1;
                    vp.neighbor_var_nums.push_back(lc.var_num);
                }
            }
        }
        for (int neighbor_var_num : vp.neighbor_var_nums)
          flag[neighbor_var_num] = 0;
    }
}

bool LSSolver::local_search(long long int _mems_limit , const char* prefix) {
    bool result = false;
    for (int t = 0; t < max_tries; t++) {
        initialize();
        if (unsat_cls.empty()) {
            result = true;
            break;
        }

        for (step = 0; step < max_steps; step++) {
            int flipv = pick_var();
            flip(flipv);
            if (mems > _mems_limit) return false;

            int cost = unsat_cls.size();
            if (verb && (cost == 0 || (step & 0x3ffff) == 0x3ffff)) {
                cout << prefix << "[ccnr] tries: "
                << t << " steps: " << step
                << " best found: " << cost
                << endl;
            }

            if (cost == 0) {
                print_solution(1);
                result = true;
                break;
            }
        }
        if (unsat_cls.empty()) {
            result = true;
            break;
        }
    }
    return result;
}

void LSSolver::initialize() {
    unsat_cls.clear();
    unsat_vars.clear();
    for (auto &i: idx_in_unsat_cls) i = 0;
    for (auto &i: idx_in_unsat_vars) i = 0;
    for (int v = 1; v <= num_vars; v++) {
      if (!indep_map[v]) sol[v] = random_gen.next(3);
      else assert(sol[v] == 3);
    }

    //unsat_appears, will be updated when calling unsat_a_clause function.
    for (int v = 1; v <= num_vars; v++) vars[v].unsat_appear = 0;

    //initialize data structure of clauses according to init solution
    for (int cid = 0; cid < num_cls; cid++) {
        auto& cl = cls[cid];
        cl.sat_count = 0;
        cl.touched_cnt = 0;
        cl.sat_var = -1;
        cl.weight = 1;

        for (lit l: cl.lits) {
            const auto val = sol[l.var_num];
            if (val == l.sense) {
                cl.sat_count++;
                cl.sat_var = l.var_num;
            } else if (val != 3) cl.touched_cnt++;
        }
        if (cl.sat_count == 0) unsat_a_clause(cid);
        if (cl.touched_cnt) touched_cls.push_back(cid); // <- touch_a_clause(cid);
    }
    avg_cl_weight = 1;
    delta_tot_cl_weight = 0;
    initialize_variable_datas();
}

void LSSolver::initialize_variable_datas() {
    //scores
    for (int v = 1; v <= num_vars; v++) {
        auto & vp = vars[v];
        vp.score = 0;
        for (lit l: vp.literals) {
            int c = l.clause_num;
            if (cls[c].sat_count == 0) {
                vp.score += cls[c].weight;
            } else if (1 == cls[c].sat_count && l.sense == sol[l.var_num]) {
                vp.score -= cls[c].weight;
            }
        }
    }
    //last flip step
    for (int v = 1; v <= num_vars; v++) vars[v].last_flip_step = 0;

    //the virtual var 0
    auto& vp = vars[0];
    vp.score = 0;
    vp.last_flip_step = 0;
}

int LSSolver::pick_var() {
    update_clause_weights();

    assert(!unsat_cls.empty() || !touched_cls.empty());
    int cid;
    if (!unsat_cls.empty() && !touched_cls.empty()) {
      cid = unsat_cls[random_gen.next(unsat_cls.size())];
    } else if (!unsat_cls.empty()) {
      cid = touched_cls[random_gen.next(touched_cls.size())];
    } else {
      cid = unsat_cls[random_gen.next(unsat_cls.size())];
    }
    clause& cl = cls[cid];
    int best_var = cl.lits[0].var_num;
    for (size_t k = 1; k < cl.lits.size(); k++) {
        int v = cl.lits[k].var_num;
        if (vars[v].score > vars[best_var].score) {
            best_var = v;
        } else if (vars[v].score == vars[best_var].score &&
                   vars[v].last_flip_step < vars[best_var].last_flip_step) {
            best_var = v;
        }
    }
    return best_var;
}

void LSSolver::flip(int flipv) {
    sol[flipv] = 1 - sol[flipv];
    int org_flipv_score = vars[flipv].score;
    mems += vars[flipv].literals.size();

    // Go through each clause the literal is in and update status
    for (lit l: vars[flipv].literals) {
        clause& cl = cls[l.clause_num];
        if (sol[flipv] == l.sense) {
            cl.sat_count++;
            if (1 == cl.sat_count) {
                sat_a_clause(l.clause_num);
                cl.sat_var = flipv;
                for (lit lc: cl.lits) {
                    vars[lc.var_num].score -= cl.weight;
                }
            } else if (2 == cl.sat_count) {
                vars[cl.sat_var].score += cl.weight;
            }
        } else {
            cl.sat_count--;
            if (0 == cl.sat_count) {
                unsat_a_clause(l.clause_num);
                for (lit lc: cl.lits) {
                    vars[lc.var_num].score += cl.weight;
                }
            } else if (1 == cl.sat_count) {
                for (lit lc: cl.lits) {
                    if (sol[lc.var_num] == lc.sense) {
                        vars[lc.var_num].score -= cl.weight;
                        cl.sat_var = lc.var_num;
                        break;
                    }
                }
            }
        }
    }
    vars[flipv].score = -org_flipv_score;
    vars[flipv].last_flip_step = step;
}


void LSSolver::sat_a_clause(int cl_num) {
    //use the position of the clause to store the last unsat clause in stack
    int last_item = unsat_cls.back();
    unsat_cls.pop_back();
    int index = idx_in_unsat_cls[cl_num];
    if (index < (int)unsat_cls.size()) unsat_cls[index] = last_item;
    idx_in_unsat_cls[last_item] = index;

    //update unsat_appear and unsat_vars
    for (lit l: cls[cl_num].lits) {
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

void LSSolver::unsat_a_clause(int cl_id) {
    idx_in_unsat_cls[cl_id] = unsat_cls.size();
    unsat_cls.push_back(cl_id);
    //update unsat_appear and unsat_vars
    for (lit l: cls[cl_id].lits) {
        vars[l.var_num].unsat_appear++;
        if (1 == vars[l.var_num].unsat_appear) {
            idx_in_unsat_vars[l.var_num] = unsat_vars.size();
            unsat_vars.push_back(l.var_num);
        }
    }
}

void LSSolver::update_clause_weights() {
    for (int c: unsat_cls) cls[c].weight++;
    for (int v: unsat_vars) vars[v].score += vars[v].unsat_appear;

    delta_tot_cl_weight += unsat_cls.size();
    if (delta_tot_cl_weight >= num_cls) {
        avg_cl_weight += 1;
        delta_tot_cl_weight -= num_cls;
        if (avg_cl_weight > swt_thresh) {
            smooth_clause_weights();
        }
    }
}

void LSSolver::smooth_clause_weights() {
    for (int v = 1; v <= num_vars; v++) vars[v].score = 0;
    int scale_avg = avg_cl_weight * swt_q;
    avg_cl_weight = 0;
    delta_tot_cl_weight = 0;
    mems += num_cls;
    for (int c = 0; c < num_cls; ++c) {
        clause *cp = &(cls[c]);
        cp->weight = cp->weight * swt_p + scale_avg;
        if (cp->weight < 1) cp->weight = 1;
        delta_tot_cl_weight += cp->weight;
        if (delta_tot_cl_weight >= num_cls) {
            avg_cl_weight += 1;
            delta_tot_cl_weight -= num_cls;
        }
        if (0 == cp->sat_count) {
            for (lit l: cp->lits) {
                vars[l.var_num].score += cp->weight;
            }
        } else if (1 == cp->sat_count) {
            vars[cp->sat_var].score -= cp->weight;
        }
    }
}

void LSSolver::print_solution(bool need_verify) {
    if (0 == get_cost()) cout << "s SATISFIABLE" << endl;
    else cout << "s UNKNOWN" << endl;

    bool sat_flag = false;
    if (need_verify) {
        for (int c = 0; c < num_cls; c++) {
            sat_flag = false;
            for (lit l: cls[c].lits) {
                if (sol[l.var_num] == l.sense) {
                    sat_flag = true;
                    break;
                }
            }
            if (!sat_flag) {
                cout << "c Error: verify error in clause " << c << endl;
                exit(-1);
                return;
            }
        }
        cout << "c Verified." << endl;
    }

    if (verb > 0) {
        cout << "v";
        for (int v = 1; v <= num_vars; v++) {
            cout << ' ';
            if (sol[v] == 0)
                cout << '-';
            cout << v;
        }
        cout << endl;
    }
}

