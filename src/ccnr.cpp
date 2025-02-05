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
#include "common.hpp"

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
    idx_in_touched_cls.resize(num_cls);

    return true;
}

// Updates variables' neighbor_vars, ran once
void LSSolver::build_neighborhood() {
    vector<uint8_t> flag(num_vars+1, 0);
    for (int vnum = 1; vnum <= num_vars; ++vnum) {
        variable& v = vars[vnum];
        for (const auto& l: v.lits) {
            int c_id = l.cl_num;
            for (const auto& lc: cls[c_id].lits) {
                if (!flag[lc.var_num] && lc.var_num != vnum) {
                    flag[lc.var_num] = 1;
                    v.neighbor_vars.push_back(lc.var_num);
                }
            }
        }
        for (const auto& v2 : v.neighbor_vars) flag[v2] = 0;
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
            if (flipv == -1) {
              cout << prefix << "[ccnr] no var to flip, restart" << endl;
              break;
            }

            flip(flipv);
            if (mems > _mems_limit) return false;
            cout << "num unsat cls: " << unsat_cls.size() << " num cls: " << cls.size() << endl;

            int u_cost = unsat_cls.size();
            int t_cost = touched_cls.size();
            if (verb && (step & 0x3ffff) == 0x3ffff) {
                cout << prefix << "[ccnr] tries: "
                << t << " steps: " << step
                << " unsat found: " << u_cost
                << " touched found: " << t_cost
                << endl;
            }

            if (u_cost == 0) {
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
      if (!indep_map[v]) {
        sol[v] = random_gen.next(3);
        cout << "init var " << v << " to : " << (int) sol[v] << endl;
      }
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
        if (cl.sat_count == 0 && cl.touched_cnt > 0) unsat_a_clause(cid);
        if (cl.touched_cnt > 0) touch_a_clause(cid);
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
        for (lit l: vp.lits) {
            int c = l.cl_num;
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

void LSSolver::print_cl(int cid) {
    for(auto& l: cls[cid].lits) {
      cout << l << " ";
    }; cout << "0 " << " sat_cnt: " << cls[cid].sat_count << " touched_cnt: " << cls[cid].touched_cnt << endl;
}

int LSSolver::pick_var() {
    assert(!unsat_cls.empty());
    update_clause_weights();
    uint32_t tries = 0;
    bool ok = false;
    int cid;
    while (!ok && tries < 100) {
      cid = unsat_cls[random_gen.next(unsat_cls.size())];
      const clause& cl = cls[cid];
      for (auto& l: cl.lits) {
        if (!indep_map[l.var_num]) {
          ok = true;
          break;
        }
      }
      tries++;
    }
    if (!ok) return -1;
    cout << " ---------  " << endl;
    cout << "decided on cl_id: " << cid << " -- "; print_cl(cid);

    const clause& cl = cls[cid];
    int best_var = -1;
    int best_score = std::numeric_limits<int>::min();
    for (auto& l: cl.lits) {
        int v = l.var_num;
        if (indep_map[v]) {
          assert(sol[v] == 3);
          continue;
        }
        int score = vars[v].score;
        if (sol[v] == 3) score /= 10;

        if (score > best_score) {
            best_var = v;
            best_score = score;
        } else if (score == best_score &&
                   vars[v].last_flip_step < vars[best_var].last_flip_step) {
            best_var = v;
        }
    }
    assert(best_var != -1);
    cout << "decided on var: " << best_var << endl;
    return best_var;
}

void LSSolver::flip(int v) {
    assert(!indep_map[v]);

    bool touch = false;
    int old_val = sol[v] ;
    if (sol[v] == 3) {
        // set to some value
        touch = true;
        sol[v] = random_gen.next(2);
        cout << "setting var " << v << " new val: " << (int)sol[v] << endl;
    } else if (random_gen.next(7) <= 5) {
        //flip
        sol[v] = 1 - sol[v];
        cout << "flipping var " << v << " new val: " << (int)sol[v] << endl;
    } else {
        // unset
        sol[v] = 3;
        cout << "unsetting var " << v << endl;
    }

    const int orig_score = vars[v].score;
    mems += vars[v].lits.size();

    // Go through each clause the literal is in and update status
    for (const lit& l: vars[v].lits) {
        clause& cl = cls[l.cl_num];
        assert(cl.sat_count >= 0);
        assert(cl.touched_cnt >= 0);
        assert(cl.sat_count <= (int)cl.lits.size());
        assert(cl.touched_cnt <= (int)cl.lits.size());
        cout << "checking effect on cl_id: " << l.cl_num << " -- "; print_cl(l.cl_num);

        if (touch) {
          cl.touched_cnt++;
          if (cl.touched_cnt == 1) touch_a_clause(l.cl_num);
        }

        if (sol[v] == l.sense) {
            // make it sat
            cl.sat_count++;
            cout << "Here, cnt: " << cl.sat_count << endl;

            if (cl.sat_count == 1) {
                sat_a_clause(l.cl_num);
                cl.sat_var = v;
                for (const lit& lc: cl.lits) vars[lc.var_num].score -= cl.weight;
            } else if (cl.sat_count == 2) {
                vars[cl.sat_var].score += cl.weight;
            }
        } else if (sol[v] == !l.sense) {
          if (old_val != 3) {
            // make it unsat
            cl.sat_count--;
            assert(cl.sat_count >= 0);
          }

          if (cl.sat_count == 0 && cl.touched_cnt > 0) {
              unsat_a_clause(l.cl_num);
              for (const lit& lc: cl.lits) vars[lc.var_num].score += cl.weight;
          } else if (cl.sat_count == 1) {
              // Have to update the var that makes the clause satisfied
              for (const lit& lc: cl.lits) {
                  if (sol[lc.var_num] == lc.sense) {
                      vars[lc.var_num].score -= cl.weight;
                      cl.sat_var = lc.var_num;
                      break;
                  }
              }
          }
        } else if (sol[v] == 3) {
          // unset
          cls[l.cl_num].touched_cnt--;
          assert(cl.touched_cnt >= 0);
          if (cls[l.cl_num].touched_cnt == 0) untouch_a_clause(l.cl_num);

          if (old_val == l.sense) {
            cl.sat_count--;
            assert(cl.sat_count >= 0);

            // make it unsat
            if (cl.sat_count == 0 && cl.touched_cnt > 0) {
              unsat_a_clause(l.cl_num);
              for (const lit& lc: cl.lits) vars[lc.var_num].score += cl.weight;
            } else if (cl.sat_count == 1) {
              // Have to update the var that makes the clause satisfied
              for (const lit& lc: cl.lits) {
                if (sol[lc.var_num] == lc.sense) {
                  vars[lc.var_num].score -= cl.weight;
                  cl.sat_var = lc.var_num;
                  break;
                }
              }
            }
          }
        }
    }
    if (!touch) {
      vars[v].score = -orig_score;
    }
    vars[v].last_flip_step = step;
}

void LSSolver::touch_a_clause(int cl_id) {
  assert(cls[cl_id].touched_cnt > 0);
  touched_cls.push_back(cl_id);
  idx_in_touched_cls[cl_id] = touched_cls.size()-1;
}

void LSSolver::untouch_a_clause(int cl_id) {
    assert(!touched_cls.empty());
    assert(cls[cl_id].touched_cnt == 0);

    int last_item = touched_cls.back();
    touched_cls.pop_back();
    int index = idx_in_touched_cls[cl_id];
    if (index < (int)touched_cls.size()) touched_cls[index] = last_item;
    idx_in_touched_cls[last_item] = index;
}

void LSSolver::sat_a_clause(int cl_id) {
    //use the position of the clause to store the last unsat clause in stack
    int last_item = unsat_cls.back();
    unsat_cls.pop_back();
    int index = idx_in_unsat_cls[cl_id];
    if (index < (int)unsat_cls.size()) unsat_cls[index] = last_item;
    idx_in_unsat_cls[last_item] = index;

    //update unsat_appear and unsat_vars
    for (lit l: cls[cl_id].lits) {
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

