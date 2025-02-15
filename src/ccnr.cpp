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
    max_tries = 1000LL*1000LL;
    max_steps = 1*1000 * 1000;
    random_gen.seed(1337);
    verb = 0;
}

bool LSSolver::make_space() {
    if (num_vars == 0) return false;
    vars.resize(num_vars+1);
    sol.resize(num_vars+1, 2);
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

bool LSSolver::local_search(long long int mems_limit , const char* prefix) {
    bool result = false;
    for (int t = 0; t < max_tries && !result; t++) {
        initialize();
        for (step = 0; step < max_steps && !result; step++) {
            if (unsat_cls.empty() && !touched_cls.empty()) {
                cout << "YAY" << endl;
                result = true;
                break;
            }
            if (touched_cls.empty()) {
              cout << prefix << "[ccnr] no touched cls, restart" << endl;
              break;
            }

            int flipv = pick_var();
            if (flipv == -1) {
              cout << prefix << "[ccnr] no var to flip, restart" << endl;
              break;
            }

            flip(flipv);
            if (mems > mems_limit) {
              cout << "mems limit reached" << endl;
              return false;
            }
            cout << "num unsat cls: " << unsat_cls.size() << " touched_cls: " << touched_cls.size() << endl;

            int u_cost = unsat_cls.size();
            int t_cost = touched_cls.size();
            if (verb && (step & 0x3ffff) == 0x3ffff) {
                cout << prefix << "[ccnr] tries: "
                << t << " steps: " << step
                << " unsat found: " << u_cost
                << " touched found: " << t_cost
                << endl;
            }
        }
    }
    return result;
}

void LSSolver::initialize() {
    unsat_cls.clear();
    unsat_vars.clear();
    touched_cls.clear();
    for (auto &i: idx_in_unsat_cls) i = 0;
    for (auto &i: idx_in_unsat_vars) i = 0;
    for (int v = 1; v <= num_vars; v++) {
      if (!indep_map[v]) {
        if (random_gen.next(100) < 20) {
          sol[v] = random_gen.next(2);
        } else {
          sol[v] = 2;
        }

        /* cout << "init var " << v << " to : " << (int) sol[v] << endl; */
      }
      else assert(sol[v] == 2);
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
            }
            if (val != 2) cl.touched_cnt++;
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

void LSSolver::print_cl(int cid) const {
    for(auto& l: cls[cid].lits) {
      cout << l << " ";
    }; cout << "0 " << " sat_cnt: " << cls[cid].sat_count << " touched_cnt: " << cls[cid].touched_cnt << endl;
}

int LSSolver::pick_var() {
    assert(!touched_cls.empty());
    update_clause_weights();
    uint32_t tries = 0;
    bool ok = false;
    int cid;
    while (!ok && tries < 100) {
      if (!unsat_cls.empty() && random_gen.next(100) < 20) {
        cid = unsat_cls[random_gen.next(unsat_cls.size())];
      } else {
        assert(!touched_cls.empty());
        cid = touched_cls[random_gen.next(touched_cls.size())];
      }
      assert(cid < (int)cls.size());

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
    /* cout << " ---------  " << endl; */
    /* cout << "decided on cl_id: " << cid << " -- "; print_cl(cid); */

    const clause& cl = cls[cid];
    int best_var = -1;
    int best_score = std::numeric_limits<int>::min();
    for (auto& l: cl.lits) {
        int v = l.var_num;
        if (indep_map[v]) {
          assert(sol[v] == 2);
          continue;
        }
        int score = vars[v].score;
        if (sol[v] == 2 && score > 0) score *= 0.8;

        if (score > best_score) {
            best_var = v;
            best_score = score;
        } else if (score == best_score &&
                   vars[v].last_flip_step < vars[best_var].last_flip_step) {
            best_var = v;
        }
    }
    assert(best_var != -1);
    /* cout << "decided on var: " << best_var << endl; */
    return best_var;
}

void LSSolver::check_clause(int cid) const {
  int sat_cnt = 0;
  int touched_cnt = 0;
  for (const lit& l: cls[cid].lits) {
    if (sol[l.var_num] == l.sense) sat_cnt++;
    if (sol[l.var_num] != 2) touched_cnt++;
  }
  cout << "Checking cl_id: " << cid << " -- "; print_cl(cid);
  cout << "sat_cnt: " << sat_cnt << " touched_cnt: " << touched_cnt << endl;
  cout << "cls[cid].sat_count: " << cls[cid].sat_count << " cls[cid].touched_cnt: " << cls[cid].touched_cnt << endl;
  assert(cls[cid].sat_count == sat_cnt);
  assert(cls[cid].touched_cnt == touched_cnt);
  if (touched_cnt > 0) {
    bool found = false;
    for(uint32_t i = 0; i < touched_cls.size(); i++) {
      int clid = touched_cls[i];
      if (clid == cid) {
        found = true;
        break;
      }
      assert(idx_in_touched_cls[clid] == (int)i);
    }
    assert(found);
  }
  if (sat_cnt == 0 && touched_cnt > 0) {
    bool found = false;
    for(int unsat_cl : unsat_cls) {
      if (unsat_cl == cid) {
        found = true;
        break;
      }
    }
    if (!found) { cout << "NOT found in unsat cls" << endl; }
    assert(idx_in_unsat_cls[cid] < (int)unsat_cls.size());
    assert(unsat_cls[idx_in_unsat_cls[cid]] == cid);
    assert(found);
  }
}

void LSSolver::check_unsat_cls() const {
  for(uint32_t i = 0; i < unsat_cls.size(); i++) {
    uint32_t clid = unsat_cls[i];
    assert(cls[clid].sat_count == 0);
    if (idx_in_unsat_cls[clid] != (int)i) {
      cout << "bad clid: " << clid << " i: " << i << endl;
    }
    assert(idx_in_unsat_cls[clid] == (int)i);
  }
}

void LSSolver::check_interals() const {
  for (uint32_t i = 0; i < cls.size(); i++) check_clause(i);
  for(uint32_t i = 0; i < unsat_cls.size(); i++) {
    uint32_t clid = unsat_cls[i];
    assert(idx_in_unsat_cls[clid] == (int)i);
    assert(cls[clid].sat_count == 0);
  }
}

void LSSolver::unset(int v) {
  assert(sol[v] != 2);
  SLOW_DEBUG_DO(check_interals());

  const int old_val = sol[v] ;
  sol[v] = 2;
  /* cout << "CHG unsetting var " << v << " prev val: " << (int)old_val << endl; */

  mems += vars[v].lits.size();

  // Go through each clause the literal is in and update status
  for (const lit& l: vars[v].lits) {
    clause& cl = cls[l.cl_num];
    assert(cl.sat_count >= 0);
    assert(cl.touched_cnt >= 0);
    assert(cl.sat_count <= (int)cl.lits.size());
    assert(cl.touched_cnt <= (int)cl.lits.size());

    cls[l.cl_num].touched_cnt--;
    assert(cl.touched_cnt >= 0);
    if (cls[l.cl_num].touched_cnt == 0) untouch_a_clause(l.cl_num);

    // it was previously satisfied
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
    } else {
      // it was previously unsatisfied
      if (cl.touched_cnt == 0) {
        // this clause is not fully untouched
        sat_a_clause(l.cl_num);
      }
    }

#ifdef SLOW_DEBUG
    cout << "Effect on cl_id: " << l.cl_num << " -- "; print_cl(l.cl_num);
    check_clause(l.cl_num);
    check_unsat_cls();
#endif
  }
}

void LSSolver::flip(int v) {
    assert(!indep_map[v]);
    SLOW_DEBUG_DO(check_interals());

    bool touch = false;
    const int old_val = sol[v] ;
    assert(old_val <= 2);
    if (sol[v] == 2) {
        // set to some value
        touch = true;
        sol[v] = random_gen.next(2);
        /* cout << "CHG setting var " << v << " new val: " << (int)sol[v] << endl; */
    } else {
        //flip
        sol[v] = 1 - sol[v];
        /* cout << "CHG flipping var " << v << " new val: " << (int)sol[v] << endl; */
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
        /* cout << "checking effect on cl_id: " << l.cl_num << " -- "; print_cl(l.cl_num); */

        if (touch) {
          cl.touched_cnt++;
          if (cl.touched_cnt == 1) touch_a_clause(l.cl_num);
        }

        if (sol[v] == l.sense) {
            // make it sat
            cl.sat_count++;

            if (cl.sat_count == 1) {
                if (old_val == 2 && cl.touched_cnt == 1) {
                  // first time touching, and made it sat, no need to update unsat_cls
                } else sat_a_clause(l.cl_num);
                cl.sat_var = v;
                for (const lit& lc: cl.lits) vars[lc.var_num].score -= cl.weight;
            } else if (cl.sat_count == 2) {
                vars[cl.sat_var].score += cl.weight;
            }
        } else if (sol[v] == !l.sense) {
          auto prev_sat_count = cl.sat_count;
          if (old_val != 2) {
            // make it unsat
            cl.sat_count--;
            assert(cl.sat_count >= 0);
          }

          // first time touching, and made it unsat
          if ((cl.sat_count == 0 && cl.touched_cnt == 1) ||
              // made it unsat with a flip
              (prev_sat_count > 0 && cl.sat_count == 0)
          ) {
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
        } else if (sol[v] == 2) {
          assert(false);
        }

#ifdef SLOW_DEBUG
        cout << "Effect on cl_id: " << l.cl_num << " -- "; print_cl(l.cl_num);
        check_clause(l.cl_num);
        for(uint32_t i = 0; i < unsat_cls.size(); i++) {
          uint32_t clid = unsat_cls[i];
          assert(cls[clid].sat_count == 0);
          if (idx_in_unsat_cls[clid] != (int)i) {
            cout << "bad clid: " << clid << " i: " << i << endl;
          }
          assert(idx_in_unsat_cls[clid] == (int)i);
        }
#endif
    }
    if (!touch) {
      vars[v].score = -orig_score;
    }
    vars[v].last_flip_step = step;

#ifdef SLOW_DEBUG
    cout << "Done flip(). Checking all clauses" << endl;
    for (uint32_t i = 0; i < cls.size(); i++) check_clause(i);
#endif
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
    if (index < (int)touched_cls.size()) {
      touched_cls[index] = last_item;
      idx_in_touched_cls[last_item] = index;
    }
}

void LSSolver::sat_a_clause(int cl_id) {
    /* cout << "sat_a_clause: cl_id: " << cl_id << endl; */
    assert(unsat_cls.size() > 0);
    if (cls[cl_id].touched_cnt > 0) assert(cls[cl_id].sat_count == 1);
    else assert(cls[cl_id].touched_cnt == 0 && cls[cl_id].sat_count == 0);

    //use the position of the clause to store the last unsat clause in stack
    int last_item = unsat_cls.back();
    unsat_cls.pop_back();
    int index = idx_in_unsat_cls[cl_id];
    if (index < (int)unsat_cls.size()) {
      assert(unsat_cls[index] == cl_id);
      unsat_cls[index] = last_item;
      idx_in_unsat_cls[last_item] = index;
    }

    //update unsat_appear and unsat_vars
    for (lit l: cls[cl_id].lits) {
        vars[l.var_num].unsat_appear--;
        if (0 == vars[l.var_num].unsat_appear) {
            last_item = unsat_vars.back();
            unsat_vars.pop_back();
            index = idx_in_unsat_vars[l.var_num];
            if (index < (int)unsat_vars.size()) {
              unsat_vars[index] = last_item;
              idx_in_unsat_vars[last_item] = index;
            }
        }
    }
}

void LSSolver::unsat_a_clause(int cl_id) {
    /* cout << "unsat_a_clause: cl_id: " << cl_id << endl; */
    assert(cls[cl_id].sat_count == 0);
    assert(cls[cl_id].touched_cnt > 0);
    unsat_cls.push_back(cl_id);
    idx_in_unsat_cls[cl_id] = unsat_cls.size()-1;
    assert(unsat_cls[idx_in_unsat_cls[cl_id]] == cl_id);

    //update unsat_appear and unsat_vars
    for (const lit& l: cls[cl_id].lits) {
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

    assert(need_verify);
    if (need_verify) {
        uint32_t num_cls_touched = 0;
        for (int cid = 0; cid < num_cls; cid++) {
            bool sat_flag = false;
            bool touched_flag = false;
            for (lit l: cls[cid].lits) {
                if (sol[l.var_num] != 2) touched_flag = true;
                if (sol[l.var_num] == l.sense) {
                    sat_flag = true;
                    break;
                }
            }
            num_cls_touched += touched_flag;
            if (!sat_flag && touched_flag) {
                cout << "c Error: verify error in cl_id : " << cid << " -- "; print_cl(cid);
                exit(-1);
                return;
            }
        }
        cout << "c Verified. Satisfied cls: " << num_cls_touched << endl;
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

