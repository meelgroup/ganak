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

#pragma once

#include <cstdint>
#include <map>

#include "clauseallocator.hpp"
#include "common.hpp"
#include "counter_config.hpp"
#include "primitive_types.hpp"
#include "statistics.hpp"
#include "comp_management.hpp"
#include "boundedqueue.hpp"
#include "TreeDecomposition.hpp"
#include "structures.hpp"
#include "heap.hpp"
#include "mpreal.h"

#include <cryptominisat5/cryptominisat.h>

using std::pair;
using std::map;


template<typename T>
inline vector<CMSat::Lit> ganak_to_cms_cl(const T& cl) {
  vector<CMSat::Lit> cms_cl;
  cms_cl.reserve(cl.size());
  for(const auto& l: cl) cms_cl.push_back(CMSat::Lit(l.var()-1, !l.sign()));
  return cms_cl;
}

inline vector<CMSat::Lit> ganak_to_cms_cl(const Lit& l) {
  vector<CMSat::Lit> cms_cl;
  cms_cl.push_back(CMSat::Lit(l.var()-1, !l.sign()));
  return cms_cl;
}

enum class RetState {
  EXIT,
  RESOLVED,
  PROCESS_COMPONENT,
  BACKTRACK,
  GO_AGAIN
};

using enum RetState;

inline std::ostream& operator<<(std::ostream& os, const RetState& val) {
  std::stringstream s;
  switch (val) {
    case RetState::EXIT : os << "EXIT"; break;
    case RetState::RESOLVED: os << "RESOLVED"; break;
    case RetState::PROCESS_COMPONENT: os << "PROCESS_COMPONENT"; break;
    case RetState::BACKTRACK : os << "BACKTRACK"; break;
    case RetState::GO_AGAIN : os << "GO_AGAIN"; break;
  }
  return os;
}


struct VS {
  VS() = default;
  VS(uint32_t _v, double _score1, uint32_t _score2) : v(_v), score1(_score1), score2(_score2) {}
  bool operator<(const VS& other) const {
    if (score1 != other.score1) return score1 > other.score1;
    else return score2 > other.score2;
  }
  uint32_t v;
  double score1 = 0;
  uint32_t score2 = 0;
};


struct OffAbs {
  OffAbs(ClauseOfs _off, uint32_t _abs) : off(_off), abs(_abs) {}
  ClauseOfs off;
  uint32_t abs;
};

struct BinClSub {
  BinClSub() = default;
  BinClSub (Lit _lit1, Lit _lit2, bool _red) :
    red(_red) {
      lit[0] = _lit1;
      lit[1] = _lit2;
    }

  uint32_t size() const { return 2; }
  const Lit* begin() const {return &lit[0];}
  Lit* begin() {return &lit[0];}
  const Lit* end() const {return begin()+2;}
  Lit* end() {return begin()+2;}
  Lit& operator[](const uint32_t at) {return lit[at];}
  const Lit& operator[](const uint32_t at) const {return lit[at];}
  bool operator<(const BinClSub& other) const {
    if (lit[0] != other.lit[0]) return lit[0]<other.lit[0];
    if (lit[1] != other.lit[1]) return lit[1]<other.lit[1];
    if (red != other.red) return !red;
    return false;
  }
  bool operator==(const BinClSub& other) const {
    return lit[0] == other.lit[0] && lit[1] == other.lit[1] &&
      red == other.red;
  }

  Lit lit[2];
  bool red;
};

inline std::ostream& operator<<(std::ostream& os, const BinClSub& cl) {
  for(const auto& l: cl) os << l << " ";
  os << "0" << " (red: " << (int)cl.red << ")";
  return os;
}

struct VarOrderLt {
  const LiteralIndexedVector<LitWatchList>& watches;
  bool operator () (uint32_t x, uint32_t y) const {
    Lit l1 = Lit(x, 0);
    auto act1 = watches[l1].activity + watches[l1.neg()].activity;
    Lit l2 = Lit(y, 0);
    auto act2 = watches[l2].activity + watches[l2.neg()].activity;
    return act1 > act2;
  }
  VarOrderLt(const LiteralIndexedVector<LitWatchList>& _watches) : watches(_watches) { }
};

template<typename T> class ClauseAllocator;

template<typename T>
class Counter {
public:
  Counter(const CounterConfiguration& _conf);
  ~Counter();
  const CounterConfiguration& get_conf() const { return conf;}
  friend class ClauseAllocator<T>;
  struct ConflictData {
    int32_t nHighestLevel = -1;
    bool bOnlyOneLitFromHighest = false;
  };
  ConflictData find_conflict_level(Lit p);
  void new_vars(const uint32_t n);
  uint32_t get_num_low_lbds() const { return num_low_lbd_cls; }
  uint32_t get_num_long_reds() const { return long_red_cls.size(); }
  uint32_t get_num_irred_long_cls() const { return long_irred_cls.size(); }
  int val(Lit lit) const { return values[lit]; }
  int val(uint32_t var) const { return values[Lit(var,1)]; }
  bool get_is_approximate() const { return is_approximate; }

  friend class ClauseAllocator<T>;
  ClauseAllocator<T>* alloc;
  vector<ClauseOfs> long_irred_cls;
  vector<ClauseOfs> long_red_cls;
  uint32_t nVars() const { return var_data.size() - 1; }
  double get_start_time() const { return start_time;}
  const auto& get_cache() const { return comp_manager->get_cache();}
  void set_generators(const vector<map<Lit, Lit>>& _gens) { generators = _gens; }
  void end_irred_cls();
  void set_indep_support(const set<uint32_t>& indeps);
  T outer_count();
  bool add_red_cl(const vector<Lit>& lits, int lbd = -1);
  bool add_irred_cl(const vector<Lit>& lits);
  void set_optional_indep_support(const set<uint32_t>& indeps);
  int32_t decision_level() const { return decisions.get_decision_level();}
  void set_lit_weight(Lit l, const T& w);
  const T& get_weight(const Lit l) { return weights[l.raw()];}
  T get_weight(const uint32_t v) {
    Lit l(v, false);
    return weights[l.raw()]+weights[l.neg().raw()];}

  // vivif stuff
  void vivif_setup();
  bool v_propagate();
  void v_backtrack();
  void v_unset_lit(const Lit l);
  void v_enqueue(const Lit l);
  TriValue v_val(const Lit l) const;
  void v_new_lev();
  void v_backup();
  void v_restore();
protected:
  CounterConfiguration conf;
  void unset_lit(Lit lit) {
    VERBOSE_DEBUG_DO(cout << "Unsetting lit: " << std::setw(8) << lit << endl);
    SLOW_DEBUG_DO(assert(val(lit) == T_TRI));
    var(lit).ante = Antecedent();
    if constexpr (weighted) if(!sat_mode() && get_weight(lit) != 1) {
      uint64_t* at = vars_act_dec.data()+decision_level()*(nVars()+1);
      bool found = (at[0] == at[lit.var()]);
      if (found) decisions[decision_level()].include_solution(get_weight(lit));
    }
    var(lit).decision_level = INVALID_DL;
    values[lit] = X_TRI;
    values[lit.neg()] = X_TRI;
  }

  const Antecedent & get_antec(Lit lit) const {
    return var_data[lit.var()].ante;
  }

  bool is_antec_of(ClauseOfs ante_cl, Lit lit) const {
    return var(lit).ante.isAClause() && (var(lit).ante.asCl() == ante_cl);
  }

  void reduce_db();
  template<class T2> void minimize_uip_cl_with_bins(T2& cl);
  vector<Lit> tmp_minim_with_bins;
  void delete_cl(const ClauseOfs cl_ofs);
  bool red_cl_can_be_deleted(ClauseOfs cl_ofs);

  // Super-slow debug, not even used right now
  bool find_offs_in_watch(const vector<ClOffsBlckL>& ws, ClauseOfs off) const;
  void check_all_cl_in_watchlists() const;

  // the first variable that is NOT in the independent support
  uint32_t indep_support_end = std::numeric_limits<uint32_t>::max();
  uint32_t opt_indep_support_end = std::numeric_limits<uint32_t>::max();

  LiteralIndexedVector<LitWatchList> watches;
  vector<T> weights;
  vector<Lit> unit_clauses_;
  vector<VarData> var_data;
  bool num_vars_set = false;
  LiteralIndexedVector<TriValue> values;
  vector<double> tdscore;
  double td_weight = 1.0;
  int td_width = 10000;
  uint32_t lbd_cutoff;
  uint32_t num_low_lbd_cls = 0; // Last time counted low LBD clauses
  uint32_t num_used_cls = 0; // last time counted used clauses
  DataAndStatistics<T> stats;

  // Computing LBD (lbd == 2 means "glue clause")
  vector<uint64_t> lbdHelper;
  uint64_t lbdHelperFlag = 0;

  template<class T2>
  uint32_t calc_lbd(const T2& lits) {
    lbdHelperFlag++;
    uint32_t nblevels = 0;
    for(const auto& l: lits) {
      if (val(l) == X_TRI) {nblevels++;continue;}
      int lev = var(l).decision_level;
      if (lev != 0 && lbdHelper[lev] != lbdHelperFlag) {
        lbdHelper[lev] = lbdHelperFlag;
        nblevels++;
        if (nblevels >= 100) { return nblevels; }
      }
    }
    return nblevels;
  }

  bool exists_unit_cl_of(const Lit l) const {
    for (const auto& l2 : unit_clauses_) if (l == l2) return true;
    return false;
  }

  template<typename T2> void attach_cl(ClauseOfs off, const T2& lits);
  Clause* add_cl(const vector<Lit> &literals, bool red);

  // adds a UIP Conflict Clause
  // and returns it as an Antecedent to the first
  // literal stored in literals
  inline Antecedent add_uip_confl_cl(const vector<Lit> &literals);

  inline bool add_bin_cl(Lit a, Lit b, bool red);

  inline VarData &var(const Lit lit) {
    return var_data[lit.var()];
  }

  inline VarData &var(const uint32_t v) {
    return var_data[v];
  }

  inline const VarData &var(const uint32_t v) const{
    return var_data[v];
  }

  inline const VarData &var(const Lit lit) const {
    return var_data[lit.var()];
  }

  inline bool is_true(const Lit &lit) const {
    return values[lit] == T_TRI;
  }

  bool is_false(Lit lit) {
    return values[lit] == F_TRI;
  }

  string lit_val_str(Lit lit) const {
    if (values[lit] == F_TRI) return "FALSE";
    else if (values[lit] == T_TRI) return "TRUE";
    else return "UNKN";
  }

  string val_to_str(const TriValue& tri) const {
    if (tri == F_TRI) return "FALSE";
    else if (tri == T_TRI) return "TRUE";
    else return "UNKN";
  }

  bool is_unknown(Lit lit) const {
    SLOW_DEBUG_DO(assert(lit.var() <= nVars()));
    SLOW_DEBUG_DO(assert(lit.var() != 0));
    return values[lit] == X_TRI;
  }

  bool is_unknown(uint32_t var) const {
    SLOW_DEBUG_DO(assert(var != 0));
    return is_unknown(Lit(var, false));
  }

  bool is_satisfied(ClauseOfs off) {
    for (auto lt: *alloc->ptr(off)) if (is_true(lt)) return true;
    return false;
  }
  bool counted_bottom_comp = true; //when false, we MUST take suggested polarities
  vector<uint8_t> seen;

  double score_of(const uint32_t v, bool ignore_td = false) const;
  double var_act(const uint32_t v) const;

  // DNF Cube stuff
  vector<Cube<T>> mini_cubes;
  void disable_small_cubes(vector<Cube<T>>& cubes);
  void disable_smaller_cube_if_overlap(uint32_t i, uint32_t i2, vector<Cube<T>>& cubes);
  void print_and_check_cubes(vector<Cube<T>>& cubes);
  void disable_cubes_if_overlap(vector<Cube<T>>& cubes);
  void extend_cubes(vector<Cube<T>>& cubes);
  int cube_try_extend_by_lit(const Lit torem, const Cube<T>& c);
  T check_count_norestart(const Cube<T>& c);
  T check_count_norestart_cms(const Cube<T>& c);
  vector<Cube<T>> one_restart_count();
  bool clash_cubes(const set<Lit>& c1, const set<Lit>& c2) const;
  bool compute_cube(Cube<T>& cube, int branch);
  vector<map<Lit, Lit>> generators;
  void symm_cubes(vector<Cube<T>>& cubes);

  vector<uint32_t> common_indep_code(const set<uint32_t>& indeps);
  const DataAndStatistics<T>& get_stats() const;
  void get_unit_cls(vector<Lit>& units) const;
  int32_t dec_level() const { return decisions.get_decision_level(); }
  void fill_cl(const Antecedent& ante, Lit*& c, uint32_t& size, Lit p) const;

private:
  T check_count(const bool also_incl_curr_and_later_dec = false);
  static constexpr bool weighted = std::is_same<T, mpfr::mpreal>::value || std::is_same<T, mpq_class>::value;
  void init_activity_scores();
  vector<vector<Lit>> v_backup_cls;
  vector<vector<ClOffsBlckL>> v_backup_watches;
#ifdef SLOW_DEBUG
  vector<vector<Lit>> debug_irred_cls;
#endif
  bool remove_duplicates(vector<Lit>& lits);
  CMSat::SATSolver* sat_solver = nullptr;
  bool ok = true;
  bool isindependent = true;

  // Needed to know what variables were active in given decision levels
  // It's needed for weighted counting to know what variable was active in
  // that component (learnt clauses can propagate vars that were not active,
  // and on backtrack, we'd multiply by them, wrongly)
  vector<uint64_t> vars_act_dec;
  uint64_t vars_act_dec_num = 0;


  // Temporaries, used during recordLastUIPClause
  mutable vector<Lit> tmp_lit; //used in recoredLastUIPClause
  vector<uint32_t> to_clear;

  vector<Lit> tmp_cl_minim; // Used during minimize_uip_cl

  double start_time;
  std::mt19937_64 mtrand;
  volatile bool appmc_timeout_fired = false;
  bool is_approximate = false;

  DecisionStack<T> decisions;
  vector<Lit> trail;
  uint32_t qhead = 0;
  CompManager<T>* comp_manager = nullptr;
  uint64_t last_reducedb_confl = 0;
  uint64_t last_reducedb_dec = 0;

  void simple_preprocess();
  bool is_implied(const vector<Lit>& cp);
  void check_implied(const vector<Lit>& cl);

  void count_loop();
  bool decide_lit();
  uint32_t find_best_branch(bool ignore_td = false);
  template<class T2> bool clause_falsified(const T2& cl) const;
  bool clause_asserting(const vector<Lit>& cl) const;
  template<class T2> bool clause_satisfied(const T2& cl) const;
  void compute_score(TWD::TreeDecomposition& tdec, bool print = true);
  void td_decompose();
  TWD::TreeDecomposition td_decompose_component(double mult = 1);
  double td_lookahead_score(const uint32_t v, const uint32_t base_comp_tw);
  void recomp_td_weight();
  void vsads_readjust();

  // Actual SAT solver.
  bool use_sat_solver(RetState& state);
  int32_t sat_start_dec_level = -1;
  Heap<VarOrderLt> order_heap;
  inline bool sat_mode() const {
    return sat_start_dec_level != -1 && decision_level() >= sat_start_dec_level;
  }
  vector<int> sat_solution;
  void check_sat_solution() const;
  // this is the actual BCP algorithm
  // starts propagating all literal in trail_
  // beginning at offset start_at_trail_ofs
  bool propagate(bool out_of_order = false);
  void get_maxlev_maxind(ClauseOfs ofs, int32_t& maxlev, uint32_t& maxind);
  bool check_watchlists() const;
  template<class T2> void check_cl_propagated_conflicted(T2& cl, uint32_t off = 0) const;
  void check_all_propagated_conflicted() const;

  void print_all_levels();
  bool restart_if_needed();
  RetState backtrack();
  void print_dec_info() const;
  template<class T2> void print_cl(const T2& cl) const;
  template<class T2> void v_print_cl(const T2& cl) const;
  void print_cl(const Lit* c, uint32_t size) const;
  void check_cl_unsat(Lit* c, uint32_t size) const;
  void print_conflict_info() const;
  void print_comp_stack_info() const;

  // AppMC
  mpz_class do_appmc_count();

  // BDD
  bool do_buddy_count(const Comp* c);
  uint64_t buddy_count();
  vector<uint32_t> vmap;
  vector<uint32_t> vmap_rev;

  // if on the current decision level
  // a second branch can be visited, RESOLVED is returned
  // otherwise returns BACKTRACK
  RetState resolve_conflict();
  void go_back_to(int32_t backj);
  uint32_t find_lev_to_set(int32_t implied_lit_lev);
  size_t find_backtrack_level_of_learnt();
  void print_trail(bool check_entail = true, bool check_anything = true) const;
  void check_trail(bool check_entail = true) const;
  bool chrono_work();
  void reduce_db_if_needed();
  void inc_act(const Lit lit);
  void set_lit(const Lit lit, int32_t dec_lev, Antecedent ant = Antecedent());

  void set_confl_state(Lit a, Lit b) {
    confl_lit = a;
    confl = Antecedent(b);
  }

  void set_confl_state(Clause* cl) {
    if (cl->red && cl->lbd > this->lbd_cutoff) {
      cl->set_used();
      /* cl->update_lbd(this->calc_lbd(*cl)); */
    }
    confl = Antecedent(this->alloc->get_offset(cl));
    confl_lit = NOT_A_LIT;
  }

  // The literals that have been set in this decision level
  vector<Lit>::const_iterator top_declevel_trail_begin() const {
    return trail.begin() + this->var(decisions.top().var).sublevel;
  }
  vector<Lit>::iterator top_declevel_trail_begin() {
    return trail.begin() + this->var(decisions.top().var).sublevel;
  }

  void init_decision_stack() {
    decisions.clear();
    trail.clear();
    // initialize the stack to contain at least level zero
    decisions.push_back(StackLevel<T>(
          1, // super comp
          2)); //comp stack offset

    // I guess this is needed so the system later knows it's fully counted
    // since this is only a dummy.
    decisions.back().change_to_right_branch();
  }

  const Lit &top_dec_lit() const { return *top_declevel_trail_begin(); }
  uint32_t trail_at_dl(uint32_t dl) const { return this->var_data[decisions.at(dl).var].sublevel; }
  uint32_t trail_at_top() const { return this->var_data[decisions.top().var].sublevel; }
  void reactivate_comps_and_backtrack_trail([[maybe_unused]] bool check_ws = true);


  /////////////////////////////////////////////
  //  Conflict analysis below
  /////////////////////////////////////////////

  // if the state name is CONFLICT,
  // then violated_clause contains the clause determining the conflict;
  Lit confl_lit = NOT_A_LIT;
  Antecedent confl;
  vector<Lit> uip_clause;

  void create_uip_cl();
  void minimize_uip_cl();
  uint32_t abst_level(const uint32_t x) const;
  bool lit_redundant(Lit p, uint32_t abstract_levels);
  vector<Lit> analyze_stack;
  void recursive_cc_min();
  bool get_polarity(const uint32_t var) const;
  bool standard_polarity(const uint32_t var) const;

  // Vivification
  int64_t v_tout;
  vector<Lit> v_tmp;
  vector<Lit> v_tmp2;
  vector<Lit> v_cl;
  uint64_t last_confl_vivif = 0;
  struct SavedCl {
    SavedCl (Lit _first, Lit _second, bool _currently_propagating) :
      first(_first), second(_second), currently_propagating(_currently_propagating)
    {
      blk1 = NOT_A_LIT;
      blk2 = NOT_A_LIT;
    }
    SavedCl () = default;

    Lit first;
    Lit second;
    Lit blk1;
    Lit blk2;
    bool currently_propagating = false;
  };
  map<ClauseOfs, SavedCl> off_to_lit12;
  void v_cl_toplevel_repair(vector<ClauseOfs>& offs);
  void v_cl_repair(ClauseOfs off);
  void vivify_cls(vector<ClauseOfs>& cls);
  void vivify_all(bool force = false, bool only_irred = false);
  bool vivify_cl(const ClauseOfs off);
  void v_shrink(Clause& c);
  template<class T2> bool v_unsat(const T2& lits);
  template<class T2> bool v_satisfied(const T2& lits);
  void v_fix_watch(Clause& cl, uint32_t i);
  template<class T2> bool propagating_cl(T2& cl) const;
  template<class T2> bool conflicting_cl(T2& cl) const;
  template<class T2> bool propagation_correctness_of_vivified(const T2& cl) const;
  template<class T2> bool currently_propagating_cl(T2& cl) const;
  template<class T2> bool v_cl_satisfied(const T2& cl) const;
  void vivif_backtrack();
  vector<Lit> v_trail;
  uint32_t v_qhead;
  uint32_t v_lev;
  vector<int32_t> v_levs;
  uint32_t v_backtrack_to;
  LiteralIndexedVector<TriValue> v_values;

  // Toplevel stuff
  void subsume_all();
  void attach_occ(vector<ClauseOfs>& offs, bool sort_and_clear);
  inline uint32_t abst_var(const uint32_t v) {return 1UL << (v % 29);}
  template <class T2> uint32_t calc_abstr(const T2& ps) {
    uint32_t abs = 0;
    if (ps.size() > 50) return ~((uint32_t)(0ULL));
    for (auto l: ps) abs |= abst_var(l.var());
    return abs;
  }
  inline bool subset_abstr(const uint32_t a, const uint32_t b) { return ((a & ~b) == 0); }

  template<class T1, class T2> bool subset(const T1& a, const T2& b);
  vector<vector<OffAbs>> occ;
  vector<ClauseOfs> clauses;
  void backw_susume_cl(ClauseOfs off);
  void backw_susume_cl_with_bin(BinClSub& b);
  void toplevel_full_probe();
  vector<Lit> bothprop_toset;

  void print_stat_line();
  uint64_t next_print_stat_cache = 4ULL*1000LL*1000LL;
  uint64_t next_print_stat_confl = 100LL*1000LL;

  // indicates if we have called end_irred_cls()
  bool ended_irred_cls = false;
};

template<typename T>
inline void Counter<T>::print_cl(const Lit* c, uint32_t size) const {
  for(uint32_t i = 0; i < size; i++) {
    Lit l = c[i];
    cout << std::setw(5) << l
      << " lev: " << std::setw(3) << var(l).decision_level
      << " ante: " << std::setw(8) << var(l).ante
      << " val: " << std::setw(7) << lit_val_str(l)
      << " sublev: " << std::setw(3) << var(l).sublevel
      << endl;
  }
}

template<typename T>
template<class T2>
void Counter<T>::print_cl(const T2& cl) const {
  for(uint32_t i = 0; i < cl.size(); i ++) {
    const auto l = cl[i];
    cout << std::left << std::setw(5) << l
      << " lev: " << std::setw(4) << var(l).decision_level
      << " ante: " << std::setw(5) << std::left << var(l).ante
      << " val: " << lit_val_str(l)
      << " sublev: " << std::setw(3) << var(l).sublevel
      << endl;
  }
}

template<typename T>
template<class T2>
void Counter<T>::v_print_cl(const T2& cl) const {
  for(uint32_t i = 0; i < cl.size(); i ++) {
    const auto l = cl[i];
    cout << std::setw(5) << l
      << " lev: " << std::setw(4) << v_levs[l.var()]
      << " val: " << val_to_str(v_val(l)) << endl;
  }
}

template<typename T>
template<class T2> bool Counter<T>::conflicting_cl(T2& cl) const {
  for(const auto&l: cl) {
    if (val(l) == T_TRI || val(l) == X_TRI) return false;
  }
  return true;
}

template<typename T>
template<class T2> bool Counter<T>::propagating_cl(T2& cl) const {
  uint32_t unk = 0;
  for(const auto&l: cl) {
    if (val(l) == T_TRI) return false;
    if (val(l) == X_TRI) {unk++; if (unk>1) break;}
  }
  return unk == 1;
}

template<typename T>
template<class T2> bool Counter<T>::currently_propagating_cl(T2& cl) const {
  uint32_t tru = 0;
  for(const auto&l: cl) {
    if (val(l) == T_TRI) {tru++; if (tru>1) return false;}
    if (val(l) == X_TRI) return false;
  }
  return tru == 1;
}

template<typename T>
inline void Counter<T>::check_cl_unsat(Lit* c, uint32_t size) const {
  bool all_false = true;
  for(uint32_t i = 0; i < size; i++) {
    if (val(c[i]) != F_TRI) {all_false = false; break;}
  }
  if (all_false) return;

  cout << "ERROR: clause is not falsified." << endl;
  print_cl(c, size);
  assert(false);
}

// this is ONLY entered, if seen[lit.var()] is false, hence this is ALWAYS a single bump
// to each variable during analysis
template<typename T>
void inline Counter<T>::inc_act(const Lit lit) {
  watches[lit].activity += 1.0;
  if (sat_mode() && order_heap.inHeap(lit.var())) order_heap.increase(lit.var());
}

template<typename T>
template<class T1, class T2> bool Counter<T>::subset(const T1& a, const T2& b) {
#ifdef VERBOSE_DEBUG
  cout << "A:" << a << endl;
  for(size_t i = 1; i < a.size(); i++) assert(a[i-1] < a[i]);
  cout << "B:" << b << endl;
  for(size_t i = 1; i < b.size(); i++) assert(b[i-1] < b[i]);
#endif
  bool ret;
  uint32_t i = 0;
  uint32_t i2;
  Lit last_b = NOT_A_LIT;
  for (i2 = 0; i2 < b.size(); i2++) {
    if (last_b != NOT_A_LIT) assert(last_b < b[i2]);
    last_b = b[i2];
    //Literals are ordered
    if (a[i] < b[i2]) {
        ret = false;
        goto end;
    }
    else if (a[i] == b[i2]) {
      i++;
      //went through the whole of A now, so A subsumes B
      if (i == a.size()) {
          ret = true;
          goto end;
      }
    }
  }
  ret = false;

  end:
  return ret;
}

template<typename T>
Antecedent Counter<T>::add_uip_confl_cl(const vector<Lit> &literals) {
  Antecedent ante;
  stats.num_clauses_learned_++;
  Clause* cl = add_cl(literals, true);
  if (cl) {
    auto off = alloc->get_offset(cl);
    long_red_cls.push_back(off);
    cl->lbd = calc_lbd(*cl);
    if (cl->lbd <= (lbd_cutoff+1)) {
      for(const auto& l: literals) inc_act(l);
    }
    ante = Antecedent(off);
  } else if (literals.size() == 2){
    ante = Antecedent(literals.back());
    stats.num_binary_red_clauses_++;
  }
  return ante;
}

template<typename T>
bool Counter<T>::add_bin_cl(Lit a, Lit b, bool red) {
   watches[a].add_bin(b, red);
   watches[b].add_bin(a, red);
   return true;
}

template<typename T>
template<class T2>
void Counter<T>::minimize_uip_cl_with_bins(T2& cl) {
  SLOW_DEBUG_DO(for(const auto& s: seen) assert(s == 0););
  uint32_t rem = 0;
  assert(cl.size() > 0);
  tmp_minim_with_bins.clear();
  for(const auto& l: cl) { seen[l.raw()] = 1; tmp_minim_with_bins.push_back(l);}
  for(const auto& l: cl) {
  /* { */
    /* Lit l = tmp_minim_with_bins[0]; */
    if (!seen[l.raw()]) continue;
    const auto& w = watches[l].binaries;
    for(const auto& bincl: w) {
      const auto& l2 = bincl.lit();
      assert(l.var() != l2.var());
      if (seen[(l2.neg()).raw()]) { seen[(l2.neg()).raw()] = 0; rem++; }
    }
  }
  cl.clear(); cl.push_back(tmp_minim_with_bins[0]);
  seen[tmp_minim_with_bins[0].raw()] = 0;
  for(uint32_t i = 1; i < tmp_minim_with_bins.size(); i++) {
    Lit l = tmp_minim_with_bins[i];
    if (seen[l.raw()]) {
      cl.push_back(l);
      seen[l.raw()] = 0;
    }
  }
  stats.rem_lits_with_bins+=rem;
  stats.rem_lits_tried++;
}

template<typename T>
template<class T2> void Counter<T>::attach_cl(ClauseOfs off, const T2& lits) {
  Lit blck_lit = lits[lits.size()/2];
  watches[lits[0]].add_cl(off, blck_lit);
  watches[lits[1]].add_cl(off, blck_lit);
}

class OuterCounter {
public:
  OuterCounter(const CounterConfiguration& conf, bool weighted, bool precise = false) {
    if (!weighted) unw_counter = new Counter<mpz_class>(conf);
    else {
      if (!precise) w_counter = new Counter<mpfr::mpreal>(conf);
      else wq_counter = new Counter<mpq_class>(conf);
    }
  }
  ~OuterCounter() {
    delete unw_counter;
    delete w_counter;
    delete wq_counter;
  }
  void set_generators(const vector<map<Lit, Lit>>& _gens) {
    if (unw_counter) unw_counter->set_generators(_gens);
    if (w_counter) w_counter->set_generators(_gens);
    if (wq_counter) wq_counter->set_generators(_gens);
  }
  void end_irred_cls() {
    if (unw_counter) unw_counter->end_irred_cls();
    if (w_counter) w_counter->end_irred_cls();
    if (wq_counter) wq_counter->end_irred_cls();
  }
  void set_indep_support(const set<uint32_t>& indeps) {
    if (unw_counter) unw_counter->set_indep_support(indeps);
    if (w_counter) w_counter->set_indep_support(indeps);
    if (wq_counter) wq_counter->set_indep_support(indeps);
  }
  mpz_class unw_outer_count() { release_assert(unw_counter); return unw_counter->outer_count();}
  mpfr::mpreal w_outer_count() { release_assert(w_counter); return w_counter->outer_count();}
  mpq_class wq_outer_count() { release_assert(wq_counter); return wq_counter->outer_count();}
  bool add_red_cl(const vector<Lit>& lits, int lbd = -1) {
    if (unw_counter) return unw_counter->add_red_cl(lits, lbd);
    if (w_counter) return w_counter->add_red_cl(lits, lbd);
    if (wq_counter) return wq_counter->add_red_cl(lits, lbd);
    release_assert(false);
  }
  bool get_is_approximate() const { release_assert(unw_counter); return unw_counter->get_is_approximate();}
  bool add_irred_cl(const vector<Lit>& lits) {
    if (unw_counter) return unw_counter->add_irred_cl(lits);
    if (w_counter) return w_counter->add_irred_cl(lits);
    if (wq_counter) return wq_counter->add_irred_cl(lits);
    release_assert(false);
  }
  void set_optional_indep_support(const set<uint32_t>& indeps) {
    if (unw_counter) unw_counter->set_optional_indep_support(indeps);
    if (w_counter) w_counter->set_optional_indep_support(indeps);
    if (wq_counter) wq_counter->set_optional_indep_support(indeps);
  }
  void set_lit_weight(const Lit l, const mpq_class& w) {
    release_assert(w_counter || wq_counter);
    if (w_counter) w_counter->set_lit_weight(l, w.get_mpq_t());
    if (wq_counter) wq_counter->set_lit_weight(l, w);
  }
  void new_vars(const uint32_t n) {
    if (unw_counter) unw_counter->new_vars(n);
    if (w_counter) w_counter->new_vars(n);
    if (wq_counter) wq_counter->new_vars(n);
  }
private:
  Counter<mpz_class>* unw_counter = nullptr;
  Counter<mpfr::mpreal>* w_counter = nullptr;
  Counter<mpq_class>* wq_counter = nullptr;
};
