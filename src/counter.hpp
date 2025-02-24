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
#include "comp_management.hpp"
#include "statistics.hpp"
#include "comp_management.hpp"
#include "TreeDecomposition.hpp"
#include "structures.hpp"
#include "heap.hpp"
#ifdef BUDDY_ENABLED
#include "bdd.h"
#endif

#include <cryptominisat5/cryptominisat.h>

using std::pair;
using std::map;

namespace GanakInt {

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
  void new_vars(const uint32_t n);
  void set_indep_support(const set<uint32_t>& indeps);
  void set_optional_indep_support(const set<uint32_t>& indeps);
  void print_indep_distrib() const;
  bool add_irred_cl(const vector<Lit>& lits);
  bool add_red_cl(const vector<Lit>& lits, int lbd = -1);
  void end_irred_cls();
  void set_lit_weight(Lit l, const T& w);
  T outer_count();
  uint32_t get_tstamp() const { return tstamp; }
  uint32_t get_tstamp(int32_t lev) const {
    if (dec_level() < lev) return UINT_MAX;
    return decisions[lev].tstamp;
  }

  const CounterConfiguration& get_conf() const { return conf;}
  uint32_t get_num_low_lbds() const { return num_low_lbd_cls; }
  uint32_t get_num_long_red_cls() const { return long_red_cls.size(); }
  uint32_t get_num_irred_long_cls() const { return long_irred_cls.size(); }
  bool get_is_approximate() const { return is_approximate; }

  uint32_t nVars() const { return var_data.size() - 1; }
  double get_start_time() const { return start_time;}
  const auto& get_cache() const { return comp_manager->get_cache();}
  void set_generators(const vector<map<Lit, Lit>>& _gens) { generators = _gens; }

  const T& get_weight(const Lit l) { return weights[l.raw()];}
  T get_weight(const uint32_t v) {
    Lit l(v, false);
    return weights[l.raw()]+weights[l.neg().raw()];}
  bool weight_larger_than(const T& a, const T&b) const {
    if constexpr (!weighted) return a > b;
    else if constexpr (!cpx) return a > b;
    else return a.real()*a.imag() > b.real()*b.imag();
  }
  auto get_indep_support_end() const { return indep_support_end; }
  auto get_opt_indep_support_end() const { return opt_indep_support_end; }
  const auto& get_var_data(uint32_t v) const { return var_data[v]; }
  auto dec_level() const { return decisions.get_decision_level(); }
  void print_trail(bool check_entail = true, bool check_anything = true) const;
  void set_stamp() { decisions[dec_level()].tstamp = ++tstamp; }

private:
  CounterConfiguration conf;
  DataAndStatistics<T> stats;
  bool num_vars_set = false;
  std::mt19937_64 mtrand;
  CMSat::SATSolver* sat_solver = nullptr;
  CompManager<T>* comp_manager = nullptr;
  T count_using_cms();

  // ReduceDB
  bool is_antec_of(ClauseOfs ante_cl, Lit lit) const {
    return var(lit).ante.isAClause() && (var(lit).ante.as_cl() == ante_cl);
  }
  void reduce_db();
  template<class T2> void minimize_uip_cl_with_bins(T2& cl);
  vector<Lit> tmp_minim_with_bins;
  void delete_cl(const ClauseOfs cl_ofs);
  bool red_cl_can_be_deleted(ClauseOfs cl_ofs);
  uint32_t lbd_cutoff;
  uint32_t num_low_lbd_cls = 0; // Last time counted low LBD clauses
  uint32_t num_used_cls = 0; // last time counted used clauses
  uint64_t last_reducedb_confl = 0;
  uint64_t last_reducedb_cls_added = 0;
  uint64_t last_reducedb_dec = 0;

  // Computing LBD (lbd == 2 means "glue clause")
  vector<uint64_t> lbd_helper;
  uint64_t lbd_helper_flag = 0;
  template<class T2> uint32_t calc_lbd(const T2& lits);

  // Clause adding
  void simple_preprocess();
  vector<Lit> unit_cls;
  bool remove_duplicates(vector<Lit>& lits);
  bool exists_unit_cl_of(const Lit l) const {
    for (const auto& l2 : unit_cls) if (l == l2) return true;
    return false;
  }
  template<typename T2> void attach_cl(ClauseOfs off, const T2& lits);
  Clause* add_cl(const vector<Lit> &literals, bool red);
  inline bool add_bin_cl(Lit a, Lit b, bool red);
  bool ended_irred_cls = false; // indicates if we have called end_irred_cls()

  // DNF Cube stuff
  bool restart_if_needed();
  vector<Cube<T>> mini_cubes;
  uint32_t disable_small_cubes(vector<Cube<T>>& cubes);
  void disable_smaller_cube_if_overlap(uint32_t i, uint32_t i2, vector<Cube<T>>& cubes);
  void print_and_check_cubes(vector<Cube<T>>& cubes);
  void disable_cubes_if_overlap(vector<Cube<T>>& cubes);
  void extend_cubes(vector<Cube<T>>& cubes);
  int cube_try_extend_by_lit(const Lit torem, const Cube<T>& c);
  T check_count_norestart(const Cube<T>& c);
  T check_count_norestart_cms(const Cube<T>& c);
  vector<Cube<T>> one_restart_count();
  bool clash_cubes(const set<Lit>& c1, const set<Lit>& c2) const;
  bool compute_cube(Cube<T>& cube, const int side);
  vector<map<Lit, Lit>> generators;
  void symm_cubes(vector<Cube<T>>& cubes);

  const DataAndStatistics<T>& get_stats() const;
  void fill_cl(const Antecedent& ante, Lit*& c, uint32_t& size, Lit p) const;

  //Debug stuff
  T check_count(const bool also_incl_curr_and_later_dec = false);
  bool is_implied(const vector<Lit>& cp);
  void check_implied(const vector<Lit>& cl);
  template<class T2> bool clause_falsified(const T2& cl) const;
  void check_sat_solution() const;
  bool check_watchlists() const;
  template<class T2> void check_cl_propagated_conflicted(T2& cl, uint32_t off = 0) const;
  void check_all_propagated_conflicted() const;
  void print_all_levels();
  void check_cl_unsat(Lit* c, uint32_t size) const;
  void check_trail(bool check_entail = true) const;
  bool find_offs_in_watch(const vec<ClOffsBlckL>& ws, ClauseOfs off) const;
  void check_all_cl_in_watchlists() const;
#ifdef SLOW_DEBUG
  vector<vector<Lit>> debug_irred_cls;
#endif

  // Weights
  uint64_t vars_act_dec_num = 0;
  static constexpr bool cpx = std::is_same<T, complex<mpq_class>>::value;
  static constexpr bool weighted = std::is_same<T, mpq_class>::value || std::is_same<T, complex<mpq_class>>::value;
  vector<T> weights;
  /** Needed to know what variables were active in given decision levels
  / It's needed for weighted counting to know what variable was active in
  / that component (learnt clauses can propagate vars that were not active,
  / and on backtrack, we'd multiply by them, wrongly) **/
  vector<uint64_t> vars_act_dec;

  // Temporaries
  mutable vector<Lit> tmp_lit; //used as temporary binary cl
  vector<uint32_t> to_clear;
  vector<uint8_t> seen;

  // Switch to approxmc
  double start_time;
  volatile bool appmc_timeout_fired = false;
  bool is_approximate = false;
  mpz_class do_appmc_count();

  // SAT solver
  bool ok = true;
  uint32_t qhead = 0;
  vector<VarData> var_data;
  LiteralIndexedVector<TriValue> values;
  int val(Lit lit) const { return values[lit]; }
  int val(uint32_t var) const { return values[Lit(var,1)]; }
  vector<Lit> trail;
  friend class ClauseAllocator<T>;
  ClauseAllocator<T>* alloc;
  vector<ClauseOfs> long_irred_cls;
  vector<ClauseOfs> long_red_cls;
  bool run_sat_solver(RetState& state);
  int32_t sat_start_dec_level = -1;
  inline bool sat_mode() const {
    return sat_start_dec_level != -1 && dec_level() >= sat_start_dec_level;
  }
  vector<int> sat_solution;
  bool propagate(bool out_of_order = false);
  void get_maxlev_maxind(ClauseOfs ofs, int32_t& maxlev, uint32_t& maxind);
  RetState backtrack();
  bool chrono_check();
  int32_t find_lev_to_set(int32_t implied_lit_lev);
  int32_t find_backtrack_level_of_learnt();
  void reduce_db_if_needed();
  void set_lit(const Lit lit, int32_t dec_lev, Antecedent ant = Antecedent());
  void unset_lit(Lit lit);
  LiteralIndexedVector<LitWatchList> watches;
  void go_back_to(int32_t backj);
  void reactivate_comps_and_backtrack_trail([[maybe_unused]] bool check_ws = true);
  inline VarData& var(const Lit lit) { return var_data[lit.var()]; }
  inline VarData& var(const uint32_t v) { return var_data[v]; }
  inline const VarData &var(const uint32_t v) const{ return var_data[v]; }
  inline const VarData &var(const Lit lit) const { return var_data[lit.var()]; }
  inline bool is_true(const Lit &lit) const { return values[lit] == T_TRI; }
  inline bool is_false(Lit lit) { return values[lit] == F_TRI; }
  inline bool is_unknown(Lit lit) const;
  inline bool is_unknown(uint32_t var) const;
  void set_confl_state(Lit a, Lit b);
  void set_confl_state(Clause* cl);

  // Decisions
  void init_decision_stack();
  void init_activity_scores();
  double var_act(const uint32_t v) const {
    return watches[Lit(v, false)].activity + watches[Lit(v, true)].activity; }
  DecisionStack<T> decisions;
  void decide_lit();
  uint32_t find_best_branch(const bool ignore_td = false, const bool also_nonindep = false);
  double score_of(const uint32_t v, bool ignore_td = false) const;
  void vsads_readjust();
  void compute_score(TWD::TreeDecomposition& tdec, uint32_t nodes, bool print = true);
  void td_decompose();
  TWD::TreeDecomposition td_decompose_component(double mult = 1);
  double td_lookahead_score(const uint32_t v, const uint32_t base_comp_tw);
  void recomp_td_weight();
  void inc_act(const Lit lit);
  Heap<VarOrderLt> order_heap; // Only active during SAT solver mode
  bool standard_polarity(const uint32_t var) const;
  bool get_polarity(const uint32_t var) const;
  vector<double> tdscore;
  double td_weight = 1.0;
  int td_width = 10000;
  uint64_t tstamp = 0;
  const Lit &top_dec_lit() const { return *top_declevel_trail_begin(); }
  vector<Lit>::const_iterator top_declevel_trail_begin() const;
  vector<Lit>::iterator top_declevel_trail_begin();
  vector<uint32_t> common_indep_code(const set<uint32_t>& indeps);

  bool is_indep = true; //< We are currently in indep mode
  // the first variable that's NOT in the indep support
  uint32_t indep_support_end = std::numeric_limits<uint32_t>::max();
  // the first variable that's NOT in the opt indep support
  uint32_t opt_indep_support_end = std::numeric_limits<uint32_t>::max();

  // Printing
  string lit_val_str(Lit lit) const;
  string val_to_str(const TriValue& tri) const;
  void print_dec_info() const;
  template<class T2> void print_cl(const T2& cl) const;
  template<class T2> void v_print_cl(const T2& cl) const;
  void print_cl(const Lit* c, uint32_t size) const;
  void print_conflict_info() const;
  void print_comp_stack_info() const;
  void print_stat_line();
  uint64_t next_print_stat_cache = 4LL*1000LL*1000LL;
  uint64_t next_print_stat_confl = 100LL*1000LL;

  // BDD
  bool should_do_buddy_count() const;
  bool do_buddy_count();
#ifdef BUDDY_ENABLED
  uint64_t buddy_count();
  vector<uint32_t> vmap;
  vector<uint32_t> vmap_rev;
  bdd mybdd_two_or(Lit l, Lit r);
#endif

  //  Conflict analysis below
  RetState resolve_conflict();
  struct ConflictData {
    int32_t nHighestLevel = -1;
    bool bOnlyOneLitFromHighest = false;
  };
  ConflictData find_conflict_level(Lit p);
  Lit confl_lit = NOT_A_LIT;
  Antecedent confl;
  vector<Lit> uip_clause;
  void create_uip_cl();
  void minimize_uip_cl();
  vector<Lit> tmp_cl_minim; // Used during minimize_uip_cl
  uint32_t abst_level(const uint32_t x) const;
  bool lit_redundant(Lit p, uint32_t abstract_levels);
  vector<Lit> analyze_stack;
  void recursive_cc_min();
  inline Antecedent add_uip_confl_cl(const vector<Lit> &literals);

  // Vivification
  void vivif_setup();
  bool v_propagate();
  void v_backtrack();
  void v_unset_lit(const Lit l);
  void v_enqueue(const Lit l);
  TriValue v_val(const Lit l) const;
  void v_new_lev();
  void v_backup();
  void v_restore();
  vector<vector<Lit>> v_backup_cls;
  vector<vector<ClOffsBlckL>> v_backup_watches;
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
  void count_loop();
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
  vector<ClauseOfs> occ_cls;
  void backw_susume_cl(ClauseOfs off);
  void backw_susume_cl_with_bin(BinClSub& b);
  void toplevel_full_probe();
  vector<Lit> bothprop_toset;
};

template<typename T>
void Counter<T>::unset_lit(Lit lit) {
    VERBOSE_DEBUG_DO(cout << "Unsetting lit: " << std::setw(8) << lit << endl);
    SLOW_DEBUG_DO(assert(val(lit) == T_TRI));
    var(lit).ante = Antecedent();
    if constexpr (weighted) if(!sat_mode() && get_weight(lit) != get_default_weight<T>()) {
      uint64_t* at = vars_act_dec.data()+dec_level()*(nVars()+1);
      bool in_comp = (at[0] == at[lit.var()]);
      if (in_comp) decisions[dec_level()].include_solution(get_weight(lit));
    }
    var(lit).decision_level = INVALID_DL;
    values[lit] = X_TRI;
    values[lit.neg()] = X_TRI;
  }

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
  if (sat_mode() && order_heap.in_heap(lit.var())) order_heap.increase(lit.var());
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
  stats.learnt_cls_added++;
  Antecedent ante;
  Clause* cl = add_cl(literals, true);
  bool bump_again = false;
  if (cl) {
    auto off = alloc->get_offset(cl);
    long_red_cls.push_back(off);
    cl->lbd = calc_lbd(*cl);
    bump_again = cl->lbd <= (lbd_cutoff+1);
    ante = Antecedent(off);
  } else if (literals.size() == 2){
    ante = Antecedent(literals.back());
    bump_again = true;
    stats.num_bin_red_cls++;
  }
  if (bump_again) for(const auto& l: literals) inc_act(l);
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
  uint32_t orig_size = cl.size();
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
      if (seen[(l2.neg()).raw()]) seen[(l2.neg()).raw()] = 0;
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
  stats.rem_lits_with_bins += orig_size - cl.size();
  stats.rem_lits_tried++;
}

template<typename T>
template<class T2> void Counter<T>::attach_cl(ClauseOfs off, const T2& lits) {
  Lit blck_lit = lits[lits.size()/2];
  watches[lits[0]].add_cl(off, blck_lit);
  watches[lits[1]].add_cl(off, blck_lit);
}

template<typename T>
bool Counter<T>::is_unknown(Lit lit) const {
    SLOW_DEBUG_DO(assert(lit.var() <= nVars()));
    SLOW_DEBUG_DO(assert(lit.var() != 0));
    return values[lit] == X_TRI;
}

template<typename T>
bool Counter<T>::is_unknown(uint32_t var) const {
    SLOW_DEBUG_DO(assert(var != 0));
    return is_unknown(Lit(var, false));
}

template<typename T>
void Counter<T>::set_confl_state(Lit a, Lit b) {
  confl_lit = a;
  confl = Antecedent(b);
}

template<typename T>
void Counter<T>::set_confl_state(Clause* cl) {
  if (cl->red && cl->lbd > this->lbd_cutoff) {
    cl->set_used();
    /* cl->update_lbd(this->calc_lbd(*cl)); */
  }
  confl = Antecedent(this->alloc->get_offset(cl));
  confl_lit = NOT_A_LIT;
}

template<typename T>
vector<Lit>::const_iterator Counter<T>::top_declevel_trail_begin() const {
  return trail.begin() + this->var(decisions.top().var).sublevel;
}

template<typename T>
vector<Lit>::iterator Counter<T>::top_declevel_trail_begin() {
  return trail.begin() + this->var(decisions.top().var).sublevel;
}

template<typename T>
template<class T2>
uint32_t Counter<T>::calc_lbd(const T2& lits) {
  lbd_helper_flag++;
  uint32_t nblevels = 0;
  for(const auto& l: lits) {
    if (val(l) == X_TRI) {nblevels++;continue;}
    int lev = var(l).decision_level;
    if (lev != 0 && lbd_helper[lev] != lbd_helper_flag) {
      lbd_helper[lev] = lbd_helper_flag;
      nblevels++;
      if (nblevels >= 100) { return nblevels; }
    }
  }
  return nblevels;
}

class OuterCounter {
public:
  OuterCounter(const CounterConfiguration& conf, bool weighted, bool cpx = false) {
    if (!weighted) unw_counter = new Counter<mpz_class>(conf);
    else {
      if (!cpx) wq_counter = new Counter<mpq_class>(conf);
      else cpx_counter = new Counter<complex<mpq_class>>(conf);
    }
    release_assert(unw_counter || wq_counter || cpx_counter);
  }
  ~OuterCounter() {
    delete cpx_counter;
    delete unw_counter;
    delete wq_counter;
  }
  void set_generators(const vector<map<Lit, Lit>>& _gens) {
    if (unw_counter) unw_counter->set_generators(_gens);
    if (wq_counter) wq_counter->set_generators(_gens);
    if (cpx_counter) cpx_counter->set_generators(_gens);
  }
  void end_irred_cls() {
    if (unw_counter) unw_counter->end_irred_cls();
    if (wq_counter) wq_counter->end_irred_cls();
    if (cpx_counter) cpx_counter->end_irred_cls();
  }
  void set_indep_support(const set<uint32_t>& indeps) {
    if (unw_counter) unw_counter->set_indep_support(indeps);
    if (wq_counter) wq_counter->set_indep_support(indeps);
    if (cpx_counter) cpx_counter->set_indep_support(indeps);
  }
  complex<mpq_class> cpx_outer_count() { release_assert(cpx_counter); return cpx_counter->outer_count();}
  mpz_class unw_outer_count() { release_assert(unw_counter); return unw_counter->outer_count();}
  mpq_class wq_outer_count() { release_assert(wq_counter); return wq_counter->outer_count();}
  bool add_red_cl(const vector<Lit>& lits, int lbd = -1) {
    if (unw_counter) return unw_counter->add_red_cl(lits, lbd);
    if (wq_counter) return wq_counter->add_red_cl(lits, lbd);
    if (cpx_counter) return cpx_counter->add_red_cl(lits, lbd);
    release_assert(false);
  }
  bool get_is_approximate() const { release_assert(unw_counter); return unw_counter->get_is_approximate();}
  bool add_irred_cl(const vector<Lit>& lits) {
    if (unw_counter) return unw_counter->add_irred_cl(lits);
    if (wq_counter) return wq_counter->add_irred_cl(lits);
    if (cpx_counter) return cpx_counter->add_irred_cl(lits);
    release_assert(false);
  }
  void set_optional_indep_support(const set<uint32_t>& indeps) {
    if (unw_counter) unw_counter->set_optional_indep_support(indeps);
    if (wq_counter) wq_counter->set_optional_indep_support(indeps);
    if (cpx_counter) return cpx_counter->set_optional_indep_support(indeps);
  }
  void set_lit_weight(const Lit l, const complex<mpq_class>& w) {
    release_assert(wq_counter || cpx_counter);
    if (wq_counter) release_assert(w.imag() == 1);
    if (wq_counter) wq_counter->set_lit_weight(l, w.real());
    if (cpx_counter) return cpx_counter->set_lit_weight(l, w);
  }
  void new_vars(const uint32_t n) {
    if (unw_counter) unw_counter->new_vars(n);
    if (wq_counter) wq_counter->new_vars(n);
    if (cpx_counter) cpx_counter->new_vars(n);
  }
  void print_indep_distrib() const {
    if (unw_counter) unw_counter->print_indep_distrib();
    if (wq_counter) wq_counter->print_indep_distrib();
    if (cpx_counter) cpx_counter->print_indep_distrib();
  }
private:
  Counter<complex<mpq_class>>* cpx_counter = nullptr;
  Counter<mpz_class>* unw_counter = nullptr;
  Counter<mpq_class>* wq_counter = nullptr;
};

}
