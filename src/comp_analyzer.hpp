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

#include "common.hpp"
#include "clauseallocator.hpp"
#include "counter_config.hpp"
#include "statistics.hpp"
#include "comp_types/comp.hpp"
#include "comp_types/comp_archetype.hpp"

#include <climits>
#include <cstdint>
#include <map>
#include <gmpxx.h>
#include "containers.hpp"
#include "stack.hpp"
#include "structures.hpp"

using std::map;
using std::pair;

namespace GanakInt {

template<typename T> class ClauseAllocator;
template<typename T> class Counter;

/* #define USE_DIRTY */

struct ClData {
  uint32_t id;
  uint32_t off;
  Lit blk_lit;
  bool operator==(const ClData& other) const {
    return id == other.id && off == other.off && blk_lit == other.blk_lit;
  }
  Lit get_lit1() const { return Lit::toLit(off); }
  Lit get_lit2() const { return blk_lit; }
  bool operator<(const ClData& other) const { return id < other.id; }
};
struct MemData {
  MemData(uint32_t _sz_bin, uint32_t _sz_long) :
    sz_bin(_sz_bin), sz_long(_sz_long)  {}
  MemData() = default;
  uint32_t sz_bin = UINT_MAX;
  uint32_t sz_long = UINT_MAX;
};

struct MyHolder {
  MyHolder () = default;
  ~MyHolder() { delete data;}
  uint32_t* data;
  // start_bin, sz_bin, start_long, sz_long, start_bin, sz_bin.... data...data.... data ... data...
  // start is number of uint32_t-s! not ClData. not bytes.
  //
  ClData& back_long(uint32_t v) {
    return (begin_long(v))[size_long(v)-1];
  }
  ClData* begin_long(uint32_t v) {
    auto start = data[v*4+2];
    return (ClData*) (data + start);
  }
  uint32_t size_long(uint32_t v) { return data[v*4+3];}
  void pop_back_long(uint32_t v) {
    data[v*4+3]--;
  }
  void resize_long(uint32_t v, uint32_t sz) {
    data[v*4+3] = sz;
  }

  //bin
  uint32_t& back_bin(uint32_t v) {
    return (begin_bin(v))[size_bin(v)-1];
  }
  uint32_t* begin_bin(uint32_t v) {
    auto start = data[v*4];
    return (uint32_t*) (data + start);
  }
  uint32_t size_bin(uint32_t v) { return data[v*4+1];}
  void pop_back_bin(uint32_t v) {
    data[v*4+1]--;
  }
  void resize_bin(uint32_t v, uint32_t sz) {
    data[v*4+1] = sz;
  }
};

// There is exactly ONE of this, inside CompManager, which is inside counter
template<typename T>
class CompAnalyzer {
public:
  CompAnalyzer(
        const LiteralIndexedVector<TriValue> & lit_values,
        Counter<T>* _counter);

  auto freq_score_of(uint32_t v) const { return var_freq_scores[v]; }
  void un_bump_score(uint32_t v) { var_freq_scores[v] --; }
  inline void bump_freq_score(uint32_t v) { var_freq_scores[v] ++; }
  const CompArchetype<T>& current_archetype() const { return archetype; }

  void initialize(const LiteralIndexedVector<LitWatchList> & literals,
      const ClauseAllocator<T>* alloc, const vector<ClauseOfs>& long_irred_cls);

  bool var_unvisited_sup_comp(const uint32_t v) const {
    SLOW_DEBUG_DO(assert(v <= max_var));
    return archetype.var_unvisited_in_sup_comp(v);
  }

  // manages the literal whenever it occurs in comp analysis
  // returns true iff the underlying variable was unvisited before
  bool manage_occ_of(const uint32_t v){
    if (archetype.var_unvisited_in_sup_comp(v)) {
      comp_vars.push_back(v);
      archetype.set_var_visited(v);
      __builtin_prefetch(holder.begin_bin(v));
      return true;
    }
    return false;
  }

  bool manage_occ_and_score_of(uint32_t v){
    if (is_unknown(v)) bump_freq_score(v);
    return manage_occ_of(v);
  }

  void setup_analysis_context(StackLevel<T>& top, const Comp & super_comp){
    archetype.re_initialize(top,super_comp);

    debug_print("Setting VAR/CL_SUP_COMP_unvisited for unset vars");
    all_vars_in_comp(super_comp, vt)
      if (is_unknown(*vt)) {
        archetype.set_var_in_sup_comp_unvisited_raw(*vt);
        var_freq_scores[*vt] = 0;
      }

    all_cls_in_comp(super_comp, it) {
      /* if (*it > max_tri_clid) break; */
      archetype.set_clause_in_sup_comp_unvisited(*it);
    }
  }

  bool explore_comp(const uint32_t v, int32_t dec_lev, const uint32_t sup_comp_cls, const uint32_t sup_comp_vars);

  // explore_comp has been called already
  // which set up search_stack, seen[] etc.
  inline Comp *make_comp_from_archetype(){
    auto p =  archetype.make_comp(comp_vars.size());
    return p;
  }

  uint32_t get_max_clid() const { return max_clid; }
  uint32_t get_max_var() const { return max_var; }
  CompArchetype<T>& get_archetype() { return archetype; }

private:
  // the id of the last clause
  // note that clause ID is the clause number,
  // different from the offset of the clause in the literal pool
  uint32_t max_clid = 0;
  uint32_t max_tri_clid = 0;
  uint32_t max_var = 0;
  int32_t backtracked = INT_MAX;
  int32_t last_declev = 0;
  static constexpr bool weighted = std::is_same<T, mpfr::mpreal>::value || std::is_same<T, mpq_class>::value;

  MyHolder holder;
  struct HolderOrigSz {uint32_t sz_bin; uint32_t sz_long;};
  vector<HolderOrigSz> holder_orig_szs;
  vector<uint64_t> last_dirty;
  vector<Lit> long_clauses_data;
  /* vector<vector<MemData>> sz_declevs; */
  vector<int32_t> last_seen;
  const LiteralIndexedVector<TriValue> & values;

  const CounterConfiguration& conf;
  const uint32_t indep_support_end;
  vector<uint32_t> var_freq_scores;
  CompArchetype<T> archetype;
  Counter<T>* counter = nullptr;

  // Quick lookup of cl based on ID

  // Used to figure out which vars are in a component
  // used in  record_comp
  // its size is the number of variables in the component
  vector<uint32_t> comp_vars;


  bool is_false(const Lit lit) const { return values[lit] == F_TRI; }
  bool is_true(const Lit lit) const { return values[lit] == T_TRI; }
  bool is_unknown(const Lit lit) const { return values[lit] == X_TRI; }
  bool is_unknown(const uint32_t v) const { return values[Lit(v, true)] == X_TRI; }
  void bump_var_occs(const uint32_t v);

  // stores all information about the comp of var
  // in variables_seen_, clauses_seen_ and
  // comp_search_stack
  // we have an isolated variable iff
  // after execution comp_search_stack.size()==1
  void record_comp(const uint32_t var, const int32_t declev, const uint32_t sup_comp_cls, const uint32_t sup_comp_vars);

  void get_cl(vector<uint32_t> &tmp, const Clause& cl, const Lit & omit_lit) {
    tmp.clear();
    for (const auto&l: cl) {
      if (l.var() != omit_lit.var()) tmp.push_back(l.raw());
    }
  }

  // This is called from record_comp, i.e. during figuring out what
  // belongs to a component. It's called on every long clause.
  bool search_clause(uint32_t v2, ClData& d, Lit const* cl_start) {
    bool sat = false;
    const auto it_v_end = comp_vars.end();

    for (auto it_l = cl_start; *it_l != SENTINEL_LIT; it_l++) {
      const uint32_t v = it_l->var();
      assert(v <= max_var);
      if (v == v2) continue;

      if (!archetype.var_nil(v)) {
        /* if (archetype.var_visited(v)) cout << "var visited" << endl; */
        /* if (archetype.var_unvisited_in_sup_comp(v)) cout << "var unvisited in sup (must be unknown)" << endl; */
        /* if (archetype.var_in_peer_comp(v)) cout << "var in peer comp" << endl; */
        manage_occ_and_score_of(v);
      }
      else {
        assert(!is_unknown(*it_l));
        if (is_false(*it_l)) continue;
        d.blk_lit = *it_l;

        //accidentally entered a satisfied clause: undo the search process
        /* cout << "satisfied clause due to: " << *it_l << endl; */
        sat = true;
        while (comp_vars.end() != it_v_end) {
          assert(comp_vars.back() <= max_var);
          archetype.set_var_in_sup_comp_unvisited(comp_vars.back());
          comp_vars.pop_back();
        }
        archetype.clear_cl(d.id);
        while(*it_l != SENTINEL_LIT)
          if(is_unknown(*(--it_l))) un_bump_score(it_l->var());
        break;
      }
    }

    if (!archetype.clause_nil(d.id)) {
      bump_freq_score(v2);
      archetype.set_clause_visited(d.id);
    }
    return sat;
  }
};

}
