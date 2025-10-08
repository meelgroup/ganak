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

class ClauseAllocator;
class Counter;

constexpr uint32_t hstride = 10; // stamp(2), lev, 2*(start, size, orig_size) which is 9, but we round up to 10 for alignment
/* #define ANALYZE_DEBUG */

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
  ~MyHolder() { delete [] data;}
  uint32_t* data;
  // start_bin, sz_bin, start_long, sz_long, start_bin, sz_bin.... data...data.... data ... data...
  // start is number of uint32_t-s! not ClData. not bytes.

  uint64_t tstamp(uint32_t v) const {return *(uint64_t*)(data+(v*hstride));}
  uint64_t& tstamp(uint32_t v) {return *(uint64_t*)(data+(v*hstride));}
  int32_t lev(uint32_t v) const {return (int32_t)data[v*hstride+2];}
  void set_lev(uint32_t v, int32_t lev) {data[v*hstride+2] = lev;}
  static constexpr uint32_t offset = 3;

  //
  // bin
  uint32_t* begin_bin(uint32_t v) {
    auto start = data[v*hstride+offset+0];
    return (uint32_t*) (data+start);
  }
  uint32_t& back_bin(uint32_t v) {
    return (begin_bin(v))[size_bin(v)-1];
  }
  uint32_t size_bin(uint32_t v) const { return data[v*hstride+offset+1];}
  uint32_t& size_bin(uint32_t v) { return data[v*hstride+offset+1];}
  uint32_t orig_size_bin(uint32_t v) const { return data[v*hstride+offset+2];}
  uint32_t& orig_size_bin(uint32_t v) { return data[v*hstride+offset+2];}
  void pop_back_bin(uint32_t v) { size_bin(v)--; }
  void resize_bin(uint32_t v, uint32_t sz) { size_bin(v) = sz; }

  // long
  ClData* begin_long(uint32_t v) {
    auto start = data[v*hstride+offset+3];
    return (ClData*) (data + start);
  }
  ClData& back_long(uint32_t v) {
    return (begin_long(v))[size_long(v)-1];
  }
  uint32_t size_long(uint32_t v) const { return data[v*hstride+offset+4];}
  uint32_t& size_long(uint32_t v) { return data[v*hstride+offset+4];}
  uint32_t orig_size_long(uint32_t v) const { return data[v*hstride+offset+5];}
  uint32_t& orig_size_long(uint32_t v) { return data[v*hstride+offset+5];}
  void pop_back_long(uint32_t v) { size_long(v)--; }
  void resize_long(uint32_t v, uint32_t sz) { size_long(v) = sz; }
};

// There is exactly ONE of this, inside CompManager, which is inside counter
class CompAnalyzer {
public:
  CompAnalyzer(
        const LiteralIndexedVector<TriValue> & lit_values,
        Counter* _counter);

  auto freq_score_of(uint32_t v) const { return var_freq_scores[v]; }
  inline void bump_freq_score(uint32_t v) { var_freq_scores[v] ++; }
  const CompArchetype& current_archetype() const { return archetype; }

  void initialize(const LiteralIndexedVector<LitWatchList> & literals,
      const ClauseAllocator* alloc, const vector<ClauseOfs>& long_irred_cls);

  bool var_unvisited_in_sup_comp(const uint32_t v) const {
    SLOW_DEBUG_DO(assert(v <= max_var));
    return archetype.var_unvisited_in_sup_comp(v);
  }

  // manages the literal whenever it occurs in comp analysis
  // returns true iff the underlying variable was unvisited before
  void manage_occ_of(const uint32_t v){
    // if unvisited, it MUST be unknown
    if (archetype.var_unvisited_in_sup_comp(v)) {
      __builtin_prefetch(holder.begin_bin(v));
      comp_vars.push_back(v);
      archetype.set_var_visited(v);
    }
  }

  // Sometimes, it's cheaper to look up the lit than the variable,
  // because we have already looked up the literal, so it's in the cache
  void manage_occ_and_score_of(Lit l) {
    if (is_unknown(l)) {
      bump_freq_score(l.var());
      manage_occ_of(l.var());
    }
  }

  void setup_analysis_context(StackLevel& top, const Comp& super_comp){
    archetype.re_initialize(top,super_comp);
    debug_print("Setting VAR/CL_SUP_COMP_unvisited for unset vars");
    all_vars_in_comp(super_comp, vt) if (is_unknown(*vt)) {
      archetype.set_var_in_sup_comp_unvisited(*vt);
      var_freq_scores[*vt] = 0;
    }
    all_cls_in_comp(super_comp, it) archetype.set_clause_in_sup_comp_unvisited(*it);
  }

  bool explore_comp(const uint32_t v, const uint32_t sup_comp_long_cls, const uint32_t sup_comp_bin_cls);

  // explore_comp has been called already
  // which set up search_stack, seen[] etc.
  inline Comp *make_comp_from_archetype(){
    SLOW_DEBUG_DO(for (auto&v: comp_vars) assert(is_unknown(v)));
    auto p = archetype.make_comp(comp_vars.size());
    return p;
  }

  uint32_t get_max_clid() const { return max_clid; }
  uint32_t get_bin_cls() const { return archetype.num_bin_cls; }
  uint32_t get_max_var() const { return max_var; }
  CompArchetype& get_archetype() { return archetype; }

private:
  // the id of the last clause
  // note that clause ID is the clause number,
  // different from the offset of the clause in the literal pool
  uint32_t max_clid = 0;
  uint32_t max_tri_clid = 0;
  uint32_t max_var = 0;
  DataAndStatistics& stats;

  MyHolder holder;
  vector<Lit> long_clauses_data;
  const LiteralIndexedVector<TriValue> & values;

  const CounterConfiguration& conf;
  const uint32_t indep_support_end;
  vector<uint32_t> var_freq_scores;
  CompArchetype archetype;
  Counter* counter = nullptr;

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
  void record_comp(const uint32_t var, const uint32_t sup_comp_cls, const uint32_t sup_comp_bin_cls);

  void get_cl(vector<uint32_t> &tmp, const Clause& cl, const Lit & omit_lit) {
    tmp.clear();
    for (const auto&l: cl) {
      if (l.var() != omit_lit.var()) tmp.push_back(l.raw());
    }
  }

  // This is called from record_comp, i.e. during figuring out what
  // belongs to a component. It's called on every long clause.
  // The clause is _definitely_ in the supercomponent
  bool search_clause(ClData& d, Lit const* cl_start) {
    bool sat = false;
    for (auto it_l = cl_start; *it_l != SENTINEL_LIT; it_l++) {
        if (is_true(*it_l)) {sat = true; break;}
    }

    if (sat) {
      archetype.set_cl_clear(d.id);
      return true;
    }

    for (auto it_l = cl_start; *it_l != SENTINEL_LIT; it_l++) {
      SLOW_DEBUG_DO(
          const uint32_t v = it_l->var();
          assert(v <= max_var);
          assert(is_false(*it_l) || archetype.var_unvisited_in_sup_comp(v) || archetype.var_visited(v)));
      manage_occ_and_score_of(*it_l);
    }

    archetype.set_clause_visited(d.id);
    return false;
  }
};

}
