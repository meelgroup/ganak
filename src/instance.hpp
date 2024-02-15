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

#include <cassert>
#include <utility>
#include <cryptominisat5/cryptominisat.h>

#include "clauseallocator.hpp"
#include "primitive_types.hpp"
#include "statistics.hpp"
#include "structures.hpp"
#include "containers.hpp"
#include "counter_config.hpp"

class ClauseAllocator;

class Instance {
public:
  Instance(const CounterConfiguration& _conf);
  ~Instance();
  void new_vars(const uint32_t n);
  uint32_t get_num_low_lbds() const { return num_low_lbd_cls; }
  uint32_t get_num_long_reds() const { return longRedCls.size(); }
  uint32_t get_num_irred_long_cls() const { return long_irred_cls.size(); }
  int val(Lit lit) const { return values[lit]; }
  int val(uint32_t var) const { return values[Lit(var,1)]; }


  friend class ClauseAllocator;
  ClauseAllocator* alloc;
  vector<ClauseOfs> long_irred_cls;
  vector<ClauseOfs> longRedCls;
#ifdef SLOW_DEBUG
  vector<vector<Lit>> debug_irred_cls;
#endif
protected:
  CounterConfiguration conf;
  void unSet(Lit lit) {
    VERBOSE_DEBUG_DO(cout << "Unsetting lit: " << lit << endl);
    var(lit).ante = Antecedent();
    var(lit).decision_level = INVALID_DL;
    values[lit] = X_TRI;
    values[lit.neg()] = X_TRI;
  }

  const Antecedent & getAntecedent(Lit lit) const {
    return variables_[lit.var()].ante;
  }

  bool hasAntecedent(Lit lit) const {
    return variables_[lit.var()].ante.isAnt();
  }

  bool isAntecedentOf(ClauseOfs ante_cl, Lit lit) const {
    return var(lit).ante.isAClause() && (var(lit).ante.asCl() == ante_cl);
  }

  void reduceDB();
  template<class T> void minimize_uip_cl_with_bins(T& cl);
  vector<Lit> tmp_minim_with_bins;
  void markClauseDeleted(const ClauseOfs cl_ofs);
  bool red_cl_can_be_deleted(ClauseOfs cl_ofs);


  bool findOfsInWatch(const vector<ClOffsBlckL>& ws, ClauseOfs off) const;
  void checkWatchLists() const;

  // Gives ACTUAL number of variables
  uint32_t nVars() const {
    return variables_.size() - 1;
  }

  DataAndStatistics stats;

  // the first variable that is NOT in the independent support
  uint32_t indep_support_end = std::numeric_limits<uint32_t>::max();

  LiteralIndexedVector<LitWatchList> watches;
  double max_activity = 0;
  vector<Lit> unit_clauses_;
  vector<VarData> variables_;
  bool num_vars_set = false;
  LiteralIndexedVector<TriValue> values;
  vector<double> tdscore;
  vector<double> tdscore2;
  double act_inc = 1.0;
  uint32_t lbd_cutoff = 2;
  uint32_t num_low_lbd_cls = 0; // Last time counted low LBD clauses
  uint32_t num_used_cls = 0; // last time counted used clauses

  // Computing LBD (lbd == 2 means "glue clause")
  vector<uint64_t> lbdHelper;
  uint64_t lbdHelperFlag = 0;

  template<class T>
  uint32_t calc_lbd(const T& lits) {
    lbdHelperFlag++;
    uint32_t nblevels = 0;
    for(const auto& l: lits) {
      if (val(l) == X_TRI) {nblevels++;continue;}
      int lev = var(l).decision_level;
      if (lev != 0 && lbdHelper[lev] != lbdHelperFlag) {
        lbdHelper[lev] = lbdHelperFlag;
        nblevels++;
        if (nblevels >= 250) { return nblevels; }
      }
    }
    return nblevels;
  }

  bool existsUnitClauseOf(const Lit l) const {
    for (const auto& l2 : unit_clauses_) if (l == l2) return true;
    return false;
  }

  template<class T> void attach_cl(ClauseOfs off, const T& lits);
  Clause* addClause(const vector<Lit> &literals, bool red);

  // adds a UIP Conflict Clause
  // and returns it as an Antecedent to the first
  // literal stored in literals
  inline Antecedent addUIPConflictClause(const vector<Lit> &literals);

  inline bool add_bin_cl(Lit a, Lit b, bool red);

  inline VarData &var(const Lit lit) {
    return variables_[lit.var()];
  }

  inline VarData &var(const uint32_t v) {
    return variables_[v];
  }

  inline const VarData &var(const Lit lit) const {
    return variables_[lit.var()];
  }

  LitWatchList & litWatchList(Lit lit) {
    return watches[lit];
  }

  const LitWatchList & litWatchList(Lit lit) const {
    return watches[lit];
  }

  inline bool is_true(const Lit &lit) const {
    return values[lit] == T_TRI;
  }

  bool is_false(Lit lit) {
    return values[lit] == F_TRI;
  }

  string lit_val_str(Lit lit) const {
    if (values[lit] == F_TRI)
      return "FALSE";
    else if (values[lit] == T_TRI)
      return "TRUE";
    else return "UNKN";
  }

  string val_to_str(const TriValue& tri) const {
    if (tri == F_TRI)
      return "FALSE";
    else if (tri == T_TRI)
      return "TRUE";
    else return "UNKN";
  }

  bool is_unknown(Lit lit) const {
    return values[lit] == X_TRI;
  }

  bool is_unknown(uint32_t var) const {
    return is_unknown(Lit(var, false));
  }

  bool is_satisfied(ClauseOfs off) {
    for (auto lt: *alloc->ptr(off)) if (is_true(lt)) return true;
    return false;
  }
protected:
  bool counted_bottom_comp = true; //when false, we MUST take suggested polarities
  vector<uint8_t> seen;
  vector<Cube> mini_cubes;

private:
  void parse_with_cms(const std::string& filename);

};

Antecedent Instance::addUIPConflictClause(const vector<Lit> &literals) {
  Antecedent ante;
  stats.num_clauses_learned_++;
  Clause* cl = addClause(literals, true);
  if (cl) {
    auto off = alloc->get_offset(cl);
    longRedCls.push_back(off);
    cl->lbd = calc_lbd(*cl);
    ante = Antecedent(off);
  } else if (literals.size() == 2){
    ante = Antecedent(literals.back());
    stats.num_binary_red_clauses_++;
  }
  return ante;
}

bool Instance::add_bin_cl(Lit a, Lit b, bool red) {
   litWatchList(a).addBinLinkTo(b, red);
   litWatchList(b).addBinLinkTo(a, red);
   return true;
}

template<class T>
void Instance::minimize_uip_cl_with_bins(T& cl) {
  SLOW_DEBUG_DO(for(const auto& s: seen) assert(s == 0););
  uint32_t rem = 0;
  assert(cl.size() > 0);
  tmp_minim_with_bins.clear();
  for(const auto& l: cl) { seen[l.toPosInt()] = 1; tmp_minim_with_bins.push_back(l);}
  for(const auto& l: cl) {
  /* { */
    /* Lit l = tmp_minim_with_bins[0]; */
    if (!seen[l.toPosInt()]) continue;
    const auto& w = watches[l].binary_links_;
    for(const auto& bincl: w) {
      const auto& l2 = bincl.lit();
      assert(l.var() != l2.var());
      if (seen[(l2.neg()).toPosInt()]) { seen[(l2.neg()).toPosInt()] = 0; rem++; }
    }
  }
  cl.clear(); cl.push_back(tmp_minim_with_bins[0]);
  seen[tmp_minim_with_bins[0].toPosInt()] = 0;
  for(uint32_t i = 1; i < tmp_minim_with_bins.size(); i++) {
    Lit l = tmp_minim_with_bins[i];
    if (seen[l.toPosInt()]) {
      cl.push_back(l);
      seen[l.toPosInt()] = 0;
    }
  }
  stats.rem_lits_with_bins+=rem;
  stats.rem_lits_tried++;
}

template<class T> void Instance::attach_cl(ClauseOfs off, const T& lits) {
  Lit blck_lit = lits[lits.size()/2];
  litWatchList(lits[0]).addWatchLinkTo(off, blck_lit);
  litWatchList(lits[1]).addWatchLinkTo(off, blck_lit);
}
