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
#include <cryptominisat5/cryptominisat.h>

#include "primitive_types.h"
#include "statistics.h"
#include "structures.h"
#include "containers.h"
#include "solver_config.h"

class Instance {
public:
  Instance(const CounterConfiguration& _config) : config_(_config), stats (this, config_) { }
  void new_vars(const uint32_t n);
  void add_irred_cl(const vector<Lit>& lits);
  uint32_t get_num_low_lbds() const { return num_low_lbd_cls; }
  uint32_t get_num_long_reds() const { return red_cls.size(); }
  uint32_t get_num_irred_long_cls() const { return stats.num_long_irred_clauses_; }
protected:
  CounterConfiguration config_;
  void unSet(Lit lit) {
    VERBOSE_DEBUG_DO(cout << "Unsetting lit: " << lit << endl);
    var(lit).ante = Antecedent(NOT_A_CLAUSE);
    var(lit).decision_level = INVALID_DL;
    lit_values_[lit] = X_TRI;
    lit_values_[lit.neg()] = X_TRI;
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
  size_t minimize_cl_with_bins(ClauseOfs off);
  template<class T> void minimize_uip_cl_with_bins(T& cl);
  vector<Lit> tmp_minim_with_bins;
  void markClauseDeleted(ClauseOfs cl_ofs);
  bool red_cl_can_be_deleted(ClauseOfs cl_ofs);


  // Compact the literal pool erasing all the clause
  // information from deleted clauses
  void compactConflictLiteralPool();
  bool findOfsInWatch(const vector<ClOffsBlckL>& ws, ClauseOfs off) const;
  void checkWatchLists() const;

  // Gives ACTUAL number of variables
  uint32_t nVars() const {
    return variables_.size() - 1;
  }

  DataAndStatistics stats;

  /*  lit_pool_: the literals of all clauses are stored here
   *   This includes both IRRED + RED clauses
   *   irred clauses are until irred_lit_pool_size_
   *   INVARIANT: first and last entries of lit_pool_ are a SENTINEL_LIT
   *
   *   Clauses begin with a ClHeader structure followed by the literals
   *   terminated by SENTINEL_LIT
   */
  vector<Lit> lit_pool_;

  // the first variable that is NOT in the independent support
  uint32_t indep_support_end = std::numeric_limits<uint32_t>::max();
  bool perform_projected_counting = false;


  // this is to determine the starting offset of
  // conflict clauses
  uint32_t irred_lit_pool_size_;

  LiteralIndexedVector<vector<ClauseOfs> > occ_lists_; // used ONLY to figure out which
                                                       // literals should we probe
  LiteralIndexedVector<LitWatchList> watches_; // watches
  vector<ClauseOfs> red_cls;
  vector<Lit> unit_clauses_;
  vector<Variable> variables_;
  LiteralIndexedVector<TriValue> lit_values_;
  vector<double> tdscore; // treewidth-decomposition score
  double act_inc = 1.0;
  uint32_t lbd_cutoff = 2;
  uint32_t num_low_lbd_cls = 0; // Last time counted low LBD clauses
  uint32_t num_used_cls = 0; // last time counted used clauses

  // Computing LBD (lbd == 2 means "glue clause")
  vector<uint64_t> lbdHelper;
  uint64_t lbdHelperFlag = 0;

  uint32_t calc_lbd(ClauseOfs offs) {
    lbdHelperFlag++;
    uint32_t nblevels = 0;
    for (auto it = beginOf(offs); *it != SENTINEL_LIT; it++) {
      int lev = var(*it).decision_level;
      if (lev == -1) {nblevels++; continue;} // can happen in weird conflict...
      if (lev != 1 && lbdHelper[lev] != lbdHelperFlag) {
        lbdHelper[lev] = lbdHelperFlag;
        nblevels++;
        if (nblevels >= 1000) { return nblevels; }
      }
    }
    return nblevels;
  }

  template<class T> uint32_t calc_lbd(T& lits) {
    lbdHelperFlag++;
    uint32_t nblevels = 0;
    for(const auto& l: lits) {
      int lev = var(l).decision_level;
      if (lev != 1 && lbdHelper[lev] != lbdHelperFlag) {
        lbdHelper[lev] = lbdHelperFlag;
        nblevels++;
        if (nblevels >= 1000) { return nblevels; }
      }
    }
    return nblevels;
  }

  void increaseActivity(ClauseOfs offs) {
    for (auto it = beginOf(offs); *it != SENTINEL_LIT; it++)
      increaseActivity(*it);
  }

  void inline increaseActivity(const Lit lit) {
    if (config_.do_single_bump && tmp_seen[lit.var()]) return;
    watches_[lit].activity += act_inc;
    if (watches_[lit].activity > 1e100) {
      //rescale
      act_inc *= 1e-90;
      for(auto& v: watches_) v.activity*=1e-90;
    }
  }

  bool isUnitClause(const Lit lit) {
    for (const auto& l : unit_clauses_) if (l == lit) return true;
    return false;
  }

  bool existsUnitClauseOf(VariableIndex v) {
    for (auto l : unit_clauses_) if (l.var() == v) return true;
    return false;
  }

  bool existsUnitClauseOf(Lit l) {
    for (auto l2 : unit_clauses_) if (l == l2) return true;
    return false;
  }

  inline ClauseIndex addClause(const vector<Lit> &literals, bool red);

  // adds a UIP Conflict Clause
  // and returns it as an Antecedent to the first
  // literal stored in literals
  inline Antecedent addUIPConflictClause(const vector<Lit> &literals);

  inline bool add_bin_cl(Lit litA, Lit litB, bool red);

  inline Variable &var(const Lit lit) {
    return variables_[lit.var()];
  }

  inline const Variable &var(const Lit lit) const {
    return variables_[lit.var()];
  }

  LitWatchList & litWatchList(Lit lit) {
    return watches_[lit];
  }

  const LitWatchList & litWatchList(Lit lit) const {
    return watches_[lit];
  }

  inline bool isTrue(const Lit &lit) const {
    return lit_values_[lit] == T_TRI;
  }

  int val(Lit lit) const { return lit_values_[lit]; }
  bool isFalse(Lit lit) {
    return lit_values_[lit] == F_TRI;
  }

  string lit_val_str(Lit lit) const {
    if (lit_values_[lit] == F_TRI)
      return "FALSE";
    else if (lit_values_[lit] == T_TRI)
      return "TRUE";
    else return "UNKN";
  }

  bool isUnknown(Lit lit) const {
    return lit_values_[lit] == X_TRI;
  }

  bool isUnknown(uint32_t var) const {
    return lit_values_[Lit(var, false)] == X_TRI;
  }

  Lit const* beginOf(ClauseOfs cl_ofs) const {
    return lit_pool_.data() + cl_ofs;
  }
  Lit* beginOf(ClauseOfs cl_ofs) {
    return lit_pool_.data() + cl_ofs;
  }

  decltype(lit_pool_.begin()) conflict_clauses_begin() {
     return lit_pool_.begin() + irred_lit_pool_size_;
   }

  const ClHeader &getHeaderOf(ClauseOfs cl_ofs) const {
    return *reinterpret_cast<const ClHeader *>(
        &lit_pool_[cl_ofs - ClHeader::overheadInLits()]);
  }

  ClHeader &getHeaderOf(ClauseOfs cl_ofs) {
    return *reinterpret_cast<ClHeader *>(
        &lit_pool_[cl_ofs - ClHeader::overheadInLits()]);
  }

  bool isSatisfied(ClauseOfs cl_ofs) {
    for (auto lt = beginOf(cl_ofs); *lt != SENTINEL_LIT; lt++)
      if (isTrue(*lt)) return true;
    return false;
  }
protected:
  bool counted_bottom_comp = true; //when false, we MUST take suggested polarities
  vector<uint8_t> target_polar;
  vector<uint8_t> tmp_seen;

  // Cubes
  vector<Lit> largest_cube;
  mpz_class largest_cube_val = 0;
  uint32_t largest_cube_level = 0;

private:
  void parseWithCMS(const std::string& filename);

};

ClauseIndex Instance::addClause(const vector<Lit> &literals, bool red) {
  if (literals.size() == 1) {
    // TODO Deal properly with the situation that opposing unit clauses are learned
    // NOTE that currently this cannot happen, we always check
    //      for at least one solution
    /* assert(!isUnitClause(literals[0].neg())); */
    unit_clauses_.push_back(literals[0]);
    return 0;
  }

  if (literals.size() == 2) {
    add_bin_cl(literals[0], literals[1], red);
    return 0;
  }

  for (uint32_t i = 0; i < ClHeader::overheadInLits(); i++) lit_pool_.push_back(Lit());
  ClauseOfs cl_ofs = lit_pool_.size();

  for (auto l : literals) { lit_pool_.push_back(l); }
  lit_pool_.push_back(SENTINEL_LIT);
  Lit blckLit = literals[literals.size()/2];
  litWatchList(literals[0]).addWatchLinkTo(cl_ofs, blckLit);
  litWatchList(literals[1]).addWatchLinkTo(cl_ofs, blckLit);
  auto& header = getHeaderOf(cl_ofs);
  header = ClHeader(0, red);
  return cl_ofs;
}

Antecedent Instance::addUIPConflictClause(const vector<Lit> &literals) {
  Antecedent ante(NOT_A_CLAUSE);
  stats.num_clauses_learned_++;
  ClauseOfs cl_ofs = addClause(literals, true);
  if (cl_ofs != NOT_A_CLAUSE) {
    red_cls.push_back(cl_ofs);
    auto& header = getHeaderOf(cl_ofs);
    header = ClHeader(calc_lbd(cl_ofs), true);
    ante = Antecedent(cl_ofs);
  } else if (literals.size() == 2){
    ante = Antecedent(literals.back());
    stats.num_binary_red_clauses_++;
  } else if (literals.size() == 1)
    stats.num_unit_red_clauses_++;
  return ante;
}

bool Instance::add_bin_cl(Lit litA, Lit litB, bool red) {
   if (litWatchList(litA).hasBinaryLinkTo(litB)) return false;
   litWatchList(litA).addBinLinkTo(litB, red);
   litWatchList(litB).addBinLinkTo(litA, red);
   return true;
}

template<class T>
void Instance::minimize_uip_cl_with_bins(T& cl) {
  SLOW_DEBUG_DO(for(const auto& s: tmp_seen) assert(s == 0););
  uint32_t rem = 0;
  assert(cl.size() > 0);
  tmp_minim_with_bins.clear();
  for(const auto& l: cl) { tmp_seen[l.toPosInt()] = 1; tmp_minim_with_bins.push_back(l);}
  for(const auto& l: cl) {
  /* { */
    /* Lit l = tmp_minim_with_bins[0]; */
    if (!tmp_seen[l.toPosInt()]) continue;
    const auto& w = watches_[l].binary_links_;
    for(const auto& l2: w) {
      assert(l.var() != l2.var());
      if (tmp_seen[(l2.neg()).toPosInt()]) { tmp_seen[(l2.neg()).toPosInt()] = 0; rem++; }
    }
  }
  cl.clear(); cl.push_back(tmp_minim_with_bins[0]);
  tmp_seen[tmp_minim_with_bins[0].toPosInt()] = 0;
  for(uint32_t i = 1; i < tmp_minim_with_bins.size(); i++) {
    Lit l = tmp_minim_with_bins[i];
    if (tmp_seen[l.toPosInt()]) {
      cl.push_back(l);
      tmp_seen[l.toPosInt()] = 0;
    }
  }
  stats.rem_lits_with_bins+=rem;
  stats.rem_lits_tried++;
}
