/*
 * instance.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#pragma once

#include <set>
#include <cassert>
#include <cryptominisat5/cryptominisat.h>

#include "statistics.h"
#include "structures.h"
#include "containers.h"
#include "solver_config.h"

using std::set;

class Instance {
public:
  Instance(const CounterConfiguration& _config) : config_(_config), stats (this) { }
  void new_vars(const uint32_t n);
  void add_irred_cl(const vector<Lit>& lits);
  size_t num_conflict_clauses() const { return conflict_clauses_.size(); }
  uint32_t num_conflict_clauses_compacted() const { return num_conflict_clauses_compacted_; }

protected:
  CounterConfiguration config_;
  void unSet(Lit lit) {
    var(lit).ante = Antecedent(NOT_A_CLAUSE);
    var(lit).bprop = false;
    var(lit).decision_level = INVALID_DL;
    lit_values_[lit] = X_TRI;
    lit_values_[lit.neg()] = X_TRI;
  }

  const Antecedent & getAntecedent(Lit lit) const {
    return variables_[lit.var()].ante;
  }

  bool antedecentBProp(Lit lit) const {
    return variables_[lit.var()].bprop;
  }

  bool hasAntecedent(Lit lit) const {
    return variables_[lit.var()].bprop || variables_[lit.var()].ante.isAnt();
  }

  bool isAntecedentOf(ClauseOfs ante_cl, Lit lit) {
    return var(lit).ante.isAClause() && (var(lit).ante.asCl() == ante_cl);
  }

  bool deleteConflictClauses();
  bool markClauseDeleted(ClauseOfs cl_ofs);

  // Compact the literal pool erasing all the clause
  // information from deleted clauses
  void compactConflictLiteralPool();

  uint32_t nVars() const {
    return variables_.size() - 1;
  }

  DataAndStatistics stats;

  /*  lit_pool_: the literals of all clauses are stored here
   *   This includes both IRRED + RED clauses
   *   irred clauses are until irred_lit_pool_size_
   *   INVARIANT: first and last entries of lit_pool_ are a SENTINEL_LIT
   *
   *   Clauses begin with a ClauseHeader structure followed by the literals
   *   terminated by SENTINEL_LIT
   */
  vector<Lit> lit_pool_;

  // the first variable that is NOT in the independent support
  uint32_t indep_support_end = std::numeric_limits<uint32_t>::max();


  // this is to determine the starting offset of
  // conflict clauses
  uint32_t irred_lit_pool_size_;

  LiteralIndexedVector<vector<ClauseOfs> > occ_lists_; // used ONLY to figure out which
                                                       // literals should we probe
  LiteralIndexedVector<LitWatchList> watches_; // watches
  vector<ClauseOfs> conflict_clauses_;
  uint32_t num_conflict_clauses_compacted_ = 0;
  vector<Lit> unit_clauses_;
  vector<Variable> variables_;
  LiteralIndexedVector<TriValue> lit_values_;
  vector<double> exscore;
  double act_inc = 1.0;

  void decayActivities(bool also_watches) {
    if (also_watches) for (auto& w: watches_) w.activity *= 0.5;
    for(auto clause_ofs: conflict_clauses_)
        getHeaderOf(clause_ofs).decayScore();
  }

  void updateActivities(ClauseOfs clause_ofs ) {
    getHeaderOf(clause_ofs).increaseScore();
    for (auto it = beginOf(clause_ofs); *it != SENTINEL_LIT; it++) {
      increaseActivity(*it);
    }
  }

  void inline increaseActivity(const Lit lit)
  {
    if (config_.do_single_bump && tmp_seen[lit.var()]) return;
    watches_[lit].activity += act_inc;
    if (watches_[lit].activity > 1e100) {
      //rescale
      act_inc *= 1e-90;
      for(auto& v: watches_) v.activity*=1e-90;
    }
  }

  bool isUnitClause(const Lit lit) {
    for (const auto& l : unit_clauses_)
      if (l == lit) return true;
    return false;
  }

  bool existsUnitClauseOf(VariableIndex v) {
    for (auto l : unit_clauses_)
      if (l.var() == v)
        return true;
    return false;
  }

  inline ClauseIndex addClause(const vector<Lit> &literals, bool irred);

  // adds a UIP Conflict Clause
  // and returns it as an Antecedent to the first
  // literal stored in literals
  inline Antecedent addUIPConflictClause(const vector<Lit> &literals);

  inline bool add_bin_cl(Lit litA, Lit litB, bool irred);

  inline Variable &var(const Lit lit) {
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

  bool isFalse(Lit lit) {
    return lit_values_[lit] == F_TRI;
  }

  bool isUnknown(Lit lit) const {
    return lit_values_[lit] == X_TRI;
  }

  bool isUnknown(uint32_t var) const {
    return lit_values_[Lit(var, false)] == X_TRI;
  }

  vector<Lit>::const_iterator beginOf(ClauseOfs cl_ofs) const {
    return lit_pool_.begin() + cl_ofs;
  }
  vector<Lit>::iterator beginOf(ClauseOfs cl_ofs) {
    return lit_pool_.begin() + cl_ofs;
  }

  decltype(lit_pool_.begin()) conflict_clauses_begin() {
     return lit_pool_.begin() + irred_lit_pool_size_;
   }

  ClauseHeader &getHeaderOf(ClauseOfs cl_ofs) {
    return *reinterpret_cast<ClauseHeader *>(
        &lit_pool_[cl_ofs - ClauseHeader::overheadInLits()]);
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

ClauseIndex Instance::addClause(const vector<Lit> &literals, bool irred) {
  if (literals.size() == 1) {
    //TODO Deal properly with the situation that opposing unit clauses are learned
    // assert(!isUnitClause(literals[0].neg()));
    unit_clauses_.push_back(literals[0]);
    return 0;
  }

  if (literals.size() == 2) {
    add_bin_cl(literals[0], literals[1], irred);
    return 0;
  }

  for (uint32_t i = 0; i < ClauseHeader::overheadInLits(); i++) lit_pool_.push_back(Lit());
  ClauseOfs cl_ofs = lit_pool_.size();

  for (auto l : literals) {
    lit_pool_.push_back(l);
  }
  lit_pool_.push_back(SENTINEL_LIT);
  Lit blckLit = literals[literals.size()/2];
  litWatchList(literals[0]).addWatchLinkTo(cl_ofs, blckLit);
  litWatchList(literals[1]).addWatchLinkTo(cl_ofs, blckLit);
  getHeaderOf(cl_ofs).set_creation_time(stats.num_conflicts_);
  return cl_ofs;
}

Antecedent Instance::addUIPConflictClause(const vector<Lit> &literals) {
    Antecedent ante(NOT_A_CLAUSE);
    stats.num_clauses_learned_++;
    ClauseOfs cl_ofs = addClause(literals, false);
    if (cl_ofs != 0) {
      conflict_clauses_.push_back(cl_ofs);
      getHeaderOf(cl_ofs).set_length(literals.size());
      ante = Antecedent(cl_ofs);
    } else if (literals.size() == 2){
      /* cout << "Binary learnt: " << literals[0] << " " << literals[1] << endl; */
      ante = Antecedent(literals.back());
      stats.num_binary_red_clauses_++;
    } else if (literals.size() == 1)
      /* cout << "Unit learnt: " << literals[0] << endl; */
      stats.num_unit_red_clauses_++;
    return ante;
}

bool Instance::add_bin_cl(Lit litA, Lit litB, bool irred) {
   if (litWatchList(litA).hasBinaryLinkTo(litB)) return false;
   litWatchList(litA).addBinLinkTo(litB, irred);
   litWatchList(litB).addBinLinkTo(litA, irred);
   return true;
}
