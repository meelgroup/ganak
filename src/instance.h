/*
 * instance.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#ifndef INSTANCE_H_
#define INSTANCE_H_

#include "statistics.h"
#include "structures.h"
#include "containers.h"
#include <set>
#include <cassert>
#include <cryptominisat5/cryptominisat.h>

using std::set;

class Instance {
public:
  Instance() : stats (this) { }
  bool get_indep_support_given() const { return indep_support_given; }
  const set<uint32_t>& get_indep_support() const { return indep_support_; }
  uint32_t get_must_mult_exp2() const { return must_mult_exp2; }
protected:

  void unSet(Lit lit) {
    var(lit).ante = Antecedent(NOT_A_CLAUSE);
    var(lit).decision_level = INVALID_DL;
    lit_values_[lit] = X_TRI;
    lit_values_[lit.neg()] = X_TRI;
  }

  Antecedent & getAntecedent(Lit lit) {
    return variables_[lit.var()].ante;
  }

  bool hasAntecedent(Lit lit) const {
    return variables_[lit.var()].ante.isAnt();
  }

  bool isAntecedentOf(ClauseOfs ante_cl, Lit lit) {
    return var(lit).ante.isAClause() && (var(lit).ante.asCl() == ante_cl);
  }

  bool isolated(VariableIndex v) {
    Lit lit(v, false);
    return (litWatchList(lit).binary_links_.size() <= 1)
        && occ_lists_[lit].empty()
        && (litWatchList(lit.neg()).binary_links_.size() <= 1)
        && occ_lists_[lit.neg()].empty();
  }

  bool free(VariableIndex v) {
    return isolated(v) && isUnknown(v);
  }

  bool deleteConflictClauses();
  bool markClauseDeleted(ClauseOfs cl_ofs);

  // Compact the literal pool erasing all the clause
  // information from deleted clauses
  void compactConflictLiteralPool();

  uint32_t num_conflict_clauses() const {
    return conflict_clauses_.size();
  }

  uint32_t nVars() {
    return variables_.size() - 1;
  }

  bool createfromFile(const std::string &file_name);
  DataAndStatistics stats;

  /** lit_pool_: the literals of all clauses are stored here
   *   INVARIANT: first and last entries of lit_pool_ are a SENTINEL_LIT
   *
   *   Clauses begin with a ClauseHeader structure followed by the literals
   *   terminated by SENTINEL_LIT
   */
  vector<Lit> lit_pool_;

  set<uint32_t> indep_support_;
  bool indep_support_given = false;
  uint32_t must_mult_exp2 = 0;


  // this is to determine the starting offset of
  // conflict clauses
  uint32_t original_lit_pool_size_;

  LiteralIndexedVector<LitWatchList> literals_;
  LiteralIndexedVector<vector<ClauseOfs> > occ_lists_;
  vector<ClauseOfs> conflict_clauses_;
  vector<Lit> unit_clauses_;
  vector<Variable> variables_;
  LiteralIndexedVector<TriValue> lit_values_;

  void decayActivities() {
    for (auto l_it = literals_.begin(); l_it != literals_.end(); l_it++)
      l_it->activity_score_ *= 0.5;

    for(auto clause_ofs: conflict_clauses_)
        getHeaderOf(clause_ofs).decayScore();

  }

  void updateActivities(ClauseOfs clause_ofs) {
    getHeaderOf(clause_ofs).increaseScore();
    for (auto it = beginOf(clause_ofs); *it != SENTINEL_LIT; it++) {
      litWatchList(*it).increaseActivity();
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

  // addUnitClause checks whether lit or lit.neg() is already a
  // unit clause
  // a negative return value implied that the Instance is UNSAT
  bool addUnitClause(const Lit lit) {
    for (auto l : unit_clauses_) {
      if (l == lit)
        return true;
      if (l == lit.neg())
        return false;
    }
    unit_clauses_.push_back(lit);
    return true;
  }

  inline ClauseIndex addClause(const vector<Lit> &literals);

  // adds a UIP Conflict Clause
  // and returns it as an Antecedent to the first
  // literal stored in literals
  inline Antecedent addUIPConflictClause(vector<Lit> &literals);

  inline bool addBinaryClause(Lit litA, Lit litB);

  /////////////////////////////////////////////////////////
  // BEGIN access to variables, literals, clauses
  /////////////////////////////////////////////////////////

  inline Variable &var(const Lit lit) {
    return variables_[lit.var()];
  }

  LitWatchList & litWatchList(Lit lit) {
    return literals_[lit];
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

  vector<Lit>::const_iterator beginOf(ClauseOfs cl_ofs) const {
    return lit_pool_.begin() + cl_ofs;
  }
  vector<Lit>::iterator beginOf(ClauseOfs cl_ofs) {
    return lit_pool_.begin() + cl_ofs;
  }

  decltype(lit_pool_.begin()) conflict_clauses_begin() {
     return lit_pool_.begin() + original_lit_pool_size_;
   }

  ClauseHeader &getHeaderOf(ClauseOfs cl_ofs) {
    return *reinterpret_cast<ClauseHeader *>(&lit_pool_[cl_ofs
        - ClauseHeader::overheadInLits()]);
  }

  bool isSatisfied(ClauseOfs cl_ofs) {
    for (auto lt = beginOf(cl_ofs); *lt != SENTINEL_LIT; lt++)
      if (isTrue(*lt))
        return true;
    return false;
  }
protected:
  CMSat::SATSolver satSolver;
  bool counted_bottom_comp = true; //when false, we MUST take suggested polarities
  vector<uint8_t> target_polar;
  vector<Lit> largest_cube;
  mpz_class largest_cube_val = 0;

private:
  void parseWithCMS(const std::string& filename);

};

ClauseIndex Instance::addClause(const vector<Lit> &literals) {
  if (literals.size() == 1) {
    //TODO Deal properly with the situation that opposing unit clauses are learned
    // assert(!isUnitClause(literals[0].neg()));
    unit_clauses_.push_back(literals[0]);
    return 0;
  }
  if (literals.size() == 2) {
    addBinaryClause(literals[0], literals[1]);
    return 0;
  }
  for (uint32_t i = 0; i < ClauseHeader::overheadInLits(); i++)
    lit_pool_.push_back(0);
  ClauseOfs cl_ofs = lit_pool_.size();

  for (auto l : literals) {
    lit_pool_.push_back(l);
    litWatchList(l).increaseActivity(1);
  }
  // make an end: SENTINEL_LIT
  lit_pool_.push_back(SENTINEL_LIT);
  litWatchList(literals[0]).addWatchLinkTo(cl_ofs);
  litWatchList(literals[1]).addWatchLinkTo(cl_ofs);
  getHeaderOf(cl_ofs).set_creation_time(stats.num_conflicts_);
  return cl_ofs;
}


Antecedent Instance::addUIPConflictClause(vector<Lit> &literals) {
    Antecedent ante(NOT_A_CLAUSE);
    stats.num_clauses_learned_++;
    ClauseOfs cl_ofs = addClause(literals);
    if (cl_ofs != 0) {
      conflict_clauses_.push_back(cl_ofs);
      getHeaderOf(cl_ofs).set_length(literals.size());
      ante = Antecedent(cl_ofs);
    } else if (literals.size() == 2){
      ante = Antecedent(literals.back());
      stats.num_binary_conflict_clauses_++;
    } else if (literals.size() == 1)
      stats.num_unit_clauses_++;
    return ante;
}

bool Instance::addBinaryClause(Lit litA, Lit litB) {
   if (litWatchList(litA).hasBinaryLinkTo(litB))
     return false;
   litWatchList(litA).addBinLinkTo(litB);
   litWatchList(litB).addBinLinkTo(litA);
   litWatchList(litA).increaseActivity();
   litWatchList(litB).increaseActivity();
   return true;
}


#endif /* INSTANCE_H_ */
