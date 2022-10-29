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
protected:

  void unSet(LiteralID lit) {
    var(lit).ante = Antecedent(NOT_A_CLAUSE);
    var(lit).decision_level = INVALID_DL;
    literal_values_[lit] = X_TRI;
    literal_values_[lit.neg()] = X_TRI;
  }

  Antecedent & getAntecedent(LiteralID lit) {
    return variables_[lit.var()].ante;
  }

  bool hasAntecedent(LiteralID lit) {
    return variables_[lit.var()].ante.isAnt();
  }

  bool isAntecedentOf(ClauseOfs ante_cl, LiteralID lit) {
    return var(lit).ante.isAClause() && (var(lit).ante.asCl() == ante_cl);
  }

  bool isolated(VariableIndex v) {
    LiteralID lit(v, false);
    return (litWatchList(lit).binary_links_.size() <= 1)
        && occurrence_lists_[lit].empty()
        && (litWatchList(lit.neg()).binary_links_.size() <= 1)
        && occurrence_lists_[lit.neg()].empty();
  }

  bool free(VariableIndex v) {
    return isolated(v) && isUnknown(v);
  }

  bool deleteConflictClauses();
  bool markClauseDeleted(ClauseOfs cl_ofs);

  // Compact the literal pool erasing all the clause
  // information from deleted clauses
  void compactConflictLiteralPool();

  // we assert that the formula is consistent
  // and has not been found UNSAT yet
  // hard wires all assertions in the literal stack into the formula
  // removes all set variables and essentially reinitiallizes all
  // further data
  void compactClauses();
  void compactVariables();
  void cleanClause(ClauseOfs cl_ofs);

  /////////////////////////////////////////////////////////
  // END access to variables and literals
  /////////////////////////////////////////////////////////


  unsigned int num_conflict_clauses() const {
    return conflict_clauses_.size();
  }

  unsigned int num_variables() {
    return variables_.size() - 1;
  }

  bool createfromFile(const std::string &file_name);
  DataAndStatistics statistics_;

  /** literal_pool_: the literals of all clauses are stored here
   *   INVARIANT: first and last entries of literal_pool_ are a SENTINEL_LIT
   *
   *   Clauses begin with a ClauseHeader structure followed by the literals
   *   terminated by SENTINEL_LIT
   */
  vector<LiteralID> literal_pool_;

  set <unsigned> independent_support_;

  vector<unsigned> var_map;


  // this is to determine the starting offset of
  // conflict clauses
  unsigned original_lit_pool_size_;

  LiteralIndexedVector<LitWatchList> literals_;
  LiteralIndexedVector<vector<ClauseOfs> > occurrence_lists_;
  vector<ClauseOfs> conflict_clauses_;
  vector<LiteralID> unit_clauses_;
  vector<Variable> variables_;
  LiteralIndexedVector<TriValue> literal_values_;

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

  bool isUnitClause(const LiteralID lit) {
    for (auto l : unit_clauses_)
      if (l == lit)
        return true;
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
  bool addUnitClause(const LiteralID lit) {
    for (auto l : unit_clauses_) {
      if (l == lit)
        return true;
      if (l == lit.neg())
        return false;
    }
    unit_clauses_.push_back(lit);
    return true;
  }

  inline ClauseIndex addClause(vector<LiteralID> &literals);

  // adds a UIP Conflict Clause
  // and returns it as an Antecedent to the first
  // literal stored in literals
  inline Antecedent addUIPConflictClause(vector<LiteralID> &literals);

  inline bool addBinaryClause(LiteralID litA, LiteralID litB);

  /////////////////////////////////////////////////////////
  // BEGIN access to variables, literals, clauses
  /////////////////////////////////////////////////////////

  inline Variable &var(const LiteralID lit) {
    return variables_[lit.var()];
  }

  LitWatchList & litWatchList(LiteralID lit) {
    return literals_[lit];
  }

  inline bool isTrue(const LiteralID &lit) const {
    return literal_values_[lit] == T_TRI;
  }

  bool isFalse(LiteralID lit) {
    return literal_values_[lit] == F_TRI;
  }

  bool isUnknown(LiteralID lit) const {
    return literal_values_[lit] == X_TRI;
  }

  vector<LiteralID>::const_iterator beginOf(ClauseOfs cl_ofs) const {
    return literal_pool_.begin() + cl_ofs;
  }
  vector<LiteralID>::iterator beginOf(ClauseOfs cl_ofs) {
    return literal_pool_.begin() + cl_ofs;
  }

  decltype(literal_pool_.begin()) conflict_clauses_begin() {
     return literal_pool_.begin() + original_lit_pool_size_;
   }

  ClauseHeader &getHeaderOf(ClauseOfs cl_ofs) {
    return *reinterpret_cast<ClauseHeader *>(&literal_pool_[cl_ofs
        - ClauseHeader::overheadInLits()]);
  }

  bool isSatisfied(ClauseOfs cl_ofs) {
    for (auto lt = beginOf(cl_ofs); *lt != SENTINEL_LIT; lt++)
      if (isTrue(*lt))
        return true;
    return false;
  }
protected:
  CMSat::SATSolver solver;
  bool counted_bottom_component = false; //when false, we MUST take suggested polarities
  vector<unsigned char> target_polar;
  vector<int> largest_cube;

private:
  void parseWithCMS(const std::string& filename);
};

ClauseIndex Instance::addClause(vector<LiteralID> &literals) {
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
  for (unsigned i = 0; i < ClauseHeader::overheadInLits(); i++)
    literal_pool_.push_back(0);
  ClauseOfs cl_ofs = literal_pool_.size();

  for (auto l : literals) {
    literal_pool_.push_back(l);
    litWatchList(l).increaseActivity(1);
  }
  // make an end: SENTINEL_LIT
  literal_pool_.push_back(SENTINEL_LIT);
  litWatchList(literals[0]).addWatchLinkTo(cl_ofs);
  litWatchList(literals[1]).addWatchLinkTo(cl_ofs);
  getHeaderOf(cl_ofs).set_creation_time(statistics_.num_conflicts_);
  return cl_ofs;
}


Antecedent Instance::addUIPConflictClause(vector<LiteralID> &literals) {
    Antecedent ante(NOT_A_CLAUSE);
    statistics_.num_clauses_learned_++;
    ClauseOfs cl_ofs = addClause(literals);
    if (cl_ofs != 0) {
      conflict_clauses_.push_back(cl_ofs);
      getHeaderOf(cl_ofs).set_length(literals.size());
      ante = Antecedent(cl_ofs);
    } else if (literals.size() == 2){
      ante = Antecedent(literals.back());
      statistics_.num_binary_conflict_clauses_++;
    } else if (literals.size() == 1)
      statistics_.num_unit_clauses_++;
    return ante;
  }

bool Instance::addBinaryClause(LiteralID litA, LiteralID litB) {
   if (litWatchList(litA).hasBinaryLinkTo(litB))
     return false;
   litWatchList(litA).addBinLinkTo(litB);
   litWatchList(litB).addBinLinkTo(litA);
   litWatchList(litA).increaseActivity();
   litWatchList(litB).increaseActivity();
   return true;
 }


#endif /* INSTANCE_H_ */
