/*
 * instance.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#include "instance.h"
#include "structures.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <sys/stat.h>
#include <cryptominisat5/cryptominisat.h>
#include <cryptominisat5/dimacsparser.h>
#include <cryptominisat5/streambuffer.h>

using CMSat::StreamBuffer;
using CMSat::DimacsParser;
using CMSat::SATSolver;

void Instance::cleanClause(ClauseOfs cl_ofs) {
  bool satisfied = false;
  for (auto it = beginOf(cl_ofs); *it != SENTINEL_LIT; it++)
    if (isTrue(*it)) {
      satisfied = true;
      break;
    }
  // mark the clause as empty if satisfied
  if (satisfied) {
    *beginOf(cl_ofs) = SENTINEL_LIT;
    return;
  }
  auto jt = beginOf(cl_ofs);
  auto it = beginOf(cl_ofs);
  // from now, all inactive literals are resolved
  for (; *it != SENTINEL_LIT; it++, jt++) {
    while (*jt != SENTINEL_LIT && !isUnknown(*jt))
      jt++;
    *it = *jt;
    if (*jt == SENTINEL_LIT)
      break;
  }
  unsigned length = it - beginOf(cl_ofs);
  // if it has become a unit clause, it should have already been asserted
  if (length == 1) {
    *beginOf(cl_ofs) = SENTINEL_LIT;
    // if it has become binary, transform it to binary and delete it
  } else if (length == 2) {
    addBinaryClause(*beginOf(cl_ofs), *(beginOf(cl_ofs) + 1));
    *beginOf(cl_ofs) = SENTINEL_LIT;
  }
}

void Instance::compactClauses() {
  vector<ClauseOfs> clause_ofs;
  clause_ofs.reserve(statistics_.num_long_clauses_);

  // clear watch links and occurrence lists
  for (auto it_lit = literal_pool_.begin(); it_lit != literal_pool_.end();
      it_lit++) {
    if (*it_lit == SENTINEL_LIT) {
      if (it_lit + 1 == literal_pool_.end())
        break;
      it_lit += ClauseHeader::overheadInLits();
      clause_ofs.push_back(1 + it_lit - literal_pool_.begin());
    }
  }

  for (auto ofs : clause_ofs)
    cleanClause(ofs);

  for (auto &l : literals_)
    l.resetWatchList();

  occurrence_lists_.clear();
  occurrence_lists_.resize(variables_.size());

  vector<Lit> tmp_pool = literal_pool_;
  literal_pool_.clear();
  literal_pool_.push_back(SENTINEL_LIT);
  ClauseOfs new_ofs;
  unsigned num_clauses = 0;
  for (auto ofs : clause_ofs) {
    auto it = (tmp_pool.begin() + ofs);
    if (*it != SENTINEL_LIT) {
      for (unsigned i = 0; i < ClauseHeader::overheadInLits(); i++)
        literal_pool_.push_back(0);
      new_ofs = literal_pool_.size();
      litWatchList(*it).addWatchLinkTo(new_ofs);
      litWatchList(*(it + 1)).addWatchLinkTo(new_ofs);
      num_clauses++;
      for (; *it != SENTINEL_LIT; it++) {
        literal_pool_.push_back(*it);
        occurrence_lists_[*it].push_back(new_ofs);
      }
      literal_pool_.push_back(SENTINEL_LIT);
    }
  }

  vector<Lit> tmp_bin;
  unsigned bin_links = 0;
  for (auto &l : literals_) {
    tmp_bin.clear();
    for (auto it = l.binary_links_.begin(); *it != SENTINEL_LIT; it++)
      if (isUnknown(*it))
        tmp_bin.push_back(*it);
    bin_links += tmp_bin.size();
    tmp_bin.push_back(SENTINEL_LIT);
    l.binary_links_ = tmp_bin;
  }
  statistics_.num_long_clauses_ = num_clauses;
  statistics_.num_binary_clauses_ = bin_links >> 1;
}

void Instance::compactVariables() {
  var_map.resize(variables_.size(), 0);
  unsigned last_ofs = 0;
  unsigned num_isolated = 0;
  unsigned num_pisolated = 0;
  LiteralIndexedVector<vector<Lit> > _tmp_bin_links(1);
  LiteralIndexedVector<TriValue> _tmp_values = literal_values_;

  for (auto l : literals_)
    _tmp_bin_links.push_back(l.binary_links_);

  assert(_tmp_bin_links.size() == literals_.size());
  for (unsigned v = 1; v < variables_.size(); v++)
    if (isUnknown(v)) {
      if (isolated(v)) {
        if (independent_support_.find(v) != independent_support_.end())
          num_pisolated ++;
        num_isolated++;
        continue;
      }
      last_ofs++;
      var_map[v] = last_ofs;
    }
  vector <unsigned> temp;
  for (auto it=independent_support_.begin(); it!=independent_support_.end(); ++it){
    if(var_map[*it] != 0){
      temp.push_back(var_map[*it]);
    }
  }
  independent_support_.clear();
  for (auto it=temp.begin(); it!=temp.end(); ++it){
    independent_support_.insert(*it);
  }
  variables_.clear();
  variables_.resize(last_ofs + 1);
  occurrence_lists_.clear();
  occurrence_lists_.resize(variables_.size());
  literals_.clear();
  literals_.resize(variables_.size());
  literal_values_.clear();
  literal_values_.resize(variables_.size(), X_TRI);

  Lit newlit;
  for (auto l = Lit(0, false); l != _tmp_bin_links.end_lit(); l.inc()) {
    if (var_map[l.var()] != 0) {
      newlit = Lit(var_map[l.var()], l.sign());
      for (auto it = _tmp_bin_links[l].begin(); *it != SENTINEL_LIT; it++) {
        assert(var_map[it->var()] != 0);
        literals_[newlit].addBinLinkTo(
            Lit(var_map[it->var()], it->sign()));
      }
    }
  }

  vector<ClauseOfs> clause_ofs;
  clause_ofs.reserve(statistics_.num_long_clauses_);
  for (auto it_lit = literal_pool_.begin(); it_lit != literal_pool_.end();
      it_lit++) {
    if (*it_lit == SENTINEL_LIT) {
      if (it_lit + 1 == literal_pool_.end())
        break;
      it_lit += ClauseHeader::overheadInLits();
      clause_ofs.push_back(1 + it_lit - literal_pool_.begin());
    }
  }

  for (auto ofs : clause_ofs) {
    litWatchList(Lit(var_map[beginOf(ofs)->var()], beginOf(ofs)->sign())).addWatchLinkTo(
        ofs);
    litWatchList(Lit(var_map[(beginOf(ofs) + 1)->var()],
            (beginOf(ofs) + 1)->sign())).addWatchLinkTo(ofs);
    for (auto it_lit = beginOf(ofs); *it_lit != SENTINEL_LIT; it_lit++) {
      *it_lit = Lit(var_map[it_lit->var()], it_lit->sign());
      occurrence_lists_[*it_lit].push_back(ofs);
    }
  }

  literal_values_.clear();
  literal_values_.resize(variables_.size(), X_TRI);
  unit_clauses_.clear();

  statistics_.num_variables_ = variables_.size() - 1 + num_isolated;

  statistics_.num_used_variables_ = num_variables();
  statistics_.num_free_variables_ = num_isolated;
  statistics_.num_free_projected_variables_ = num_pisolated;
}

void Instance::compactConflictLiteralPool(){
  auto write_pos = conflict_clauses_begin();
  vector<ClauseOfs> tmp_conflict_clauses = conflict_clauses_;
  conflict_clauses_.clear();
  for(auto clause_ofs: tmp_conflict_clauses){
    auto read_pos = beginOf(clause_ofs) - ClauseHeader::overheadInLits();
    for(unsigned i = 0; i < ClauseHeader::overheadInLits(); i++)
      *(write_pos++) = *(read_pos++);
    ClauseOfs new_ofs =  write_pos - literal_pool_.begin();
    conflict_clauses_.push_back(new_ofs);
    // first substitute antecedent if clause_ofs implied something
    if(isAntecedentOf(clause_ofs, *beginOf(clause_ofs)))
      var(*beginOf(clause_ofs)).ante = Antecedent(new_ofs);

    // now redo the watches
    litWatchList(*beginOf(clause_ofs)).replaceWatchLinkTo(clause_ofs,new_ofs);
    litWatchList(*(beginOf(clause_ofs)+1)).replaceWatchLinkTo(clause_ofs,new_ofs);
    // next, copy clause data
    assert(read_pos == beginOf(clause_ofs));
    while(*read_pos != SENTINEL_LIT)
      *(write_pos++) = *(read_pos++);
    *(write_pos++) = SENTINEL_LIT;
  }
  literal_pool_.erase(write_pos,literal_pool_.end());
}

bool Instance::deleteConflictClauses() {
  statistics_.times_conflict_clauses_cleaned_++;
  vector<ClauseOfs> tmp_conflict_clauses = conflict_clauses_;
  conflict_clauses_.clear();
  vector<double> tmp_ratios;
  double score;
  for(auto clause_ofs: tmp_conflict_clauses){
    score = getHeaderOf(clause_ofs).score();
    tmp_ratios.push_back(score);

  }
  vector<double> tmp_ratiosB = tmp_ratios;

  sort(tmp_ratiosB.begin(), tmp_ratiosB.end());

  double cutoff = tmp_ratiosB[tmp_ratiosB.size()/2];

  for(unsigned i = 0; i < tmp_conflict_clauses.size(); i++){
    if(tmp_ratios[i] < cutoff){
      if(!markClauseDeleted(tmp_conflict_clauses[i]))
        conflict_clauses_.push_back(tmp_conflict_clauses[i]);
    } else
      conflict_clauses_.push_back(tmp_conflict_clauses[i]);
  }
  return true;
}


bool Instance::markClauseDeleted(ClauseOfs cl_ofs){
  // only first literal may possibly have cl_ofs as antecedent
  if(isAntecedentOf(cl_ofs, *beginOf(cl_ofs)))
    return false;

  litWatchList(*beginOf(cl_ofs)).removeWatchLinkTo(cl_ofs);
  litWatchList(*(beginOf(cl_ofs)+1)).removeWatchLinkTo(cl_ofs);
  return true;
}

static Lit cmsLitToG(const CMSat::Lit& l) {
  return Lit(l.var()+1, !l.sign());
}

void Instance::parseWithCMS(const std::string& filename) {
  unsigned verb = 0;
  #ifndef USE_ZLIB
  FILE * in = fopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<FILE*, CMSat::FN>, SATSolver> parser(&solver, NULL, verb);
  #else
  gzFile in = gzopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<gzFile, CMSat::GZ>, SATSolver> parser(&solver, NULL, verb);
  #endif
  if (in == NULL) {
      std::cout << "ERROR! Could not open file '" << filename
      << "' for reading: " << strerror(errno) << endl;
      std::exit(-1);
  }
  if (!parser.parse_DIMACS(in, true)) exit(-1);
  #ifndef USE_ZLIB
  fclose(in);
  #else
  gzclose(in);
  #endif

  if (parser.sampling_vars_found) {
    for(const auto& lit: parser.sampling_vars) independent_support_.insert(lit+1);
  } else {
    for(uint32_t i = 1; i < solver.nVars()+1; i++) independent_support_.insert(i);
  }
}

bool Instance::createfromFile(const std::string &filename) {
  // The solver is empty
  assert(variables_.empty());
  assert(literal_values_.empty());
  assert(occurrence_lists_.empty());
  assert(literals_.empty());
  assert(literal_pool_.empty());
  assert(independent_support_.empty());
  assert(unit_clauses_.empty());
  assert(conflict_clauses_.empty());

  parseWithCMS(filename);

  literal_pool_.push_back(SENTINEL_LIT);
  variables_.push_back(Variable());
  variables_.resize(solver.nVars() + 1);
  literal_values_.resize(solver.nVars() + 1, X_TRI);
  occurrence_lists_.resize(solver.nVars() + 1);
  literals_.resize(solver.nVars() + 1);
  target_polar.resize(solver.nVars() + 1);

  solver.start_getting_small_clauses(
      std::numeric_limits<uint32_t>::max(),
      std::numeric_limits<uint32_t>::max(),
      false);

  statistics_.num_original_clauses_ = 0;
  vector<CMSat::Lit> cms_cl;
  vector<Lit> literals;
  while(solver.get_next_small_clause(cms_cl)) {
    literals.clear();
    for(const auto&l: cms_cl) literals.push_back(cmsLitToG(l));
    statistics_.num_original_clauses_++;
    statistics_.incorporateClauseData(literals);
    ClauseOfs cl_ofs = addClause(literals);
    if (literals.size() >= 3)
      for (const auto& l : literals)
        occurrence_lists_[l].push_back(cl_ofs);
  }
  solver.end_getting_small_clauses();

  statistics_.num_variables_ = solver.nVars();
  statistics_.num_used_variables_ = num_variables();
  statistics_.num_free_variables_ = solver.nVars() - num_variables();
  statistics_.num_unit_clauses_ = unit_clauses_.size();
  original_lit_pool_size_ = literal_pool_.size();

  return solver.okay();
}

