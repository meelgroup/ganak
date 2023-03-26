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

void Instance::compactConflictLiteralPool(){
  auto write_pos = conflict_clauses_begin();
  vector<ClauseOfs> tmp_conflict_clauses = conflict_clauses_;
  conflict_clauses_.clear();
  for(auto clause_ofs: tmp_conflict_clauses){
    auto read_pos = beginOf(clause_ofs) - ClauseHeader::overheadInLits();
    for(uint32_t i = 0; i < ClauseHeader::overheadInLits(); i++)
      *(write_pos++) = *(read_pos++);
    ClauseOfs new_ofs =  write_pos - lit_pool_.begin();
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
  lit_pool_.erase(write_pos,lit_pool_.end());
}

bool Instance::deleteConflictClauses() {
  stats.times_conflict_clauses_cleaned_++;
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

  for(uint32_t i = 0; i < tmp_conflict_clauses.size(); i++){
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
  uint32_t verb = 0;
  #ifndef USE_ZLIB
  FILE * in = fopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<FILE*, CMSat::FN>, SATSolver> parser(&satSolver, NULL, verb);
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

  indep_support_given = parser.sampling_vars_found;
  if (parser.sampling_vars_found) {
    for(const auto& lit: parser.sampling_vars) indep_support_.insert(lit+1);
  } else {
    for(uint32_t i = 1; i < satSolver.nVars()+1; i++) indep_support_.insert(i);
  }
  must_mult_exp2 = parser.must_mult_exp2;
}

bool Instance::createfromFile(const std::string &filename) {
  // The solver is empty
  assert(variables_.empty());
  assert(lit_values_.empty());
  assert(occ_lists_.empty());
  assert(literals_.empty());
  assert(lit_pool_.empty());
  assert(indep_support_.empty());
  assert(unit_clauses_.empty());
  assert(conflict_clauses_.empty());

  parseWithCMS(filename);

  lit_pool_.push_back(SENTINEL_LIT);
  variables_.push_back(Variable());
  variables_.resize(satSolver.nVars() + 1);
  lit_values_.resize(satSolver.nVars() + 1, X_TRI);
  occ_lists_.resize(satSolver.nVars() + 1);
  literals_.resize(satSolver.nVars() + 1);
  target_polar.resize(satSolver.nVars() + 1);
  if (!satSolver.okay()) return satSolver.okay();

  satSolver.start_getting_small_clauses(
      std::numeric_limits<uint32_t>::max(),
      std::numeric_limits<uint32_t>::max(),
      false);

  stats.num_original_clauses_ = 0;
  vector<CMSat::Lit> cms_cl;
  vector<Lit> literals;
  while(satSolver.get_next_small_clause(cms_cl)) {
    literals.clear();
    for(const auto&l: cms_cl) literals.push_back(cmsLitToG(l));
    stats.num_original_clauses_++;
    stats.incorporateClauseData(literals);
    ClauseOfs cl_ofs = addClause(literals);
    if (literals.size() >= 3)
      for (const auto& l : literals)
        occ_lists_[l].push_back(cl_ofs);
  }
  satSolver.end_getting_small_clauses();

  stats.nVars_ = satSolver.nVars();
  stats.num_unit_clauses_ = unit_clauses_.size();
  original_lit_pool_size_ = lit_pool_.size();

  return satSolver.okay();
}
