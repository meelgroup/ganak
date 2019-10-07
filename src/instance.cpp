/*
 * instance.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#include "instance.h"

#include <algorithm>
#include <fstream>
#include <sys/stat.h>


using namespace std;

void Instance::cleanClause(ClauseOfs cl_ofs) {
  bool satisfied = false;
  for (auto it = beginOf(cl_ofs); *it != SENTINEL_LIT; it++)
    if (isSatisfied(*it)) {
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
    while (*jt != SENTINEL_LIT && !isActive(*jt))
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
    if (!isoCC){
    *beginOf(cl_ofs) = SENTINEL_LIT;
    }
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

  // cout << "[Instance::compactClauses()]"<< endl;
  vector<LiteralID> tmp_pool = literal_pool_;
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
      literal(*it).addWatchLinkTo(new_ofs);
      literal(*(it + 1)).addWatchLinkTo(new_ofs);
      num_clauses++;
      for (; *it != SENTINEL_LIT; it++) {
        // cout << (*it).val()<< " ";
        literal_pool_.push_back(*it);
        occurrence_lists_[*it].push_back(new_ofs);
      }
      // cout << "0 "<< endl;
      literal_pool_.push_back(SENTINEL_LIT);
    }
  }

  vector<LiteralID> tmp_bin;
  unsigned bin_links = 0;
  for (auto &l : literals_) {
    tmp_bin.clear();
    for (auto it = l.binary_links_.begin(); *it != SENTINEL_LIT; it++)
      if (isActive(*it)){
        // cout << "BIN "<< (*it).val() << " 0"<< endl;  
        tmp_bin.push_back(*it);
      }
    bin_links += tmp_bin.size();
    tmp_bin.push_back(SENTINEL_LIT);
    l.binary_links_ = tmp_bin;
  }
  statistics_.num_long_clauses_ = num_clauses;
  statistics_.num_binary_clauses_ = bin_links >> 1;
}

void Instance::compactVariables(vector<unsigned> *var_map, vector<unsigned> *rev_map) {
  
  vector<Variable> temp_variables_;
  for (auto v: variables_){
    temp_variables_.push_back(v);
  }
  
  var_map->resize(variables_.size(), 0);
  rev_map->resize(variables_.size(), 0);
  unsigned last_ofs = 0;
  unsigned num_isolated = 0;
  unsigned num_isolated_special = 0;
  LiteralIndexedVector<vector<LiteralID> > _tmp_bin_links(1);
  LiteralIndexedVector<TriValue> _tmp_values = literal_values_;

  for (auto l : literals_)
    _tmp_bin_links.push_back(l.binary_links_);

  assert(_tmp_bin_links.size() == literals_.size());
  for (unsigned v = 1; v < variables_.size(); v++)
    if (isActive(v)) {
      if (isolated(v)) {
        if(!mpf_cmp_d(variables_[v].get_weight().get_mpf_t(), 2)){
          num_isolated++;
        }
        else{
          num_isolated_special++;
        }
        continue;
      }
      last_ofs++;
      (*var_map)[v] = last_ofs;
      (*rev_map)[last_ofs] = v;
    }
  vector <unsigned> temp;
  for (auto it=statistics_.independent_support_.begin(); it!=statistics_.independent_support_.end(); ++it){
    if((*var_map)[*it] != 0){
      temp.push_back((*var_map)[*it]);
    }
  }
  statistics_.independent_support_.clear();
  for (auto it=temp.begin(); it!=temp.end(); ++it){
    statistics_.independent_support_.insert(*it);
  }
  variables_.clear();
  variables_.push_back(Variable());
  for (int i = 1; i <= last_ofs;i++){
    variables_.push_back(temp_variables_[(*rev_map)[i]]);
    // cout << "var map "<< i << " " << (*rev_map)[i]<<endl; 
  }
  // cout << endl;
  // variables_.resize(last_ofs + 1);
  // int i =0;
  // for (auto v : variables_){
    // cout<< "var weight: " << i <<" "<< variables_[i].get_weight()<< endl;
    // i++;
  // }
  occurrence_lists_.clear();
  occurrence_lists_.resize(variables_.size());
  literals_.clear();
  literals_.resize(variables_.size());
  literal_values_.clear();
  literal_values_.resize(variables_.size(), X_TRI);

  unsigned bin_links = 0;
  LiteralID newlit;
  for (auto l = LiteralID(0, false); l != _tmp_bin_links.end_lit(); l.inc()) {
    if ((*var_map)[l.var()] != 0) {
      newlit = LiteralID((*var_map)[l.var()], l.sign());
      for (auto it = _tmp_bin_links[l].begin(); *it != SENTINEL_LIT; it++) {
        assert((*var_map)[it->var()] != 0);
        // cout << "[Compact Instance::addBinaryClause] ";
        // literal_pool_.push_back(newlit);
        // literal_pool_.push_back(LiteralID((*var_map)[it->var()], it->sign()));
        // literal_pool_.push_back(SENTINEL_LIT);
        // cout << newlit.val() << " " << LiteralID((*var_map)[it->var()], it->sign()).val() << endl;
        literals_[newlit].addBinLinkTo(
            LiteralID((*var_map)[it->var()], it->sign()));
      }
      bin_links += literals_[newlit].binary_links_.size() - 1;
    }
  }

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

  for (auto ofs : clause_ofs) {
    literal(LiteralID((*var_map)[beginOf(ofs)->var()], beginOf(ofs)->sign())).addWatchLinkTo(
        ofs);
    literal(LiteralID((*var_map)[(beginOf(ofs) + 1)->var()],
            (beginOf(ofs) + 1)->sign())).addWatchLinkTo(ofs);
    for (auto it_lit = beginOf(ofs); *it_lit != SENTINEL_LIT; it_lit++) {
      *it_lit = LiteralID((*var_map)[it_lit->var()], it_lit->sign());
      occurrence_lists_[*it_lit].push_back(ofs);
    }
  }

  literal_values_.clear();
  literal_values_.resize(variables_.size(), X_TRI);
  unit_clauses_.clear();
  statistics_.num_variables_ = variables_.size() - 1 + num_isolated + num_isolated_special;

  statistics_.num_used_variables_ = num_variables();
  statistics_.num_free_unweighted_variables_ = num_isolated;
  statistics_.num_free_weighted_variables_ = num_isolated_special;
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
    literal(*beginOf(clause_ofs)).replaceWatchLinkTo(clause_ofs,new_ofs);
    literal(*(beginOf(clause_ofs)+1)).replaceWatchLinkTo(clause_ofs,new_ofs);
    // next, copy clause data
    assert(read_pos == beginOf(clause_ofs));
    while(*read_pos != SENTINEL_LIT)
      *(write_pos++) = *(read_pos++);
    *(write_pos++) = SENTINEL_LIT;
  }
  literal_pool_.erase(write_pos,literal_pool_.end());
}



//bool Instance::deleteConflictClauses() {
//  statistics_.times_conflict_clauses_cleaned_++;
//  vector<ClauseOfs> tmp_conflict_clauses = conflict_clauses_;
//  conflict_clauses_.clear();
//  vector<double> tmp_ratios;
//  double score, lifetime;
//  for(auto clause_ofs: tmp_conflict_clauses){
//    score = getHeaderOf(clause_ofs).score();
//    lifetime = statistics_.num_conflicts_ - getHeaderOf(clause_ofs).creation_time();
//    tmp_ratios.push_back(score/lifetime/(getHeaderOf(clause_ofs).length()));
//  }
//  vector<double> tmp_ratiosB = tmp_ratios;
//
//  sort(tmp_ratiosB.begin(), tmp_ratiosB.end());
//
//  double cutoff = tmp_ratiosB[tmp_ratiosB.size()/2];
//
//  for(unsigned i = 0; i < tmp_conflict_clauses.size(); i++){
//    if(tmp_ratios[i] < cutoff){
//      if(!markClauseDeleted(tmp_conflict_clauses[i]))
//        conflict_clauses_.push_back(tmp_conflict_clauses[i]);
//    } else
//      conflict_clauses_.push_back(tmp_conflict_clauses[i]);
//  }
//  return true;
//}

bool Instance::deleteConflictClauses() {
  statistics_.times_conflict_clauses_cleaned_++;
  vector<ClauseOfs> tmp_conflict_clauses = conflict_clauses_;
  conflict_clauses_.clear();
  vector<double> tmp_ratios;
  double score;
  // double lifetime;
  for(auto clause_ofs: tmp_conflict_clauses){
    score = getHeaderOf(clause_ofs).score();
    // lifetime = statistics_.num_conflicts_ - getHeaderOf(clause_ofs).creation_time();
   // tmp_ratios.push_back(score/lifetime);
    tmp_ratios.push_back(score);

  }
  vector<double> tmp_ratiosB = tmp_ratios;

  sort(tmp_ratiosB.begin(), tmp_ratiosB.end());

  double cutoff = tmp_ratiosB[tmp_ratiosB.size()/2];

  for(unsigned i = 0; i < tmp_conflict_clauses.size(); i++){
    // cout << " The LBD score is " << getHeaderOf(tmp_conflict_clauses[i]).LBD_score() << endl;
    if(tmp_ratios[i] < cutoff && 
    getHeaderOf(tmp_conflict_clauses[i]).LBD_score() > 3){
      if(!markClauseDeleted(tmp_conflict_clauses[i])){
        conflict_clauses_.push_back(tmp_conflict_clauses[i]);
      }
    } 
    else{
      conflict_clauses_.push_back(tmp_conflict_clauses[i]);
    }
  }
  return true;
}


bool Instance::markClauseDeleted(ClauseOfs cl_ofs){
  // only first literal may possibly have cl_ofs as antecedent
  if(isAntecedentOf(cl_ofs, *beginOf(cl_ofs)))
    return false;

  literal(*beginOf(cl_ofs)).removeWatchLinkTo(cl_ofs);
  literal(*(beginOf(cl_ofs)+1)).removeWatchLinkTo(cl_ofs);
  return true;
}

bool Instance::create_independent_support_from_file(const string &file_name){
  char v;
  int var;
  struct stat filestatus;
  statistics_.independent_support_.clear();
  ///Begin File input
  ifstream input_file(file_name);
  if (!input_file) {
    cerr << "Cannot open file: " << file_name <<" using all the variables as Independent Support"<< endl;
    for (unsigned v = 1; v < variables_.size(); v++){
      statistics_.independent_support_.insert(v);
    }
    return 1;
  }
  stat(file_name.c_str(), &filestatus);
  if (input_file >> v && v != 'v'){
    cerr << "Invalid Independent Support file using all the variables as Independent Support" << endl;
    for (unsigned v = 1; v < variables_.size(); v++){
      statistics_.independent_support_.insert(v);
    }
    ///END NEW
    input_file.close();
  //  /// END FILE input
    return 1;
  }
  bool novarpresent = true;
  while ((input_file >> var) && var != 0){
    novarpresent = false;
    statistics_.independent_support_.insert(var);
    // independent_support_[var] = true; 
  }
  if (novarpresent){
    cerr << "Invalid Independent Support file using all the variables as Independent Support" << endl;
    for (unsigned v = 1; v < variables_.size(); v++){
      statistics_.independent_support_.insert(v);
    }
    ///END NEW
    // cout << "Reading of Independent Support finished sucessfully " << endl;
    input_file.close();
  //  /// END FILE input
    return 1;
  }
  cout << "Reading of Independent Support finished sucessfully " << endl;
  input_file.close();
  return 1;
}

bool Instance::createfromFile(const string &file_name) {
  unsigned int nVars, nCls;
  unsigned nVarsWeighted = 0;
  int lit;
  mpf_class weight;
  unsigned max_ignore = 1000000;
  unsigned clauses_added = 0;
  LiteralID llit;
  vector<LiteralID> literals;
  string idstring;
  char c;

  // clear everything
  literal_pool_.clear();
  literal_pool_.push_back(SENTINEL_LIT);

  variables_.clear();
  variables_.push_back(Variable()); //initializing the Sentinel
  literal_values_.clear();
  unit_clauses_.clear();

  ///BEGIN File input
  ifstream input_file(file_name);
  if (!input_file) {
    cerr << "Cannot open file: " << file_name << endl;
    exit(0);
  }

  struct stat filestatus;
  stat(file_name.c_str(), &filestatus);
  variables_.reserve(10000);
  variables_.resize(10000);
  literals.reserve(10000);

  while (input_file >> c && c != 'p'){
    if (c == 'w'){
      input_file >> lit;
      input_file >> weight;
      if (abs(lit) >= variables_.size()){
        int newsize = 2*variables_.size()>lit ? 2*variables_.size() : lit+1;
        variables_.resize(newsize);   // table doubling or (literal_index+1), whichever is greater
      }
      if (lit < 0){
        variables_[-1*lit].assign_weight(1-weight);
      }
      else{
        variables_[lit].assign_weight(weight);
      }
      // cout << c << " " << lit << " "<< weight<<" "<<1-weight<< endl;
    }
    input_file.ignore(max_ignore, '\n');
  }
  if (!(input_file >> idstring && idstring == "cnf" && input_file >> nVars
      && input_file >> nCls)) {
    cerr << "Invalid CNF file" << endl;
    exit(0);
  }

  variables_.resize(nVars + 1);
  literal_values_.resize(nVars + 1, X_TRI);
  literal_pool_.reserve(filestatus.st_size);
  conflict_clauses_.reserve(2*nCls);
  occurrence_lists_.clear();
  occurrence_lists_.resize(nVars + 1);

  literals_.clear();
  literals_.resize(nVars + 1);

  while ((input_file >> c) && clauses_added < nCls) {
    if (c == 'w'){
      input_file >> lit;
      input_file >> weight;
      if (lit < 0){
        variables_[-1*lit].assign_weight(1-weight);
      }
      else{
        variables_[lit].assign_weight(weight);
      }
      // cout << c << " " << lit << " "<< weight<<" "<<1-weight<< endl;
    }
    else{
      input_file.unget(); //extracted a nonspace character to determine if we have a clause, so put it back
      if ((c == '-') || isdigit(c)) {
        literals.clear();
        bool skip_clause = false;
        while ((input_file >> lit) && lit != 0) {
          bool duplicate_literal = false;
          for (auto i : literals) {
            if (i.toInt() == lit) {
              duplicate_literal = true;
              break;
            }
            if (i.toInt() == -lit) {
              skip_clause = true;
              break;
            }
          }
          if (!duplicate_literal) {
            literals.push_back(lit);
          }
        }
        if (!skip_clause) {
          assert(!literals.empty());
          clauses_added++;
          statistics_.incorporateClauseData(literals);
          ClauseOfs cl_ofs = addClause(literals);
          if (literals.size() >= 3)
            for (auto l : literals)
              occurrence_lists_[l].push_back(cl_ofs);
        }
      }
    }
    input_file.ignore(max_ignore, '\n');
  }
  input_file.unget();
  while (input_file >> c) {
    if (c == 'w'){
      input_file >> lit;
      input_file >> weight;
      if (lit < 0){
        variables_[-1*lit].assign_weight(1-weight);
      }
      else{
        variables_[lit].assign_weight(weight);
      }
      // cout << c << " " << lit << " "<< weight<<" "<<1-weight<< endl;
    }
    input_file.ignore(max_ignore, '\n');
  }
  // cout << "weights"<<endl;
  // int i = 0;
  // for (auto v: variables_){
  //   if(i > 0){
  //     cout <<i << " " << v.weight_ << endl;
  //   }
  //   i++;
  // }
  ///END NEW
  input_file.close();
  //  /// END FILE input

  statistics_.num_variables_ = statistics_.num_original_variables_ = nVars;
  statistics_.num_used_variables_ = num_variables();
  statistics_.num_free_unweighted_variables_ = nVars - num_variables();

  statistics_.num_original_clauses_ = nCls;

  statistics_.num_original_binary_clauses_ = statistics_.num_binary_clauses_;
  statistics_.num_original_unit_clauses_ = statistics_.num_unit_clauses_ =
      unit_clauses_.size();

  original_lit_pool_size_ = literal_pool_.size();
  return true;
}

