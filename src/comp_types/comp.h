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

#include <assert.h>
#include <vector>
#include "common.h"
#include "../primitive_types.h"

using std::vector;

class Component {
public:

  void reserveSpace(uint32_t nVars, uint32_t num_clauses) {
    vs_cls_data_.reserve(nVars + num_clauses + 2);
  }

  void set_id(CacheEntryID id) { id_ = id; }
  CacheEntryID id() const { return id_; }

  void addVar(const VariableIndex var) {
    // the only time a varsSENTINEL is added should be in a
    // call to closeVariableData(..)
    SLOW_DEBUG_DO(assert(var != varsSENTINEL));
    vs_cls_data_.push_back(var);
  }

  void closeVariableData() {
    vs_cls_data_.push_back(varsSENTINEL);
    clauses_ofs_ = vs_cls_data_.size();
  }

  void addCl(const ClauseIndex cl) {
    // the only time a clsSENTINEL is added should be in a
    // call to closeClauseData(..)
    SLOW_DEBUG_DO(assert(cl != clsSENTINEL));
    vs_cls_data_.push_back(cl);
  }

  void closeClauseData() {
    vs_cls_data_.push_back(clsSENTINEL);
    assert(*(clsBegin()-1) == 0);
  }

  vector<VariableIndex>::const_iterator varsBegin() const {
    return vs_cls_data_.begin();
  }

  vector<ClauseIndex>::const_iterator clsBegin() const {
    return vs_cls_data_.begin() + clauses_ofs_;
  }

  uint32_t nVars() const { return clauses_ofs_ - 1; }
  uint32_t numLongClauses() const { return vs_cls_data_.size() - clauses_ofs_ - 1; }
  bool empty() const { return vs_cls_data_.empty(); }

  // Creates the full CNF as a component, at start-up. In other words, this is called ONCE
  void createStartingComponent(uint32_t max_var_id, uint32_t max_clause_id) {
    vs_cls_data_.clear();
    clauses_ofs_ = 1;

    // Add all variables to top comp
    for (uint32_t v = 1; v <= max_var_id; v++) addVar(v);
    closeVariableData();

    //Add all clauses to top comp
    for (uint32_t clid = 1; clid <= max_clause_id; clid++) addCl(clid);
    closeClauseData();
  }

  void clear() {
    clauses_ofs_ = 0;
    vs_cls_data_.clear();
  }

private:
  // data_ stores the comp data:
  // for better cache performance the
  // clause and variable data are stored in
  // a contiguous piece of memory
  // variables SENTINEL clauses SENTINEL
  // this order has to be taken care of on filling
  // in the data!
  vector<uint32_t> vs_cls_data_;
  uint32_t clauses_ofs_ = 0;
  // id_ will identify denote the entry in the cacheable comp database,
  // where a Packed version of this comp is stored
  // yet this does not imply that the model count of this comp is already known
  // once the model count is known, a link to the packed comp will be stored
  // in the hash table
  CacheEntryID id_ = 0;
};
