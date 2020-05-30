/*
 * component.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#ifndef COMPONENT_H_
#define COMPONENT_H_

#include <assert.h>
#include <vector>

#include "../primitive_types.h"

using namespace std;

class Component {
public:

  void reserveSpace(unsigned int num_variables, unsigned int num_clauses) {
    data_.reserve(num_variables + num_clauses + 2);
  }

  void set_id(CacheEntryID id) {
    id_ = id;
  }

  CacheEntryID id() const{
    return id_;
  }

  void addVar(const VariableIndex var) {
    // the only time a varsSENTINEL is added should be in a
    // call to closeVariableData(..)
    assert(var != varsSENTINEL);
    data_.push_back(var);
  }

  void closeVariableData() {
    data_.push_back(varsSENTINEL);
    clauses_ofs_ = data_.size();
  }

  void addCl(const ClauseIndex cl) {
    // the only time a clsSENTINEL is added should be in a
    // call to closeClauseData(..)
    assert(cl != clsSENTINEL);
    data_.push_back(cl);
  }

  void closeClauseData() {
    data_.push_back(clsSENTINEL);
    assert(*(clsBegin()-1) == 0);
  }

  vector<VariableIndex>::const_iterator varsBegin() const {
    return data_.begin();
  }

  vector<ClauseIndex>::const_iterator clsBegin() const {
    return data_.begin() + clauses_ofs_;
  }

  unsigned num_variables() const {
    return clauses_ofs_ - 1;
  }

  unsigned numLongClauses() const {
    return data_.size() - clauses_ofs_ - 1;
  }

  bool empty() const {
    return data_.empty();
  }

  void createAsDummyComponent(unsigned max_var_id, unsigned max_clause_id) {
    data_.clear();
    clauses_ofs_ = 1;
    for (unsigned idvar = 1; idvar <= max_var_id; idvar++)
      addVar(idvar);
    closeVariableData();
    if (max_clause_id > 0)
      for (unsigned idcl = 1; idcl <= max_clause_id; idcl++)
        addCl(idcl);
    closeClauseData();
  }

  void clear() {
    clauses_ofs_ = 0;
    data_.clear();
  }

private:
  // data_ stores the component data:
  // for better cache performance the
  // clause and variable data are stored in
  // a contiguous piece of memory
  // variables SENTINEL clauses SENTINEL
  // this order has to be taken care of on filling
  // in the data!
  vector<unsigned> data_;
  unsigned clauses_ofs_ = 0;
  // id_ will identify denote the entry in the cacheable component database,
  // where a Packed version of this component is stored
  // yet this does not imply that the model count of this component is already known
  // once the model count is known, a link to the packed component will be stored
  // in the hash table
  CacheEntryID id_ = 0;
};



#endif /* COMPONENT_H_ */
