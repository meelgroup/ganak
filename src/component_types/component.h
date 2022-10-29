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

using std::vector;

class Component {
public:

  void reserveSpace(unsigned int num_variables, unsigned int num_clauses) {
    vs_cls_data_.reserve(num_variables + num_clauses + 2);
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
    vs_cls_data_.push_back(var);
  }

  void closeVariableData() {
    vs_cls_data_.push_back(varsSENTINEL);
    clauses_ofs_ = vs_cls_data_.size();
  }

  void addCl(const ClauseIndex cl) {
    // the only time a clsSENTINEL is added should be in a
    // call to closeClauseData(..)
    assert(cl != clsSENTINEL);
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

  unsigned num_variables() const {
    return clauses_ofs_ - 1;
  }

  unsigned numLongClauses() const {
    return vs_cls_data_.size() - clauses_ofs_ - 1;
  }

  bool empty() const {
    return vs_cls_data_.empty();
  }

  // Creates the full CNF, at startup
  void createAsDummyComponent(unsigned max_var_id, unsigned max_clause_id) {
    vs_cls_data_.clear();
    clauses_ofs_ = 1;

    // Add all variables to top component
    for (unsigned idvar = 1; idvar <= max_var_id; idvar++)
      addVar(idvar);
    closeVariableData();

    //Add all clauses to top component
    if (max_clause_id > 0)
      for (unsigned idcl = 1; idcl <= max_clause_id; idcl++)
        addCl(idcl);
    closeClauseData();
  }

  void clear() {
    clauses_ofs_ = 0;
    vs_cls_data_.clear();
  }

private:
  // data_ stores the component data:
  // for better cache performance the
  // clause and variable data are stored in
  // a contiguous piece of memory
  // variables SENTINEL clauses SENTINEL
  // this order has to be taken care of on filling
  // in the data!
  vector<unsigned> vs_cls_data_;
  unsigned clauses_ofs_ = 0;
  // id_ will identify denote the entry in the cacheable component database,
  // where a Packed version of this component is stored
  // yet this does not imply that the model count of this component is already known
  // once the model count is known, a link to the packed component will be stored
  // in the hash table
  CacheEntryID id_ = 0;
};



#endif /* COMPONENT_H_ */
