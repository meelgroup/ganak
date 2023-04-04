/*
 * difference_packed_comp.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#ifndef DIFFERENCE_PACKED_COMPONENT_H_
#define DIFFERENCE_PACKED_COMPONENT_H_

#include <functional>
#include <set>
#include <iostream>

#include "base_packed_comp.h"
#include "comp.h"
#ifdef DOPCC
#include "../clhash/clhash.h"
#else
#include "../clhash/minim.h"
#endif
#include <math.h>
#include <structures.h>

using std::set;

class DifferencePackedComponent: public BasePackedComponent {
public:

  DifferencePackedComponent() { }
  inline DifferencePackedComponent(Component &rComp, const BPCSizes& sz);
  inline DifferencePackedComponent(vector<void *> & randomseedforCLHASH, Component &rComp, const BPCSizes& sz);
  inline bool contains_any_var(const std::set<uint32_t>& vars, const BPCSizes& sz);

  bool pcc() const {
    return is_pcc;
  }

  uint32_t nVars(const BPCSizes& sz) const{
    if (is_pcc) return old_num_vars;
    uint32_t *p = (uint32_t *) data_;
    return (*p >> sz.bits_of_data_size) & sz.variable_mask;
  }

  uint32_t data_size(const BPCSizes& sz) const {
    if (is_pcc) return old_size;
    return *data_ & sz.data_size_mask;
  }

  uint32_t raw_data_byte_size(const BPCSizes& sz) const {
    if (is_pcc) return num_hash_elems* sizeof(uint64_t) + model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
    else return data_size(sz)* sizeof(uint32_t) + model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
  }

    // raw data size with the overhead
    // for the supposed 16byte alignment of malloc
  uint32_t sys_overhead_raw_data_byte_size(const BPCSizes& sz) const {
    uint32_t ds;
    if (is_pcc) ds = num_hash_elems* sizeof(uint64_t);

    ds = data_size(sz)* sizeof(uint32_t);
    uint32_t ms = model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
    uint32_t mask = 0xfffffff0;
    return (ds & mask) + ((ds & 15)?16:0) +(ms & mask) + ((ms & 15)?16:0);
  }

  bool equals(const DifferencePackedComponent &comp, const BPCSizes& sz) const {
    assert(!is_pcc);
    if(hashkey_ != comp.hashkey()) return false;
    uint32_t* p = data_;
    uint32_t* r = comp.data_;
    while(p != data_ + data_size(sz)) {
        if(*(p++) != *(r++)) return false;
    }
    return true;
  }

#ifdef DOPCC
  uint64_t *compute_clhash(){ return clhashkey_; }
  bool equals(const DifferencePackedComponent &comp, uint64_t* clhash_key) const {
    if(hashkey_ != comp.hashkey()) return false;
    for (uint32_t i=0; i<num_hash_elems;i++){
      if (clhash_key[i] != clhashkey_[i]) return false;
    }
    return true;
  }
#endif
};

bool DifferencePackedComponent::contains_any_var(const std::set<uint32_t>& vars, const BPCSizes& sz) {
  BitStufferReader bs(data_);

  bs.read_bits(sz.bits_of_data_size);
  uint32_t nvars = bs.read_bits(sz.bits_per_variable);
  uint32_t bits_per_var_diff = bs.read_bits(5);
  uint32_t var = bs.read_bits(sz.bits_per_variable);
  if (vars.count(var)) return true;

  // no variables in component??? whatever.
  if (!bits_per_var_diff) return false;

  for (uint32_t i = 1; i < nvars; i++) {
    var += bs.read_bits(bits_per_var_diff) + 1;
    /* cout << "var: " << var << endl; */
    if (vars.count(var)) return true;
  }
  return false;
}

DifferencePackedComponent::DifferencePackedComponent(Component &rComp, const BPCSizes& sz) {
  uint32_t max_var_diff = 0;
  uint32_t hashkey_vars = *rComp.varsBegin();
  for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++) {
    hashkey_vars = (hashkey_vars * 3) + *it;
    if ((*it - *(it - 1)) - 1 > max_var_diff)
      max_var_diff = (*it - *(it - 1)) - 1 ;
  }

  uint32_t hashkey_clauses = *rComp.clsBegin();
  uint32_t max_clause_diff = 0;
  if (*rComp.clsBegin()) {
    for (auto jt = rComp.clsBegin() + 1; *jt != clsSENTINEL; jt++) {
      hashkey_clauses = hashkey_clauses*3 + *jt;
      if (*jt - *(jt - 1) - 1 > max_clause_diff)
        max_clause_diff = *jt - *(jt - 1) - 1;
    }
  }

  hashkey_ = hashkey_vars + ((uint32_t) hashkey_clauses << 11) + ((uint32_t) hashkey_clauses >> 23);

  //VERIFIED the definition of bits_per_var_diff and bits_per_clause_diff
  uint32_t bits_per_var_diff = log2(max_var_diff) + 1;
  uint32_t bits_per_clause_diff = log2(max_clause_diff) + 1;

  assert(bits_per_var_diff <= 31);
  assert(bits_per_clause_diff <= 31);

  uint32_t data_size_vars = sz.bits_of_data_size + 2*sz.bits_per_variable + 5;

  data_size_vars += (rComp.nVars() - 1) * bits_per_var_diff ;

  uint32_t data_size_clauses = 0;
  if(*rComp.clsBegin())
    data_size_clauses += sz.bits_per_clause + 5
       + (rComp.numLongClauses() - 1) * bits_per_clause_diff;

  uint32_t data_size = (data_size_vars + data_size_clauses)/sz.bits_per_block;
    data_size+=  ((data_size_vars + data_size_clauses) % sz.bits_per_block)? 1 : 0;

  data_ = new uint32_t[data_size];

  assert((data_size >> sz.bits_of_data_size) == 0);
  BitStuffer<uint32_t> bs(data_);

  bs.stuff(data_size, sz.bits_of_data_size);
  bs.stuff(rComp.nVars(), sz.bits_per_variable);
  bs.stuff(bits_per_var_diff, 5);
  bs.stuff(*rComp.varsBegin(), sz.bits_per_variable); // 1st var in component
  /* cout << "ORIG data_size: " << data_size << endl; */
  /* cout << "ORIG nvars: " << rComp.nVars() << endl; */
  /* cout << "ORIG bits_per_var_diff: " << bits_per_var_diff << endl; */
  /* cout << "ORIG var: " << *rComp.varsBegin() << endl; */


  // all remaining vars in component
  if(bits_per_var_diff)
    for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++) {
      bs.stuff(*it - *(it - 1) - 1, bits_per_var_diff);
      /* cout << "ORIG var: " << *it << endl; */
    }

  if (*rComp.clsBegin()) {
    bs.stuff(bits_per_clause_diff, 5);
    bs.stuff(*rComp.clsBegin(), sz.bits_per_clause);
    if(bits_per_clause_diff)
      for (auto jt = rComp.clsBegin() + 1; *jt != clsSENTINEL; jt++)
        bs.stuff(*jt - *(jt - 1) - 1, bits_per_clause_diff);
  }

  // to check whether the "END" block of bits_per_clause()
  // many zeros fits into the current
  //bs.end_check(bits_per_clause());
  // this will tell us if we computed the data_size
  // correctly
  bs.assert_size(data_size);
}

DifferencePackedComponent::DifferencePackedComponent(vector<void *>& random, Component &rComp, const BPCSizes& sz) {
  // first, generate hashkey, and compute max diff for cls and vars
  uint32_t max_var_diff = 0;
  uint32_t hashkey_vars = *rComp.varsBegin();
  for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++) {
    hashkey_vars = (hashkey_vars * 3) + *it;
    if ((*it - *(it - 1)) - 1 > max_var_diff)
      max_var_diff = (*it - *(it - 1)) - 1 ;
  }

  uint32_t hashkey_clauses = *rComp.clsBegin();
  uint32_t max_clause_diff = 0;
  if (*rComp.clsBegin()) {
    for (auto jt = rComp.clsBegin() + 1; *jt != clsSENTINEL; jt++) {
      hashkey_clauses = hashkey_clauses*3 + *jt;
      if (*jt - *(jt - 1) - 1 > max_clause_diff)
        max_clause_diff = *jt - *(jt - 1) - 1;
    }
  }

  hashkey_ = hashkey_vars + ((uint32_t) hashkey_clauses << 11) + ((uint32_t) hashkey_clauses >> 23);

  uint32_t bits_per_var_diff = log2(max_var_diff) + 1;
  uint32_t bits_per_clause_diff = log2(max_clause_diff) + 1;

  assert(bits_per_var_diff <= 31);
  assert(bits_per_clause_diff <= 31);

  uint32_t data_size_vars = sz.bits_of_data_size + 2*sz.bits_per_variable + 5;

  data_size_vars += (rComp.nVars() - 1) * bits_per_var_diff ;
  uint32_t data_size_clauses = 0;
  if(*rComp.clsBegin())
    data_size_clauses += sz.bits_per_clause + 5
       + (rComp.numLongClauses() - 1) * bits_per_clause_diff;

  uint32_t data_size = (data_size_vars + data_size_clauses)/sz.bits_per_block;
  data_size += ((data_size_vars + data_size_clauses) % sz.bits_per_block)? 1 : 0;

  data_ = new uint32_t[data_size];
  assert((data_size >> sz.bits_of_data_size) == 0);
  BitStuffer<uint32_t> bs(data_);

  bs.stuff(data_size, sz.bits_of_data_size);
  bs.stuff(rComp.nVars(), sz.bits_per_variable);
  bs.stuff(bits_per_var_diff, 5);
  bs.stuff(*rComp.varsBegin(), sz.bits_per_variable);

  if(bits_per_var_diff)
    for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++)
      bs.stuff(*it - *(it - 1) - 1, bits_per_var_diff);

  if (*rComp.clsBegin()) {
    bs.stuff(bits_per_clause_diff, 5);
    bs.stuff(*rComp.clsBegin(), sz.bits_per_clause);
    if(bits_per_clause_diff)
     for (auto jt = rComp.clsBegin() + 1; *jt != clsSENTINEL; jt++)
      bs.stuff(*jt - *(jt - 1) - 1, bits_per_clause_diff);
  }

  // to check wheter the "END" block of bits_per_clause()
  // many zeros fits into the current
  //bs.end_check(bits_per_clause());
  // this will tell us if we computed the data_size
  // correctly
  bs.assert_size(data_size);

#ifdef DOPCC
  clhashkey_ = new uint64_t[random.size()];
  for(size_t i=0; i<random.size();i++){
    clhasher h(random[i]);
    clhashkey_[i] = h(data_, data_size);
  }
#endif

  num_hash_elems = random.size();
}

#endif /* DIFFERENCE_PACKED_COMPONENT_H_ */
