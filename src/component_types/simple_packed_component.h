/*
 * simple_packed_component.h
 *
 *  Created on: Feb 4, 2013
 *      Author: mthurley
 */

#ifndef SIMPLE_PACKED_COMPONENT_H_
#define SIMPLE_PACKED_COMPONENT_H_

#include "base_packed_component.h"
#include "component.h"


#include "../primitive_types.h"

class SimplePackedComponent : public BasePackedComponent {
public:

  SimplePackedComponent() {
  }

  inline SimplePackedComponent(Component &rComp);

  unsigned num_variables() const{
    uint64_t *p = (uint64_t *) data_;
    return (*p >> bits_of_data_size()) & (uint64_t) variable_mask();
  }

  unsigned data_size() const {
       return *data_ & _data_size_mask;
  }

  unsigned data_only_byte_size() const {
      return data_size()* sizeof(unsigned);
  }

  unsigned raw_data_byte_size() const {
        return data_size()* sizeof(unsigned)
             + model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
  }

  bool equals(const SimplePackedComponent &comp) const {
    if(hashkey_ != comp.hashkey())
      return false;
    unsigned* p = data_;
    unsigned* r = comp.data_;
    while(p != data_ + data_size()) {
        if(*(p++) != *(r++))
            return false;
    }
    return true;
  }

};

SimplePackedComponent::SimplePackedComponent(Component &rComp) {
  unsigned data_size_vars = bits_of_data_size() + bits_per_variable() + rComp.num_variables() * bits_per_variable();
  unsigned data_size_clauses = rComp.numLongClauses() * bits_per_clause();
  unsigned data_size = (data_size_vars + data_size_clauses)/bits_per_block();

  data_size+=  ((data_size_vars + data_size_clauses) % bits_per_block())? 1 : 0;

  data_ =  new unsigned[data_size];

  BitStuffer<unsigned> bs(data_);
  unsigned hashkey_vars = 0;
  unsigned hashkey_clauses = 0;

  assert((data_size >> bits_of_data_size()) == 0);

  bs.stuff(data_size, bits_of_data_size());
  bs.stuff(rComp.num_variables(),bits_per_variable());

  for (auto it = rComp.varsBegin(); *it != varsSENTINEL; it++) {
    hashkey_vars = (hashkey_vars *3) + *it;
    bs.stuff(*it, bits_per_variable());
  }

  if (*rComp.clsBegin())
    for (auto jt = rComp.clsBegin(); *jt != clsSENTINEL; jt++) {
      hashkey_clauses = (hashkey_clauses *3) + *jt;
      bs.stuff(*jt, bits_per_clause());
    }
  bs.assert_size(data_size);

  hashkey_ = hashkey_vars + (((unsigned) hashkey_clauses) << 16);
}

#endif /* SIMPLE_PACKED_COMPONENT_H_ */
