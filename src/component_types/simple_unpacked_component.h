/*
 * simple_unpacked_component.h
 *
 *  Created on: Feb 16, 2013
 *      Author: mthurley
 */

#ifndef SIMPLE_UNPACKED_COMPONENT_H_
#define SIMPLE_UNPACKED_COMPONENT_H_

#include "cacheable_component.h"
#include "component.h"


#include "../primitive_types.h"

class SimpleUnpackedComponent : public BaseComponent {
public:

    SimpleUnpackedComponent() = default;

    inline SimpleUnpackedComponent(Component &rComp);

    unsigned num_variables() const override {
        return *(data_ + 1);
    }

    unsigned data_size() const {
        return *data_;
    }

    unsigned data_only_byte_size() const override {
        return data_size() * sizeof(unsigned);
    }

    unsigned long raw_data_byte_size() const override {
        return data_size() * sizeof(unsigned)
               + model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t)
               + variables_.size() * sizeof(VariableIndex) + sizeof(cache_hit);
    }

    unsigned long sys_overhead_raw_data_byte_size() const override { //TODO: Fix
        unsigned ds = data_size() * sizeof(unsigned);
        unsigned ms = model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
//      unsigned mask = 0xfffffff8;
//      return (ds & mask) + ((ds & 7)?8:0)
//            +(ms & mask) + ((ms & 7)?8:0);
        unsigned mask = 0xfffffff0;
        return (ds & mask) + ((ds & 15) ? 16 : 0)
               + (ms & mask) + ((ms & 15) ? 16 : 0);
    }

    bool equals(const SimpleUnpackedComponent &comp) const {
        if (hashkey_ != comp.hashkey())
            return false;
        unsigned *p = data_;
        unsigned *r = comp.data_;
        while (p != data_ + data_size()) {
            if (*(p++) != *(r++))
                return false;
        }
        return true;
    }

    void clear() override {
        // before deleting the contents of this component,
        // we should make sure that this component is not present in the component stack anymore!
        assert(isDeletable());
        delete[] data_;
        data_ = nullptr;
    }

private:
    //First element contains the size of data_
    //Second element contains the number of variables
    //Next follow the relevant variables ended by varsSENTINEL
    //Next follow the relevant clauses ended by clsSENTINEL
    unsigned *data_ = nullptr;
};


SimpleUnpackedComponent::SimpleUnpackedComponent(Component &rComp) {
    unsigned data_size = rComp.num_variables() + rComp.numLongClauses() + 2;

    unsigned *p = data_ = new unsigned[data_size];

    *p = data_size;
    *(++p) = rComp.num_variables();
    unsigned hashkey_vars = 0;
    for (auto it = rComp.varsBegin(); *it != varsSENTINEL; it++) {
        *(++p) = *it;
        hashkey_vars = (hashkey_vars * 3) + *it;
    }

    unsigned hashkey_clauses = 0;

    for (auto jt = rComp.clsBegin(); *jt != clsSENTINEL; jt++) {
        *(++p) = *jt;
        hashkey_clauses = (hashkey_clauses * 3) + *jt;
    }
    //*(++p) = 0;

    hashkey_ = hashkey_vars + (((unsigned) hashkey_clauses) << 16);

    assert(p - data_ + 1 == data_size);
}


#endif /* SIMPLE_UNPACKED_COMPONENT_H_ */
