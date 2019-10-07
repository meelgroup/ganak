/*
 * difference_packed_component.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#ifndef DIFFERENCE_PACKED_COMPONENT_H_
#define DIFFERENCE_PACKED_COMPONENT_H_

#include "base_packed_component.h"
#include "component.h"

#include "../clhash/clhash.h"
#include <math.h>



class DifferencePackedComponent:public BasePackedComponent {
public:

  DifferencePackedComponent() {
  }
  inline DifferencePackedComponent(unsigned int graphhash, unsigned int num_variables);
  inline DifferencePackedComponent(void * random, Component &rComp);
  inline DifferencePackedComponent(Component &rComp);
  inline DifferencePackedComponent(void * random, vector<int>& key);

  unsigned num_variables() const{
      if (hack_cachet){
        return 0;
      }
      if (hack_deleted) {
          return old_num_vars;
      }

    uint64_t *p = (uint64_t *) data_;
    return (*p >> bits_of_data_size()) & (uint64_t) variable_mask();

  }

  unsigned data_size() const {
      if (hack_cachet){
        return 8;
      }
      if (hack_deleted)
          return old_size;

         return *data_ & _data_size_mask;
    }

  unsigned data_only_byte_size() const {
        if (hack_cachet){
          return 8;
        }
        return data_size()* sizeof(unsigned);
    }

  unsigned raw_model_count_byte_size() const{
    return model_count_.get_mpf_t()->_mp_size * sizeof(mp_limb_t);
  }
  unsigned raw_data_byte_size() const {
          if (hack_cachet){
            return 8 + model_count_.get_mpf_t()->_mp_size * sizeof(mp_limb_t); 
          }
          return data_size()* sizeof(unsigned)
               + model_count_.get_mpf_t()->_mp_size * sizeof(mp_limb_t);
    }

    // raw data size with the overhead
    // for the supposed 16byte alignment of malloc
    unsigned sys_overhead_raw_data_byte_size() const {
      if (hack_cachet){
        return 0;
      }
      unsigned ds = data_size()* sizeof(unsigned);
      unsigned ms = model_count_.get_mpf_t()->_mp_size * sizeof(mp_limb_t);
//      unsigned mask = 0xfffffff8;
//      return (ds & mask) + ((ds & 7)?8:0)
//            +(ms & mask) + ((ms & 7)?8:0);
      unsigned mask = 0xfffffff0;
            return (ds & mask) + ((ds & 15)?16:0)
                  +(ms & mask) + ((ms & 15)?16:0);
    }

  bool equals(const uint64_t clhash_key, const uint64_t num_variables) const {
    // cout << "DIFFERENCE "<< num_variables << " "<< number_variables << endl;
    // cout << "KEY" << clhash_key << " "<< clhash_key_<< endl;
    return ((clhash_key == clhash_key_)); 
    // && (number_variables == num_variables));
  }

  bool equals(const uint64_t clhash_key) const {
    return (clhash_key == clhash_key_);
  }

  bool equals(const DifferencePackedComponent &comp) const {
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

  // bool equals(const DifferencePackedComponent &comp, const unsigned char* sha1_hashkey_dash) const {
  bool equals(const DifferencePackedComponent &comp, const uint64_t clhash_key) const {
    if(hashkey_ != comp.hashkey()){
      return false;
    }

    return (clhash_key == clhash_key_);
    //  bool we_return = true;
    //  for(unsigned i = 0; i < 20; i++) {
    //      if (sha1_hashkey_dash[i] != sha1_hashkey_[i]) {
    //          we_return = false;
    //          break;
    //      }
    //  }
    //  return we_return;

    //data_ is NULL!! Cannot check
    /*unsigned* p = data_;
    unsigned* r = comp.data_;
    while(p != data_ + data_size()) {
        if(*(p++) != *(r++)) {
            if (we_return == true) {
                cout << "UNLIKELY!!! false->true" << endl;
            }
            return false;
        }
    }

    if (we_return == false) {
        cout << "UNLIKELY!!! true->false" << endl;
    }
    return true;*/
  }

  // void compute_sha1_hash(unsigned char* sha1_hashkey_dash) const
  // {
  //     SHA_CTX ctx;
  //     int ret = SHA1_Init(&ctx);
  //     assert(ret == 1);
  //     ret = SHA1_Update(&ctx, (void*)data_, sizeof(unsigned)*data_size());
  //     assert(ret == 1);
  //     ret = SHA1_Final(sha1_hashkey_dash, &ctx);
  //     assert(ret == 1);
  // }
  
uint64_t getvar(){
  return 0;
  // return number_variables;
  }

uint64_t compute_cachetclhash(){
  return clhash_key_;
  }

uint64_t compute_clhash(){
  return clhash_key_;
  // clhasher h(random_key_); 
  // return h((void*)data_, sizeof(unsigned)*data_size());
  }
private:

};

DifferencePackedComponent::DifferencePackedComponent(unsigned int graphhash,unsigned int num_variables) {
    hack_cachet = true;
    clhash_key_ = graphhash;
    // number_variables = num_variables;
    delete[] data_;
    data_ = nullptr;
}


DifferencePackedComponent::DifferencePackedComponent(void * random, vector<int>& key) {
  hack_cachet = true;
  delete[] data_;
  data_ = nullptr;
  random_key_ = random;
  clhasher h(random);
  clhash_key_ = h(key.data(), key.size());
}


DifferencePackedComponent::DifferencePackedComponent(void * random,Component &rComp) {
  unsigned max_var_diff = 0;
  unsigned hashkey_vars = *rComp.varsBegin();
  for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++) {
    hashkey_vars = (hashkey_vars * 3) + *it;
    if ((*it - *(it - 1)) - 1 > max_var_diff)
      max_var_diff = (*it - *(it - 1)) - 1 ;
  }

  unsigned hashkey_clauses = *rComp.clsBegin();
  unsigned max_clause_diff = 0;
  if (*rComp.clsBegin()) {
    for (auto jt = rComp.clsBegin() + 1; *jt != clsSENTINEL; jt++) {
      hashkey_clauses = hashkey_clauses*3 + *jt;
      if (*jt - *(jt - 1) - 1 > max_clause_diff)
        max_clause_diff = *jt - *(jt - 1) - 1;
    }
  }

  hashkey_ = hashkey_vars + ((unsigned) hashkey_clauses << 11) + ((unsigned) hashkey_clauses >> 23);

  //VERIFIED the definition of bits_per_var_diff and bits_per_clause_diff
  unsigned bits_per_var_diff = log2(max_var_diff) + 1;
  unsigned bits_per_clause_diff = log2(max_clause_diff) + 1;

  assert(bits_per_var_diff <= 31);
  assert(bits_per_clause_diff <= 31);

  unsigned data_size_vars = bits_of_data_size() + 2*bits_per_variable() + 5;

  data_size_vars += (rComp.num_variables() - 1) * bits_per_var_diff ;
  unsigned data_size_clauses = 0;
  if(*rComp.clsBegin())
    data_size_clauses += bits_per_clause() + 5
       + (rComp.numLongClauses() - 1) * bits_per_clause_diff;

  unsigned data_size = (data_size_vars + data_size_clauses)/bits_per_block();
    data_size+=  ((data_size_vars + data_size_clauses) % bits_per_block())? 1 : 0;

  data_ = new unsigned[data_size];
  assert((data_size >> bits_of_data_size()) == 0);
  BitStuffer<unsigned> bs(data_);

  bs.stuff(data_size, bits_of_data_size());
  bs.stuff(rComp.num_variables(), bits_per_variable());
  bs.stuff(bits_per_var_diff, 5);
  bs.stuff(*rComp.varsBegin(), bits_per_variable());

  if(bits_per_var_diff)
  for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++)
    bs.stuff(*it - *(it - 1) - 1, bits_per_var_diff);


  if (*rComp.clsBegin()) {
    bs.stuff(bits_per_clause_diff, 5);
    bs.stuff(*rComp.clsBegin(), bits_per_clause());
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

  // SHA_CTX ctx;
  // int ret = SHA1_Init(&ctx);
  // assert(ret == 1);
  // ret = SHA1_Update(&ctx, (void*)data_, sizeof(unsigned)*data_size);
  // assert(ret == 1);
  // ret = SHA1_Final(sha1_hashkey_, &ctx);
  // assert(ret == 1);
  // clhasher h(UINT64_C(0x23a23cf5033c3c81),UINT64_C(0xb3816f6a2c68e530));
  random_key_ = random;
  // memcpy (random_key_, random, sizeof(random));
  clhasher h(random);
  clhash_key_ = h((void*)data_, sizeof(unsigned)*data_size);
  delete[] data_;
  data_ = nullptr;
  hack_cachet = true;
  // clhash_key_ = clhash(random_,(void*)data_, sizeof(unsigned)*data_size);

}

DifferencePackedComponent::DifferencePackedComponent(Component &rComp) {

  unsigned max_var_diff = 0;
  unsigned hashkey_vars = *rComp.varsBegin();
  for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++) {
    hashkey_vars = (hashkey_vars * 3) + *it;
    if ((*it - *(it - 1)) - 1 > max_var_diff)
      max_var_diff = (*it - *(it - 1)) - 1 ;
  }

  unsigned hashkey_clauses = *rComp.clsBegin();
  unsigned max_clause_diff = 0;
  if (*rComp.clsBegin()) {
    for (auto jt = rComp.clsBegin() + 1; *jt != clsSENTINEL; jt++) {
      hashkey_clauses = hashkey_clauses*3 + *jt;
      if (*jt - *(jt - 1) - 1 > max_clause_diff)
        max_clause_diff = *jt - *(jt - 1) - 1;
    }
  }

  hashkey_ = hashkey_vars + ((unsigned) hashkey_clauses << 11) + ((unsigned) hashkey_clauses >> 23);

  //VERIFIED the definition of bits_per_var_diff and bits_per_clause_diff
  unsigned bits_per_var_diff = log2(max_var_diff) + 1;
  unsigned bits_per_clause_diff = log2(max_clause_diff) + 1;

  assert(bits_per_var_diff <= 31);
  assert(bits_per_clause_diff <= 31);

  unsigned data_size_vars = bits_of_data_size() + 2*bits_per_variable() + 5;

  data_size_vars += (rComp.num_variables() - 1) * bits_per_var_diff ;

  unsigned data_size_clauses = 0;
  if(*rComp.clsBegin())
    data_size_clauses += bits_per_clause() + 5
       + (rComp.numLongClauses() - 1) * bits_per_clause_diff;

  unsigned data_size = (data_size_vars + data_size_clauses)/bits_per_block();
    data_size+=  ((data_size_vars + data_size_clauses) % bits_per_block())? 1 : 0;

  data_ = new unsigned[data_size];

  assert((data_size >> bits_of_data_size()) == 0);
  BitStuffer<unsigned> bs(data_);

  bs.stuff(data_size, bits_of_data_size());
  bs.stuff(rComp.num_variables(), bits_per_variable());
  bs.stuff(bits_per_var_diff, 5);
  bs.stuff(*rComp.varsBegin(), bits_per_variable());

  if(bits_per_var_diff)
  for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++)
    bs.stuff(*it - *(it - 1) - 1, bits_per_var_diff);


  if (*rComp.clsBegin()) {
    bs.stuff(bits_per_clause_diff, 5);
    bs.stuff(*rComp.clsBegin(), bits_per_clause());
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
}
//DifferencePackedComponent::DifferencePackedComponent(Component &rComp) {
//
//  unsigned max_var_diff = 0;
//  unsigned hashkey_vars = *rComp.varsBegin();
//  for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++) {
//    hashkey_vars = (hashkey_vars * 3) + *it;
//    if ((*it - *(it - 1)) > max_var_diff)
//      max_var_diff = (*it - *(it - 1)) ;
//  }
//
//  unsigned hashkey_clauses = *rComp.clsBegin();
//  unsigned max_clause_diff = 0;
//  if (*rComp.clsBegin()) {
//    for (auto jt = rComp.clsBegin() + 1; *jt != clsSENTINEL; jt++) {
//      hashkey_clauses = hashkey_clauses * 3 + *jt;
//      if (*jt - *(jt - 1) > max_clause_diff)
//        max_clause_diff = *jt - *(jt - 1);
//    }
//  }
//
//  hashkey_ = hashkey_vars + (((unsigned) hashkey_clauses) << 16);
//
//  //VERIFIED the definition of bits_per_var_diff and bits_per_clause_diff
//  unsigned bits_per_var_diff = log2(max_var_diff) + 1;
//  unsigned bits_per_clause_diff = log2(max_clause_diff) + 1;
//
//  unsigned data_size_vars = 2*bits_per_variable() + 5;
//
//  data_size_vars += (rComp.num_variables() - 1) * bits_per_var_diff ;
//
//  unsigned data_size_clauses = bits_per_clause();
//  if(*rComp.clsBegin())
//    data_size_clauses += bits_per_clause() + 5
//       + (rComp.numLongClauses() - 1) * bits_per_clause_diff;
//
//  unsigned data_size = (data_size_vars + data_size_clauses)/bits_per_block();
//    data_size+=  ((data_size_vars + data_size_clauses) % bits_per_block())? 1 : 0;
//
//  data_ = new unsigned[data_size];
//
//  BitStuffer<unsigned> bs(data_);
//
//  bs.stuff(rComp.num_variables(), bits_per_variable());
//  bs.stuff(bits_per_var_diff, 5);
//  bs.stuff(*rComp.varsBegin(), bits_per_variable());
//
//  for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++)
//    bs.stuff(*it - *(it - 1), bits_per_var_diff);
//
//  if (*rComp.clsBegin()) {
//    bs.stuff(bits_per_clause_diff, 5);
//    bs.stuff(*rComp.clsBegin(), bits_per_clause());
//    for (auto jt = rComp.clsBegin() + 1; *jt != clsSENTINEL; jt++)
//      bs.stuff(*jt - *(jt - 1), bits_per_clause_diff);
//  }
//
//  // to check wheter the "END" block of bits_per_clause()
//  // many zeros fits into the current
//  bs.end_check(bits_per_clause());
//  // this will tell us if we computed the data_size
//  // correctly
//  bs.assert_size(data_size);
//}

#endif /* DIFFERENCE_PACKED_COMPONENT_H_ */
