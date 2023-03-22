/*
 * base_packed_comp.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#ifndef BASE_PACKED_COMPONENT_H_
#define BASE_PACKED_COMPONENT_H_

#include <assert.h>
#include <gmpxx.h>
#include <iostream>

using std::cout;

template <class T>
 class BitStuffer {
 public:
  BitStuffer(T *data):data_start_(data),p(data){
    *p = 0;
  }

  void stuff(const uint32_t val, const uint32_t num_bits_val){
      assert(num_bits_val > 0);
      assert((val >> num_bits_val) == 0);
      if(end_of_bits_ == 0)
        *p = 0;
      assert((*p >> end_of_bits_) == 0);
      *p |= val << end_of_bits_;
      end_of_bits_ += num_bits_val;
      if (end_of_bits_ > _bits_per_block){
        //assert(*p);
        end_of_bits_ -= _bits_per_block;
        *(++p) = val >> (num_bits_val - end_of_bits_);
        assert(!(end_of_bits_ == 0) | (*p == 0));
      }
      else if (end_of_bits_ == _bits_per_block){
        end_of_bits_ -= _bits_per_block;
        p++;
      }
  }

  void assert_size(uint32_t size){
    if(end_of_bits_ == 0)
       p--;
    assert(p - data_start_ == size - 1);
  }

 private:
  T *data_start_ = nullptr;
  T *p = nullptr;
  // in the current block
  // the bit postion just after the last bit written
  uint32_t end_of_bits_ = 0;

  static const uint32_t _bits_per_block = (sizeof(T) << 3);

};


class BasePackedComponent {
public:


  void set_hacked(const uint64_t _old_size, const uint64_t _old_num_vars) {
    old_size = _old_size;
    old_num_vars = _old_num_vars;
    delete[] data_;
    data_ = nullptr;
  }

  bool get_hacked(){
    if(old_size)
      return true;
    return false;
  }

  static uint32_t bits_per_variable() {
    return _bits_per_variable;
  }
  static uint32_t variable_mask() {
      return _variable_mask;
  }
  static uint32_t bits_per_clause() {
    return _bits_per_clause;
  }

  static uint32_t bits_per_block(){
	  return _bits_per_block;
  }

  static uint32_t bits_of_data_size(){
    return _bits_of_data_size;
  }

  static void adjustPackSize(uint32_t maxVarId, uint32_t maxClId);

  BasePackedComponent() {}
  BasePackedComponent(uint32_t creation_time): creation_time_(creation_time) {}

  ~BasePackedComponent() {
    if (data_){
      delete [] data_;
    }
#ifdef DOPCC
    if (clhashkey_) delete [] clhashkey_;
#endif
  }
  static void outbit(uint32_t v){
   for(auto i=0; i<32;i++){
     cout << ((v&2147483648)?"1":"0");
      v&=2147483648-1;
      v <<= 1;
    }
  }


  static uint32_t log2(uint32_t v){
         // taken from
         // http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup
         static const signed char LogTable256[256] =
         {
         #define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
             -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
             LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
             LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
         };

         uint32_t r;     // r will be lg(v)
         uint32_t t, tt; // temporaries

         if ((tt = (v >> 16)))
         {
           r = (t = (tt >> 8)) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
         }
         else
         {
           r = (t = (v >> 8)) ? 8 + LogTable256[t] : LogTable256[v];
         }
         return r;
  }

  uint32_t creation_time() {
    return creation_time_;
  }

  const mpz_class &model_count() const {
    return model_count_;
  }

  uint32_t alloc_of_model_count() const{
        return sizeof(mpz_class)
               + model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
  }

  void set_creation_time(uint32_t time) {
    creation_time_ = time;
  }

  void set_model_count(const mpz_class &rn, uint32_t time) {
    model_count_ = rn;
    length_solution_period_and_flags_ = (time - creation_time_) | (length_solution_period_and_flags_ & 1);
  }

  uint32_t hashkey() const  {
    return hashkey_;
  }

  bool modelCountFound(){
    return (length_solution_period_and_flags_ >> 1);
  }

 // inline bool equals(const BasePackedComponent &comp) const;

  // a cache entry is deletable
  // only if it is not connected to an active
  // comp in the comp stack
  bool isDeletable() const {
    return length_solution_period_and_flags_ & 1;
  }
  void set_deletable() {
    length_solution_period_and_flags_ |= 1;
  }

  void clear() {
    // before deleting the contents of this comp,
    // we should make sure that this comp is not present in the comp stack anymore!
    assert(isDeletable());
    if (data_) delete [] data_;
    data_ = nullptr;
#ifdef DOPCC
    if (clhashkey_) delete [] clhashkey_;
    clhashkey_ = nullptr;
#endif
  }

  static uint32_t _debug_static_val;

protected:
  // data_ contains in packed form the variable indices
  // and clause indices of the comp ordered
  // structure is
  // var var ... clause clause ...
  // clauses begin at clauses_ofs_
  uint32_t* data_ = nullptr;

#ifdef DOPCC
  uint64_t* clhashkey_ = nullptr;
#endif
  // vector <uint64_t> clhash_key_;
  uint32_t hashkey_ = 0;

  mpz_class model_count_;

  uint32_t creation_time_ = 1;
  uint32_t hack_ = 0;
  uint32_t old_size = 0;
  uint32_t old_num_vars = 0;


  // this is:  length_solution_period = length_solution_period_and_flags_ >> 1
  // length_solution_period == 0 means unsolved
  // and the first bit is "delete_permitted"
  uint32_t length_solution_period_and_flags_ = 0;

  // deletion is permitted only after
  // the copy of this comp in the stack
  // does not exist anymore


protected:
  static uint32_t _bits_per_clause, _bits_per_variable; // bitsperentry
  static uint32_t _bits_of_data_size; // number of bits needed to store the data size
  static uint32_t _data_size_mask;
  static uint32_t _variable_mask, _clause_mask;
  static const uint32_t _bits_per_block= (sizeof(uint32_t) << 3);

};

#endif /* BASE_PACKED_COMPONENT_H_ */
