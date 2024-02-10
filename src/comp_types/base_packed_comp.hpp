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

#include <cassert>
#include <gmpxx.h>
#include <iostream>
#include "common.hpp"

using std::cout;

struct BitStufferReader {
  BitStufferReader(uint32_t* _data) {data = _data;}
  uint32_t read_bits(uint32_t bits) {
    assert(bits <= 32);
    uint32_t ret = 0;

    uint32_t byte = data[at_byte];
    byte >>= at_bit % 32;
    ret += byte;
    uint32_t bits_read = 32-(at_bit)%32;
    // ret has 32-(at_bit%32) bits
    if (bits_read > bits) {
      // cut to BITS
      ret <<= (32-bits);
      ret >>= (32-bits);
      at_bit += bits;
      return ret;
    }
    at_bit += bits_read;
    assert(at_bit % 32 == 0);
    at_byte++;
    ret += read_bits(bits-bits_read) << bits_read;
    return ret;
  }

  uint32_t* data;
  uint32_t at_bit = 0;
  uint32_t at_byte = 0;
};

template <class T>
 class BitStuffer {
 public:
  BitStuffer(T *data):data_start_(data),p(data){
    *p = 0;
  }

  void stuff(const uint32_t val, const uint32_t num_bits_val){
      assert(num_bits_val > 0);
      assert((val >> num_bits_val) == 0);
      if(end_of_bits_ == 0) *p = 0;
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

  const uint32_t _bits_per_block = (sizeof(T) << 3);

};


class BasePackedComponent {
public:
  static uint32_t log2(uint32_t v) {
         // taken from
         // http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup
         const signed char LogTable256[256] =
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

  uint32_t creation_time() const {
    return creation_time_;
  }

  const mpz_class &model_count() const {
    return *model_count_;
  }

  uint32_t alloc_of_model_count() const{
    if (!model_count_) return 0;
    return sizeof(mpz_class)+sys_overhead_raw_data_byte_size();
  }

  // raw data size with the overhead
  // for the supposed 16byte alignment of malloc
  uint32_t sys_overhead_raw_data_byte_size() const {
    uint32_t ds;
    ds = 0;
    uint32_t ms = model_count_->get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
    uint32_t mask = 0xfffffff0;
    return (ds & mask)+((ds & 15)?16:0) + (ms & mask)+((ms & 15)?16:0);
  }

  void set_creation_time(uint32_t time) {
    creation_time_ = time;
  }

  void set_model_count(const mpz_class &rn, uint32_t time) {
    assert(model_count_ == NULL);
    model_count_ = new mpz_class(rn);
    length_solution_period_and_flags_ = (time - creation_time_) | (length_solution_period_and_flags_ & 1);
  }

  uint32_t get_hashkey() const  { return hashkey_; }
  bool modelCountFound(){
    return (length_solution_period_and_flags_ >> 1);
  }

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
    SLOW_DEBUG_DO(assert(isDeletable()));
  }

protected:
  // data_ contains in packed form the variable indices
  // and clause indices of the comp ordered
  // structure is
  // var var ... clause clause ...
  // clauses begin at clauses_ofs_

  uint64_t clhashkey_;
  uint32_t hashkey_ = 0;

  mpz_class* model_count_ = NULL;
  uint32_t creation_time_ = 1;


  // this is:  length_solution_period = length_solution_period_and_flags_ >> 1
  // length_solution_period == 0 means unsolved
  // and the first bit is "delete_permitted"
  uint32_t length_solution_period_and_flags_ = 0;

  // deletion is permitted only after
  // the copy of this comp in the stack
  // does not exist anymore
};
