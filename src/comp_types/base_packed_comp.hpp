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
#include <cstdint>
#include <gmpxx.h>
#include <iostream>
#include "../common.hpp"
#include "../mpreal.h"

using std::cout;

namespace GanakInt {

template<typename T>
class BaseComp {
public:
  uint32_t last_used_time() const { return last_used_time_; }
  const T& model_count() const { return *model_count_; }
  uint32_t bignum_bytes() const{
    if (!model_count_) return 0;
    return mp_data_size();
  }

  // raw data size with the overhead
  // for the supposed 16byte alignment of malloc
  uint32_t mp_data_size() const;

  // These below are to help us erase better from the cache
  void set_last_used_time(uint32_t time) { last_used_time_ = time; }
  void avg_last_used_time(uint32_t time, uint32_t div) {
    assert(time >= last_used_time_);
    last_used_time_ += (time-last_used_time_)/div;
  }

  void set_model_count(const T& rn) {
    assert(model_count_ == nullptr);
    model_count_ = new T(rn);
  }

  uint32_t get_hashkey() const  { return (uint32_t)clhashkey_; }
  bool model_count_found(){ return model_count_ != nullptr; }

  // a cache entry is deletable
  // only if it is not connected to an active
  // comp in the comp stack
  bool is_deletable() const { return delete_permitted; }
  void set_deletable() { delete_permitted = true; }

  void clear() {
    // before deleting the contents of this comp,
    // we should make sure that this comp is not present in the comp stack anymore!
    SLOW_DEBUG_DO(assert(is_deletable()));
  }

protected:
  uint64_t clhashkey_;
  T* model_count_ = nullptr;
  uint32_t last_used_time_:31 = 1; //effectively the score
  uint32_t delete_permitted:1 = false;
};

inline auto helper(mpz_ptr val) {
  return val->_mp_alloc * sizeof(mp_limb_t);
}

template<>
inline uint32_t BaseComp<mpz_class>::mp_data_size() const {
  return sizeof(mpz_class) + helper(model_count_->get_mpz_t());
}

template<>
inline uint32_t BaseComp<mpq_class>::mp_data_size() const {
  uint32_t a = helper(model_count_->get_den_mpz_t());
  uint32_t b = helper(model_count_->get_num_mpz_t());
  return sizeof(mpz_class) + a + b;
}

template<>
inline uint32_t BaseComp<mpfr::mpreal>::mp_data_size() const {
  unsigned long prec = mpfr_get_prec(model_count_->mpfr_ptr());
  unsigned long num_limbs = (prec + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
  size_t memory_usage = sizeof(mpfr::mpreal) + num_limbs * sizeof(mp_limb_t);
  return memory_usage;
}

}
