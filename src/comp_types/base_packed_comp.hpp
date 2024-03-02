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
#include "../common.hpp"

using std::cout;

class BaseComp {
public:
  uint32_t last_used_time() const { return last_used_time_; }
  const mpz_class &model_count() const { return *model_count_; }
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

  void set_last_used_time(uint32_t time) { last_used_time_ = time; }
  void set_model_count(const mpz_class &rn) {
    assert(model_count_ == nullptr);
    model_count_ = new mpz_class(rn);
    delete_permitted = 1;
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
    SLOW_DEBUG_DO(assert(isDeletable()));
  }

protected:
  uint64_t clhashkey_;
  mpz_class* model_count_ = nullptr;
  uint32_t last_used_time_ = 1; //effectively the score

  // deletion is permitted only after
  // the copy of this comp in the stack
  // does not exist anymore
  bool delete_permitted = false;
};
