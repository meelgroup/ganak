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

using std::cout;

namespace GanakInt {

class BaseComp {
public:
  BaseComp() = default;
  BaseComp& operator=(const BaseComp& b) {
    if (b.model_count_) model_count_ = b.model_count_->dup();
    else model_count_ = nullptr;
    last_used_time_ = b.last_used_time_;
    delete_permitted = b.delete_permitted;
    return *this;
  }
  BaseComp(BaseComp&& b) noexcept
    : model_count_(std::move(b.model_count_)),
      last_used_time_(b.last_used_time_),
      delete_permitted(b.delete_permitted)
  {
    b.last_used_time_ = 1;
    b.delete_permitted = false;
  }
  BaseComp(const BaseComp& b) {
    if (b.model_count_) model_count_ = b.model_count_->dup();
    else model_count_ = nullptr;
    last_used_time_ = b.last_used_time_;
    delete_permitted = b.delete_permitted;
  }
  uint32_t last_used_time() const { return last_used_time_; }
  const FF& model_count() const { return model_count_; }
  uint64_t bignum_bytes() const{
    if (!model_count_) return 0;
    return model_count_->bytes_used();
  }

  // These below are to help us erase better from the cache
  void set_last_used_time(uint32_t time) { last_used_time_ = time; }
  void avg_last_used_time(uint32_t time, uint32_t div) {
    assert(time >= last_used_time_);
    last_used_time_ += (time-last_used_time_)/div;
  }

  void set_model_count(const FF& rn) {
    assert(model_count_ == nullptr);
    model_count_ = rn->dup();
  }

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
  FF model_count_ = nullptr;
  uint32_t last_used_time_:31 = 1; //effectively the score
  uint32_t delete_permitted:1 = false;
};

}
