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

#include "base_packed_comp.hpp"
#include "comp.hpp"
#include "../clhash/clhash.h"

namespace GanakInt {

class HashedComp: public BaseComp {
public:
  HashedComp() = default;
  HashedComp(const HashedComp&) = default;
  HashedComp& operator=(const HashedComp&) = default;
  HashedComp(void* hash_seed, const Comp& comp) {
    clhasher h(hash_seed);
    clhashkey_ = h(comp.get_raw_data(), comp.get_size());
    model_count_ = nullptr;
  }
  uint64_t bignum_bytes() const { return BaseComp::bignum_bytes(); }

  uint64_t get_clhashkey() const { return BaseComp::clhashkey_; }
  bool equals_clhashkey(const HashedComp &comp) const {
    return BaseComp::clhashkey_ == comp.get_clhashkey();
  }
};

}
