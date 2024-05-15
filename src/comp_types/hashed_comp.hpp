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

template<typename T>
class HashedComp: public BaseComp<T> {
public:
  HashedComp() = default;
  inline HashedComp(void* hash_seed, const Comp& r_comp);
  uint32_t bignum_bytes() const { return BaseComp<T>::bignum_bytes(); }

  uint64_t get_clhashkey() const { return BaseComp<T>::clhashkey_; }
  bool equals_clhashkey(const HashedComp &comp) const {
    return BaseComp<T>::clhashkey_ == comp.get_clhashkey();
  }
};

template<typename T>
HashedComp<T>::HashedComp(void* hash_seed, const Comp& comp) {
  clhasher h(hash_seed);
  BaseComp<T>::clhashkey_ = h(comp.get_raw_data().data(), comp.get_raw_data().size());
}
