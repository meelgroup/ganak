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

#include <functional>
#include <set>
#include <iostream>

#include "base_packed_comp.hpp"
#include "comp.hpp"
#include "../clhash/clhash.h"
#include "primitive_types.hpp"
#include "structures.hpp"

class HashedComp: public BaseComp {
public:
  HashedComp() = default;
  inline HashedComp(void* hash_seed, Component &r_comp, uint32_t* tmp_data);
  uint32_t raw_data_byte_size() const {
    return BaseComp::alloc_of_model_count();
  }

  uint64_t get_clhashkey() const { return clhashkey_; }
  bool equals_clhashkey(const HashedComp &comp) const {
    return clhashkey_ == comp.get_clhashkey();
  }
};

HashedComp::HashedComp(
    void* hash_seed, Component &rComp, uint32_t* tmp_data) {
  auto data = tmp_data;
  uint32_t at = 0;

  data[at++] = rComp.nVars();
  if (*rComp.vars_begin()) {
    for (auto it = rComp.vars_begin(); *it != sentinel; it++) data[at++] = *it;
  }
  data[at++] = sentinel;

  if (*rComp.cls_begin()) {
    for (auto jt = rComp.cls_begin(); *jt != sentinel; jt++) data[at++] =*jt;
  }
  data[at] = at;
  at++;

  clhasher h(hash_seed);
  clhashkey_ = h(data, at);
}
