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

#include "chibihash64.h"
#include "common.hpp"
#include "comp.hpp"

namespace GanakInt {

class HashedComp  {
public:
  HashedComp() = default;
  HashedComp(const HashedComp&) = default;
  HashedComp& operator=(const HashedComp&) = default;
  HashedComp(HashedComp&&) noexcept = default;
  static uint64_t set_comp(const Comp& comp, const uint64_t hash_seed, const BPCSizes& /*bpc*/) {
    return chibihash64(comp.get_raw_data(), comp.get_size()*4, hash_seed);
  }
  bool equals(const HashedComp&) const {
    return true;
  }
  uint64_t comp_bytes() const { return 0; }
  void set_free() { }
};

}
