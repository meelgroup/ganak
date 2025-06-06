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
#include "chibihash64.h"
#include "common.hpp"
#include "comp.hpp"
#include "containers.hpp"
#include <vector>
using std::vector;

namespace GanakInt {

class HashedComp: public BaseComp {
public:
  HashedComp() = default;
  HashedComp(const HashedComp&) = default;
  HashedComp& operator=(const HashedComp&) = default;
  HashedComp(uint64_t hash_seed, const Comp& comp, vector<Lit*> long_cls,
      const LiteralIndexedVector<TriValue>& vals) {
    static vector<uint64_t> d;
    static vector<Lit> tmp;

    d.clear();
    d.reserve(comp.get_size());
    all_cls_in_comp(comp, cl_id) {
      tmp.clear();
      Lit* cl = long_cls[*cl_id];
      for(auto& l = cl; *l != SENTINEL_LIT; l++) {
        assert(vals[*l] != T_TRI); //cannot be true
        if (vals[*l] == F_TRI) continue; // skip falsified lits
        tmp.push_back(*l);
      }
      uint64_t c = chibihash64(tmp.data(), tmp.size()*sizeof(Lit), hash_seed);
      d.push_back(c);
    }
    auto sz = d.size();
    std::sort(d.begin(), d.end());
    auto last = std::unique(d.begin(), d.end());
    d.erase(last, d.end());
    if (d.size() < sz) {
      /* cerr << "WARNING: duplicate clauses in component, size reduced from " << sz */
      /*      << " to " << d.size() << endl; */
    }
    d.push_back(sentinel);
    all_vars_in_comp(comp, v) d.push_back(*v);

    clhashkey_ = chibihash64(d.data(), d.size()*sizeof(uint64_t), hash_seed);
    model_count_ = nullptr;
  }
  uint64_t bignum_bytes() const { return BaseComp::bignum_bytes(); }

  uint64_t get_clhashkey() const { return BaseComp::clhashkey_; }
  bool equals_clhashkey(const HashedComp &comp) const {
    return BaseComp::clhashkey_ == comp.get_clhashkey();
  }
};

}
