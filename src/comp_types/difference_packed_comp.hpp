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

#include "common.hpp"
#include "comp.hpp"
#include "chibihash64.h"

namespace GanakInt {

class DiffPackedComp  {
  class BitStuffer {
  public:
    BitStuffer(uint32_t *_data):data(_data),p(data){ *p = 0; }

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

    void assert_size(const uint32_t size){
      if(end_of_bits_ == 0) p--;
      assert(p - data == size - 1);
    }

  private:
    uint32_t *data = nullptr;
    uint32_t *p = nullptr;
    uint32_t end_of_bits_ = 0; // in the current block, the bit postion just after the last bit written
    const uint32_t _bits_per_block = (sizeof(uint32_t) << 3);
  };


public:
  DiffPackedComp() = default;
  ~DiffPackedComp() { delete[] data; }
  DiffPackedComp& operator=(const DiffPackedComp& other) noexcept{
    if (data != nullptr) delete[] data;
    data_size = other.data_size;
    data = new uint32_t[data_size];
    std::memcpy(data, other.data, data_size * sizeof(uint32_t));
    return *this;
  }
  DiffPackedComp(const DiffPackedComp& other) {
    data_size = other.data_size;
    data = new uint32_t[data_size];
    std::memcpy(data, other.data, data_size * sizeof(uint32_t));
  }
  DiffPackedComp(DiffPackedComp&& other) noexcept : data(other.data), data_size(other.data_size) {
    other.data = nullptr;
    other.data_size = 0;
  }

  uint64_t comp_bytes() const {
    return data_size * sizeof(uint32_t);
  }

  void set_free() {
    delete[] data;
    data = nullptr;
    data_size = 0;
  }

  bool equals(const DiffPackedComp &other) const {
    SLOW_DEBUG_DO(assert(data != nullptr));
    SLOW_DEBUG_DO(assert(other.data != nullptr));

    if (data_size != other.data_size) return false;
    return std::memcmp(data, other.data, data_size * sizeof(uint32_t)) == 0;
  }

  uint64_t set_comp(const Comp &comp, const uint64_t hash_seed, const BPCSizes& sz) {
    // first, generate hashkey, and compute max diff for cls and vars
    uint32_t max_var_diff = 0;
    uint32_t v = *comp.vars_begin();
    for (auto it = comp.vars_begin() + 1; *it != sentinel; it++) {
      v = (v * 3) + *it;
      if ((*it - *(it - 1)) - 1 > max_var_diff)
        max_var_diff = (*it - *(it - 1)) - 1 ;
    }

    uint32_t cl = *comp.cls_begin();
    uint32_t max_clause_diff = 0;
    if (*comp.cls_begin() != sentinel) {
      for (auto it = comp.cls_begin() + 1; *it != sentinel; it++) {
        cl = cl*3 + *it;
        if (*it - *(it - 1) - 1 > max_clause_diff)
          max_clause_diff = *it - *(it - 1) - 1;
      }
    }

    uint32_t bits_per_var_diff = mlog2(max_var_diff) + 1;
    uint32_t bits_per_clause_diff = mlog2(max_clause_diff) + 1;

    assert(bits_per_var_diff <= 31);
    assert(bits_per_clause_diff <= 31);

    uint32_t data_size_vars = sz.bits_of_data_size + 2*sz.bits_per_variable + 5;

    data_size_vars += (comp.nVars() - 1) * bits_per_var_diff ;
    uint32_t data_size_clauses = 0;
    if(*comp.cls_begin() != sentinel) {
      data_size_clauses += sz.bits_per_clause + 5
         + (comp.num_long_cls() - 1) * bits_per_clause_diff;
    }

    data_size = (data_size_vars + data_size_clauses)/sz.bits_per_block;
    data_size += ((data_size_vars + data_size_clauses) % sz.bits_per_block)? 1 : 0;
    data = new uint32_t[data_size];
    assert((data_size >> sz.bits_of_data_size) == 0);
    BitStuffer bs(data);

    bs.stuff(data_size, sz.bits_of_data_size);
    bs.stuff(comp.nVars(), sz.bits_per_variable);
    bs.stuff(bits_per_var_diff, 5);
    bs.stuff(*comp.vars_begin(), sz.bits_per_variable);

    if(bits_per_var_diff) {
      for (auto it = comp.vars_begin() + 1; *it != sentinel; it++)
        bs.stuff(*it - *(it - 1) - 1, bits_per_var_diff);
    }

    if (*comp.cls_begin()) {
      bs.stuff(bits_per_clause_diff, 5);
      bs.stuff(*comp.cls_begin(), sz.bits_per_clause);
      if(bits_per_clause_diff) {
        for (auto jt = comp.cls_begin() + 1; *jt != sentinel; jt++)
          bs.stuff(*jt - *(jt - 1) - 1, bits_per_clause_diff);
      }
    }

    // To check whether the "END" block of bits_per_clause()
    // many zeros fits into the current bs.end_check(bits_per_clause());
    // This will tell us if we computed the data_size correctly
    bs.assert_size(data_size);

    return chibihash64(comp.get_raw_data(), comp.get_size()*4, hash_seed);
  }

private:
  uint32_t* data = nullptr; // the packed data
  uint32_t data_size = 0; // the size of the packed data in bytes
};

}
