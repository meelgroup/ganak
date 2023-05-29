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

#include "base_packed_comp.h"
#include "comp.h"
#include "../clhash/clhash.h"
/* #include "../clhash/minim.h" */
#include "structures.h"

class DifferencePackedComponent: public BasePackedComponent {
public:
  DifferencePackedComponent() { }
  inline DifferencePackedComponent(void* randomseedforCLHASH, Component &rComp, const BPCSizes& sz, uint32_t* tmp_data);
  uint32_t raw_data_byte_size() const {
    return BasePackedComponent::alloc_of_model_count();
  }


  uint64_t get_clhashkey() const { return clhashkey_; }
  bool equals_clhashkey(const DifferencePackedComponent &comp) const {
    if (hashkey_ != comp.get_hashkey()) return false;
    if (clhashkey_ != comp.get_clhashkey()) return false;
    return true;
  }
};

DifferencePackedComponent::DifferencePackedComponent(void* randomseedforCLHASH, Component &rComp, const BPCSizes& sz, uint32_t* tmp_data) {
  // first, generate hashkey, and compute max diff for cls and vars
  uint32_t max_var_diff = 0;
  uint32_t hashkey_vars = *rComp.varsBegin();
  for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++) {
    hashkey_vars = (hashkey_vars * 3) + *it;
    if ((*it - *(it - 1)) - 1 > max_var_diff)
      max_var_diff = (*it - *(it - 1)) - 1 ;
  }

  uint32_t hashkey_clauses = *rComp.clsBegin();
  uint32_t max_clause_diff = 0;
  if (*rComp.clsBegin()) {
    for (auto jt = rComp.clsBegin() + 1; *jt != clsSENTINEL; jt++) {
      hashkey_clauses = hashkey_clauses*3 + *jt;
      if (*jt - *(jt - 1) - 1 > max_clause_diff)
        max_clause_diff = *jt - *(jt - 1) - 1;
    }
  }

  hashkey_ = hashkey_vars + ((uint32_t) hashkey_clauses << 11) +
    ((uint32_t) hashkey_clauses >> 23);
  uint32_t bits_per_var_diff = log2(max_var_diff) + 1;
  uint32_t bits_per_clause_diff = log2(max_clause_diff) + 1;

  assert(bits_per_var_diff <= 31);
  assert(bits_per_clause_diff <= 31);

  uint32_t data_size_vars = sz.bits_of_data_size + 2*sz.bits_per_variable + 5;

  data_size_vars += (rComp.nVars() - 1) * bits_per_var_diff ;
  uint32_t data_size_clauses = 0;
  if(*rComp.clsBegin())
    data_size_clauses += sz.bits_per_clause + 5
       + (rComp.numLongClauses() - 1) * bits_per_clause_diff;

  uint32_t data_size = (data_size_vars + data_size_clauses)/sz.bits_per_block;
  data_size += ((data_size_vars + data_size_clauses) % sz.bits_per_block)? 1 : 0;

  auto data_ = tmp_data;
  assert((data_size >> sz.bits_of_data_size) == 0);
  BitStuffer<uint32_t> bs(data_);

  bs.stuff(data_size, sz.bits_of_data_size);
  bs.stuff(rComp.nVars(), sz.bits_per_variable);
  bs.stuff(bits_per_var_diff, 5);
  bs.stuff(*rComp.varsBegin(), sz.bits_per_variable);

  if(bits_per_var_diff)
    for (auto it = rComp.varsBegin() + 1; *it != varsSENTINEL; it++)
      bs.stuff(*it - *(it - 1) - 1, bits_per_var_diff);

  if (*rComp.clsBegin()) {
    bs.stuff(bits_per_clause_diff, 5);
    bs.stuff(*rComp.clsBegin(), sz.bits_per_clause);
    if(bits_per_clause_diff)
     for (auto jt = rComp.clsBegin() + 1; *jt != clsSENTINEL; jt++)
      bs.stuff(*jt - *(jt - 1) - 1, bits_per_clause_diff);
  }

  // to check wheter the "END" block of bits_per_clause()
  // many zeros fits into the current
  //bs.end_check(bits_per_clause());
  // this will tell us if we computed the data_size
  // correctly
  bs.assert_size(data_size);

  clhasher h(randomseedforCLHASH);
  clhashkey_ = h(data_, data_size);
}
