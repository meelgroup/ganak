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
#include <vector>
#include "../common.hpp"
#include "structures.hpp"

using std::vector;

namespace GanakInt {

class Comp {
public:
  Comp() = delete;

  void set_id(CacheEntryID id) { id_ = id; }
  CacheEntryID id() const { return id_; }

  uint32_t* vs_cls_data() {
    uint32_t* ptr = (uint32_t*)this;
    ptr += sizeof(Comp)/sizeof(uint32_t);
    return ptr;
  }

  uint32_t const* vs_cls_data() const {
    uint32_t* ptr = (uint32_t*)this;
    ptr += sizeof(Comp)/sizeof(uint32_t);
    return ptr;
  }

  inline void add_var(const uint32_t var) {
    // the only time a sentinel is added should be in a
    // call to closeVarData(..)
    SLOW_DEBUG_DO(assert(var != sentinel));
    vs_cls_data()[size++] = var;
  }

  inline void close_vars_data() {
    vs_cls_data()[size++] = sentinel;
    clauses_offs = size;
  }

  inline void add_cl(const ClauseIndex cl) {
    // the only time a sentinel is added should be in a
    // call to closeClauseData(..)
    SLOW_DEBUG_DO(assert(cl != sentinel));
    vs_cls_data()[size++] = cl;
  }

  inline void close_cls_data() {
    vs_cls_data()[size++] = sentinel;
    assert(*(cls_begin()-1) == 0);
  }

  uint32_t* vars_begin() {
    return vs_cls_data();
  }

  uint32_t const* vars_begin() const {
    return vs_cls_data();
  }

  uint32_t* cls_begin() {
    return vs_cls_data() + clauses_offs;
  }


  uint32_t const* cls_begin() const {
    return vs_cls_data() + clauses_offs;
  }

  uint32_t nVars() const { return clauses_offs - 1; }
  uint32_t num_long_cls() const { return size - clauses_offs - 1; }
  uint32_t num_bin_cls() const { return bin_cls; }
  void set_num_bin_cls(const uint32_t n) { bin_cls = n; }
  uint32_t get_size() const { return size; }

  const uint32_t* get_raw_data() const { return vs_cls_data();}
  bool empty() const { return size == 0; }

  // Creates the full CNF as a component, at start-up. In other words, this is called ONCE
  void create_init_comp(uint32_t max_var_id, uint32_t max_cl_id, uint32_t _bin_cls) {
    clauses_offs = 1;
    bin_cls = _bin_cls;

    // Add all variables to top comp
    for (uint32_t v = 1; v <= max_var_id; v++) add_var(v);
    close_vars_data();

    //Add all clauses to top comp
    for (uint32_t clid = 1; clid <= max_cl_id; clid++) add_cl(clid);
    close_cls_data();
  }

  void clear() {
    size = 0;
    clauses_offs = 0;
    bin_cls = 0;
    id_ = 0;
  }

private:
  // data_ stores the comp data:
  // for better cache performance the
  // clause and variable data are stored in
  // a contiguous piece of memory
  // [variables SENTINEL clauses SENTINEL]
  // this order has to be taken care of on filling
  // in the data!
  //vector<uint32_t> vs_cls_data;
  uint32_t size = 0;

  uint32_t clauses_offs = 0;
  uint32_t bin_cls = 0;
  // id_ will identify denote the entry in the cacheable comp database,
  // where a Packed version of this comp is stored
  // yet this does not imply that the model count of this comp is already known
  // once the model count is known, a link to the packed comp will be stored
  // in the hash table
  CacheEntryID id_ = 0;
};

inline Comp* reserve_comp_space(uint32_t nVars, uint32_t num_clauses) {
  // vars, clauses, and the two sentinels
  uint32_t bytes_needed = (nVars + num_clauses + 2) * sizeof(uint32_t);
  bytes_needed += sizeof(Comp);
  Comp* ptr = (Comp*) malloc (bytes_needed);
  ptr->clear();
  return ptr;
}

}
