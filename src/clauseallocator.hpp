/******************************************
Copyright (C) 2009-2020 Mate Soos

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

#include <stdlib.h>
#include <cstdint>
#include <map>
#include <vector>

#include "counter_config.hpp"
#include "structures.hpp"

class Counter;

using std::vector;

class ClauseAllocator {
public:
  ClauseAllocator(const CounterConfiguration& _conf);
  ~ClauseAllocator();

  Clause* new_cl(bool _red, uint32_t sz) {
    void* mem = alloc_enough(sz);
    Clause* real = new (mem) Clause(_red, sz);
    return real;
  }

  ClauseOfs get_offset(const Clause* ptr) const;
  inline Clause* ptr(const ClauseOfs offset) const { return (Clause*)(&data_start[offset]); }
  void clause_free(Clause* c);
  void clause_free(ClauseOfs offset);
  bool consolidate(Counter* solver, const bool force = false);
  size_t mem_used() const;

private:
  void update_offsets(
    vector<ClauseOfs>& offsets,
    ClauseOfs* new_data_start,
    ClauseOfs*& new_ptr
  );

  void move_one_watchlist(vector<ClOffsBlckL>& ws, ClauseOfs* new_data_start, ClauseOfs*& new_ptr);

  ClauseOfs move_cl(
    ClauseOfs* new_data_start
    , ClauseOfs*& new_ptr
    , Clause* old
  ) const;

  uint32_t* data_start; ///<Stack starts at these positions
  uint64_t size; ///<The number of BASE_DATA_TYPE datapieces currently used in each stack
  uint64_t capacity; ///<The number of BASE_DATA_TYPE datapieces allocated
  uint64_t currently_used_sz; ///< The estimated used size of the stack
  const CounterConfiguration& conf;
  void* alloc_enough(const uint32_t num_lits);
};
