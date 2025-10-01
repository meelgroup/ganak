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
#include "comp_types/cacheable_comp.hpp"
#include "statistics.hpp"
#include <set>
#include <gmpxx.h>
#include "stack.hpp"

using std::set;

namespace GanakInt {

// There is EXACTLY ONE of this
class CompCacheIF {
public:
  CompCacheIF();
  virtual ~CompCacheIF();

  virtual CacheEntryID add_new_comp(void* comp, CacheEntryID super_comp_id) = 0;
  virtual uint32_t get_extra_bytes(void* comp) const = 0;
  virtual bool find_comp_and_incorporate_cnt(StackLevel &top, const uint32_t nvars, const void* comp) = 0;
  virtual void* create_new_comp(const Comp &comp, uint64_t hash_seed, const BPCSizes& bpc) = 0;
  virtual void free_comp(void* comp) = 0;

  virtual void make_entry_deletable(CacheEntryID id) = 0;
  virtual void init(Comp &super_comp, uint64_t hash_seed, const BPCSizes& bpc) = 0;
  virtual uint64_t get_num_entries_used() const = 0;
  virtual uint64_t get_max_num_entries() const = 0;
  virtual uint64_t compute_size_allocated() = 0;
  virtual bool cache_full(uint64_t extra_will_be_added = 0) const = 0;
  virtual bool exists(CacheEntryID id) const = 0;
  virtual void unlink(CacheEntryID id) = 0;
  virtual uint64_t clean_pollutions_involving(CacheEntryID id) = 0;
  virtual void erase(CacheEntryID id) = 0;
  virtual void store_value(const CacheEntryID id, const FF& model_count) = 0;
  virtual uint64_t calc_cutoff() const = 0;
  virtual bool delete_some_entries() = 0;
  virtual uint64_t unlink_from_tree(CacheEntryID id) = 0;
  virtual uint64_t num_descendants(CacheEntryID id) = 0;
  virtual uint64_t num_siblings(CacheEntryID id) = 0;
  virtual void test_descendantstree_consistency() = 0;
  virtual void debug_mem_data() const = 0;
};

}
