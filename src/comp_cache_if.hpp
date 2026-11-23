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
#include <gmpxx.h>
#include "stack.hpp"
#include <functional>

namespace GanakInt {

// There is EXACTLY ONE of this
class CompCacheIF {
public:
  CompCacheIF();
  virtual ~CompCacheIF();

  virtual void init(Comp &super_comp, uint32_t hash_seed, const BPCSizes& bpc,
      const vector<Comp*>* comp_stack) = 0;
  // Builds the lookup key of comp, used by the next find_probe()/add_probe()
  virtual void set_probe(const Comp &comp, uint32_t hash_seed, const BPCSizes& bpc) = 0;
  virtual bool find_probe_and_incorporate_cnt(StackLevel &top, uint32_t nvars,
      int* out_node = nullptr) = 0;
  // May evict, which moves entries: the ids of the comps on comp_stack are updated
  virtual CacheEntryID add_probe(uint32_t nvars, const Comp& super_comp) = 0;
  // d-DNNF: associate a circuit node id with a cache entry (retrieved on a hit via
  // out_node above). No-op unless compiling.
  virtual void set_compile_node(CacheEntryID, int) {}
  virtual void make_entry_deletable(CacheEntryID id) = 0;
  // Ids held outside the cache, remapped by eviction like the comp stack's
  void set_anchor_source(std::function<void(vector<CacheEntryID*>&)> f) { anchor_source = std::move(f); }
  std::function<void(vector<CacheEntryID*>&)> anchor_source;
  [[nodiscard]] virtual uint64_t get_num_entries_used() const = 0;
  [[nodiscard]] virtual uint64_t get_max_num_entries() const = 0;
  [[nodiscard]] virtual bool exists(CacheEntryID id) const = 0;
  // Removes every entry with id >= from. Entries are kept in creation order, so
  // these are the comps created since `from` was, and all their descendants.
  virtual uint64_t remove_from(CacheEntryID from) = 0;
  virtual void store_value(CacheEntryID id, const FF& model_count) = 0;
  virtual void debug_mem_data() const = 0;
#ifdef SLOW_DEBUG
  virtual void check_removal(CacheEntryID from, const vector<CacheEntryID>& roots, CacheEntryID kept_father) const = 0;
#endif
};

}
