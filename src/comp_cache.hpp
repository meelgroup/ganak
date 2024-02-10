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

#include "comp_types/base_packed_comp.hpp"
#include "statistics.hpp"
#include "counter_config.hpp"
#include "Vec.hpp"
#include <set>
#include <gmpxx.h>

#include "comp_types/comp.hpp"
#include "stack.hpp"
using std::set;

// There is EXACTLY ONE of this
class ComponentCache {
public:
  ComponentCache(DataAndStatistics &_stats, const CounterConfiguration& _conf);
  ~ComponentCache() { for(auto& c: entry_base_) c.set_free(); }

  void init(Component &super_comp, void* hash_seed);
  void delete_comps_with_vars(const set<uint32_t>& vars);
  uint64_t get_num_entries_used() const
  {
    uint64_t ret = 0;
    for (uint32_t id = 2; id < entry_base_.size(); id++)
      if (!entry_base_[id].is_free()) ret++;

    return ret;
  }

  // compute the size in bytes of the comp cache from scratch
  // the value is stored in bytes_memory_usage_
  uint64_t compute_size_allocated();
  bool cache_full() const {
    return stats.cache_full(free_entry_base_slots_.size() * sizeof(CacheableComponent));
  }

  CacheableComponent &entry(CacheEntryID id) { return entry_base_[id]; }
  const CacheableComponent &entry(CacheEntryID id) const { return entry_base_[id]; }
  CacheableComponent &entry(const Component& comp) { return entry(comp.id()); }
  const CacheableComponent &entry(const Component& comp) const { return entry(comp.id()); }
  bool hasEntry(CacheEntryID id) const { return !entry_base_[id].is_free(); }

  // removes the entry id from the hash table
  // but not from the entry base
  inline void removeFromHashTable(CacheEntryID id);

  // we delete the Component with ID id
  // and all its descendants from the cache
  inline uint64_t cleanPollutionsInvolving(CacheEntryID id);

  // creates a CCacheEntry in the entry base
  // which contains a packed copy of comp
  // returns the id of the entry created
  // stores in the entry the position of
  // comp which is a part of the comp stack
  inline CacheEntryID storeAsEntry(CacheableComponent &ccomp,
                            CacheEntryID super_comp_id);

  bool manageNewComponent(StackLevel &top,
      const uint32_t nvars, const CacheableComponent &packed_comp) {
    stats.num_cache_look_ups_++;
    uint32_t table_ofs = packed_comp.get_hashkey() & table_size_mask_;
    CacheEntryID act_id = table_[table_ofs];
    if (!act_id) return false;
    while(act_id){
      if (entry(act_id).equals_clhashkey(packed_comp)) {
        stats.incorporate_cache_hit(nvars);
        top.includeSolution(entry(act_id).model_count());
        return true;
      }
      act_id = entry(act_id).next_bucket_element();
    }
    return false;
  }

  // unchecked erase of an entry from entry_base_
  void eraseEntry(CacheEntryID id) {
    stats.incorporate_cache_erase(entry_base_[id]);
    entry_base_[id].set_free();
    free_entry_base_slots_.push_back(id);
  }


  // store the number in model_count as the model count of CacheEntryID id
  inline void storeValueOf(CacheEntryID id, const mpz_class &model_count);

  bool deleteEntries();

  // delete entries, keeping the descendants tree consistent
  inline void removeFromDescendantsTree(CacheEntryID id);

  // test function to ensure consistency of the descendant tree
  inline void test_descendantstree_consistency();
  void debug_dump_data();
private:

  void considerCacheResize(){
    if (entry_base_.size() > table_.size()) {
      reHashTable(2*table_.size());
    }
  }
  void reHashTable(const uint32_t size) {
    table_.clear();
    table_.resize(size,0);
    // we assert that table size is a power of 2
    // otherwise the table_size_mask_ doesn't work
    assert((table_.size() & (table_.size() - 1)) == 0);
    table_size_mask_ = table_.size() - 1;
    for (uint32_t id = 2; id < entry_base_.size(); id++)
      if (!entry_base_[id].is_free()) {
        entry_base_[id].set_next_bucket_element(0);
        if(entry_base_[id].modelCountFound()) {
          uint32_t table_ofs=tableEntry(id);
          entry_base_[id].set_next_bucket_element(table_[table_ofs]);
          table_[table_ofs] = id;
        }
    }
  }

  uint32_t tableEntry(CacheEntryID id) const {
    return entry(id).get_hashkey() & table_size_mask_;
  }
  void add_descendant(CacheEntryID compid, CacheEntryID descendantid) {
      assert(descendantid != entry(compid).first_descendant());
      entry(descendantid).set_next_sibling(entry(compid).first_descendant());
      entry(compid).set_first_descendant(descendantid);
    }

  void remove_firstdescendantOf(CacheEntryID compid) {
      CacheEntryID desc = entry(compid).first_descendant();
      if (desc != 0)
        entry(compid).set_first_descendant(entry(desc).next_sibling());
    }

  vec<CacheableComponent> entry_base_;
  vector<CacheEntryID> free_entry_base_slots_;

  // the actual hash table
  // by means of which the cache is accessed
  vector<CacheEntryID> table_;

  uint32_t table_size_mask_; // table is always power-of-two size

  DataAndStatistics &stats;
  const CounterConfiguration &conf;
  uint64_t my_time_ = 0;
};

CacheEntryID ComponentCache::storeAsEntry(CacheableComponent &ccomp, CacheEntryID super_comp_id){
  CacheEntryID id;

  while (cache_full()) {
    verb_print(1, "Cache full. Deleting some entries.");
    deleteEntries();
  }

  assert(!cache_full());
  ccomp.set_creation_time(my_time_++);

  if (free_entry_base_slots_.empty()) {
    /* bool at_capacity = (entry_base_.capacity() == entry_base_.size()); */
    /* if (at_capacity) entry_base_.reserve(entry_base_.capacity()*1.3); */
    entry_base_.push_back(ccomp);
    id = entry_base_.size() - 1;
  } else {
    id = free_entry_base_slots_.back();
    assert(id < entry_base_.size());
    assert(entry_base_[id].is_free());
    free_entry_base_slots_.pop_back();
    entry_base_[id] = ccomp;
  }
  compute_size_allocated();
  VERBOSE_DEBUG_DO(if (stats.total_num_cached_comps_ % 100000 == 99999) debug_dump_data());

  entry(id).set_father(super_comp_id);
  add_descendant(super_comp_id, id);
  SLOW_DEBUG_DO(assert(hasEntry(id)));
  SLOW_DEBUG_DO(assert(hasEntry(super_comp_id)));

  stats.incorporate_cache_store(entry(id), 0);

#ifdef SLOW_DEBUG
  for (uint32_t u = 2; u < entry_base_.size(); u++)
    if (!entry_base_[u].is_free()) {
      assert(entry_base_[u].father() != id);
      /* assert(entry_base_[u].first_descendant() != id); */
      assert(entry_base_[u].next_sibling() != id);
    }
#endif
  return id;
}

uint64_t ComponentCache::cleanPollutionsInvolving(const CacheEntryID id) {
  uint64_t removed = 0;
  const CacheEntryID father = entry(id).father();
  if (entry(father).first_descendant() == id) {
    entry(father).set_first_descendant(entry(id).next_sibling());
  } else {
    CacheEntryID act_sibl = entry(father).first_descendant();
    while (act_sibl) {
      CacheEntryID next_sibl = entry(act_sibl).next_sibling();
      if (next_sibl == id) {
        entry(act_sibl).set_next_sibling(entry(next_sibl).next_sibling());
        break;
      }
      act_sibl = next_sibl;
    }
  }
  CacheEntryID next_child = entry(id).first_descendant();
  entry(id).set_first_descendant(0);
  while (next_child) {
    CacheEntryID act_child = next_child;
    next_child = entry(act_child).next_sibling();
    removed+=cleanPollutionsInvolving(act_child);
  }
  removeFromHashTable(id);
  eraseEntry(id);
  return 1+removed;
}

void ComponentCache::removeFromHashTable(CacheEntryID id) {
  //assert(false);
  uint32_t act_id = table_[tableEntry(id)];
  if(act_id == id){
    table_[tableEntry(id)] = entry(act_id).next_bucket_element();
  }
  else {
  while (act_id) {
        CacheEntryID next_id = entry(act_id).next_bucket_element();
        if (next_id == id) {
          entry(act_id).set_next_bucket_element(entry(next_id).next_bucket_element());
          break;
        }
        act_id = next_id;
      }
  }
}

void ComponentCache::removeFromDescendantsTree(CacheEntryID id) {
  assert(hasEntry(id));
  // we need a father for this all to work
  assert(entry(id).father());
  assert(hasEntry(entry(id).father()));
  // two steps
  // 1. remove id from the siblings list
  CacheEntryID father = entry(id).father();
  if (entry(father).first_descendant() == id) {
    entry(father).set_first_descendant(entry(id).next_sibling());
  } else {
    CacheEntryID act_sibl = entry(father).first_descendant();
    while (act_sibl) {
      CacheEntryID next_sibl = entry(act_sibl).next_sibling();
      if (next_sibl == id) {
        entry(act_sibl).set_next_sibling(entry(next_sibl).next_sibling());
        break;
      }
      act_sibl = next_sibl;
    }
  }

  // 2. add the children of this one as
  //    siblings to the current siblings
  CacheEntryID act_child = entry(id).first_descendant();
  while (act_child) {
    CacheEntryID next_child = entry(act_child).next_sibling();
    entry(act_child).set_father(father);
    entry(act_child).set_next_sibling(entry(father).first_descendant());
    entry(father).set_first_descendant(act_child);
    act_child = next_child;
  }
}

void ComponentCache::storeValueOf(CacheEntryID id, const mpz_class &model_count) {
#ifdef CHECK_COUNT
  //we disable cache on check_count, to remove an error source
  return;
#endif
  considerCacheResize();
  uint32_t table_ofs = tableEntry(id);
  // when storing the new model count the size of the model count
  // and hence that of the comp will change
  SLOW_DEBUG_DO(assert(!entry(id).is_free()));
  /* SLOW_DEBUG_DO(assert(stats.sum_bytes_cached_comps_ > entry(id).SizeInBytes())); */
  stats.sum_bytes_cached_comps_ -= entry(id).SizeInBytes();
  entry(id).set_model_count(model_count,my_time_);
  entry(id).set_creation_time(my_time_);
  entry(id).set_next_bucket_element(table_[table_ofs]);
  table_[table_ofs] = id;
  stats.sum_bytes_cached_comps_ += entry(id).SizeInBytes();
}
