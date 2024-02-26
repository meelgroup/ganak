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
class CompCache {
public:
  CompCache(DataAndStatistics &_stats, const CounterConfiguration& _conf);
  ~CompCache() { for(auto& c: entry_base) c.set_free(); }

  void init(Comp &super_comp, void* hash_seed);
  void delete_comps_with_vars(const set<uint32_t>& vars);
  uint64_t get_num_entries_used() const
  {
    uint64_t ret = 0;
    for (uint32_t id = 2; id < entry_base.size(); id++)
      if (!entry_base[id].is_free()) ret++;

    return ret;
  }

  // compute the size in bytes of the comp cache from scratch
  // the value is stored in bytes_memory_usage_
  uint64_t compute_size_allocated();
  bool cache_full() const {
    return stats.cache_full(free_entry_base_slots.size() * sizeof(CacheableComp));
  }

  CacheableComp &entry(CacheEntryID id) { return entry_base[id]; }
  const CacheableComp &entry(CacheEntryID id) const { return entry_base[id]; }
  CacheableComp &entry(const Comp& comp) { return entry(comp.id()); }
  const CacheableComp &entry(const Comp& comp) const { return entry(comp.id()); }
  bool exists(CacheEntryID id) const { return !entry_base[id].is_free(); }

  // removes the entry id from the hash table
  // but not from the entry base
  inline void unlink(CacheEntryID id);

  // we delete the Comp with ID id
  // and all its descendants from the cache
  inline uint64_t clean_pollutions_involving(CacheEntryID id);

  // creates a CCacheEntry in the entry base
  // which contains a packed copy of comp
  // returns the id of the entry created
  // stores in the entry the position of
  // comp which is a part of the comp stack
  inline CacheEntryID new_comp(CacheableComp &ccomp, CacheEntryID super_comp_id);

  bool find_comp_and_incorporate_cnt(StackLevel &top, const uint32_t nvars, const CacheableComp &packed_comp) {
    stats.num_cache_look_ups_++;
    uint32_t table_ofs = packed_comp.get_hashkey() & tbl_size_mask;
    CacheEntryID act_id = table[table_ofs];
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

  // unchecked erase of an entry from entry_base
  void erase(CacheEntryID id) {
    stats.incorporate_cache_erase(entry_base[id]);
    entry_base[id].set_free();
    free_entry_base_slots.push_back(id);
  }


  // store the number in model_count as the model count of CacheEntryID id
  inline void store_value(const CacheEntryID id, const mpz_class &model_count);

  bool delete_some_entries();

  // delete entries, keeping the descendants tree consistent
  inline void unlink_from_tree(CacheEntryID id);

  // test function to ensure consistency of the descendant tree
  inline void test_descendantstree_consistency();
  void debug_dump_data();
private:

  void consider_cache_resize(){
    if (entry_base.size() > table.size()) rehash_table(2*table.size());
  }

  void rehash_table(const uint32_t size) {
    table.clear();
    table.resize(size,0);
    // we assert that table size is a power of 2
    // otherwise the tbl_size_mask doesn't work
    assert((table.size() & (table.size() - 1)) == 0);
    tbl_size_mask = table.size() - 1;
    for (uint32_t id = 2; id < entry_base.size(); id++)
      if (!entry_base[id].is_free()) {
        entry_base[id].set_next_bucket_element(0);
        if(entry_base[id].model_count_found()) {
          uint32_t table_ofs=table_entry(id);
          entry_base[id].set_next_bucket_element(table[table_ofs]);
          table[table_ofs] = id;
        }
    }
  }

  uint32_t table_entry(CacheEntryID id) const {
    return entry(id).get_hashkey() & tbl_size_mask;
  }

  void add_descendant(CacheEntryID compid, CacheEntryID descendantid) {
      assert(descendantid != entry(compid).first_descendant());
      entry(descendantid).set_next_sibling(entry(compid).first_descendant());
      entry(compid).set_first_descendant(descendantid);
    }

  vec<CacheableComp> entry_base;
  vector<CacheEntryID> free_entry_base_slots;

  // the actual hash table
  // by means of which the cache is accessed
  vector<CacheEntryID> table;

  uint32_t tbl_size_mask; // table is always power-of-two size

  DataAndStatistics &stats;
  const CounterConfiguration &conf;
  uint64_t my_time = 0;
};

CacheEntryID CompCache::new_comp(CacheableComp &ccomp, CacheEntryID super_comp_id){
  CacheEntryID id;

  while (cache_full()) {
    verb_print(1, "Cache full. Deleting some entries.");
    delete_some_entries();
  }

  assert(!cache_full());
  ccomp.set_last_used_time(my_time++);

  if (free_entry_base_slots.empty()) {
    /* bool at_capacity = (entry_base.capacity() == entry_base.size()); */
    /* if (at_capacity) entry_base_reserve(entry_base.capacity()*1.3); */
    entry_base.push_back(ccomp);
    id = entry_base.size() - 1;
  } else {
    id = free_entry_base_slots.back();
    assert(id < entry_base.size());
    assert(entry_base[id].is_free());
    free_entry_base_slots.pop_back();
    entry_base[id] = ccomp;
  }
  compute_size_allocated();
  VERBOSE_DEBUG_DO(if (stats.total_num_cached_comps_ % 100000 == 99999) debug_dump_data());

  entry(id).set_father(super_comp_id);
  add_descendant(super_comp_id, id);
  SLOW_DEBUG_DO(assert(exists(id)));
  SLOW_DEBUG_DO(assert(exists(super_comp_id)));


#ifdef SLOW_DEBUG
  for (uint32_t u = 2; u < entry_base.size(); u++)
    if (!entry_base[u].is_free()) {
      assert(entry_base[u].father() != id);
      /* assert(entry_base[u].first_descendant() != id); */
      assert(entry_base[u].next_sibling() != id);
    }
#endif
  return id;
}

// Recursively unlinks & removes id and its descendants
uint64_t CompCache::clean_pollutions_involving(const CacheEntryID id) {
  uint64_t removed = 0;

  // unlink id from the father's siblings list
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

  // Recursively unlink & delete all children
  CacheEntryID next_child = entry(id).first_descendant();
  entry(id).set_first_descendant(0);
  while (next_child) {
    const CacheEntryID act_child = next_child;
    next_child = entry(act_child).next_sibling();
    removed += clean_pollutions_involving(act_child);
  }

  // Finally, remove & erase the ID
  unlink(id);
  erase(id);
  return 1+removed;
}

void CompCache::unlink(CacheEntryID id) {
  uint32_t act_id = table[table_entry(id)];
  if (act_id == id){
    table[table_entry(id)] = entry(act_id).next_bucket_element();
  } else {
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

// Used only during cache freeing. Unlinks from descendants tree
void CompCache::unlink_from_tree(CacheEntryID id) {
  assert(exists(id));
  // we need a father for this all to work
  assert(entry(id).father());
  assert(exists(entry(id).father()));


  // unlink id from the father's siblings list
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

  // link the children of this one as siblings to the current siblings
  CacheEntryID act_child = entry(id).first_descendant();
  while (act_child) {
    CacheEntryID next_child = entry(act_child).next_sibling();
    entry(act_child).set_father(father);
    entry(act_child).set_next_sibling(entry(father).first_descendant());
    entry(father).set_first_descendant(act_child);
    act_child = next_child;
  }
}

void CompCache::store_value(const CacheEntryID id, const mpz_class &model_count) {
#ifdef CHECK_COUNT
  //we disable cache on check_count, to remove an error source
  return;
#endif

  consider_cache_resize();
  uint32_t table_ofs = table_entry(id);
  // when storing the new model count the size of the model count
  // and hence that of the comp will change
  SLOW_DEBUG_DO(assert(!entry(id).is_free()));
  /* SLOW_DEBUG_DO(assert(stats.sum_bytes_cached_comps_ > entry(id).size_in_bytes())); */
  stats.sum_bytes_cached_comps_ -= entry(id).size_in_bytes();
  entry(id).set_model_count(model_count);
  entry(id).set_last_used_time(my_time);
  entry(id).set_next_bucket_element(table[table_ofs]);
  table[table_ofs] = id;
  stats.sum_bytes_cached_comps_ += entry(id).size_in_bytes();
}
