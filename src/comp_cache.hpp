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
#include "counter_config.hpp"
#include "Vec.hpp"
#include <set>
#include <gmpxx.h>
#include "time_mem.hpp"

#include "comp_types/comp.hpp"
#include "stack.hpp"
using std::set;

namespace GanakInt {

// There is EXACTLY ONE of this
class CompCache {
public:
  CompCache(const uint32_t num_vars, DataAndStatistics &_stats, const CounterConfiguration& _conf);
  ~CompCache() { for(auto& c: entry_base) c.set_free(); }

  void init(Comp &super_comp, void* hash_seed);
  uint64_t get_num_entries_used() const {
    uint64_t ret = 0;
    for (uint32_t id = 2; id < entry_base.size(); id++)
      if (!entry_base[id].is_free()) ret++;

    return ret;
  }

  // compute the size in bytes of the comp cache from scratch
  // the value is stored in bytes_memory_usage_
  uint64_t compute_size_allocated();
  bool cache_full(uint64_t extra_will_be_added = 0) const {
    return stats.cache_full(free_entry_base_slots.size() * sizeof(CacheableComp),
        extra_will_be_added);
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
  inline uint64_t calc_extra_mem_after_push() const;

  bool find_comp_and_incorporate_cnt(StackLevel &top, const uint32_t nvars, const CacheableComp &packed_comp) {
    stats.num_cache_look_ups++;
    uint32_t table_ofs = packed_comp.get_hashkey() & tbl_size_mask;
    CacheEntryID act_id = table[table_ofs];
    if (!act_id) return false;
    while(act_id){
      if (entry(act_id).equals_clhashkey(packed_comp)) {
        stats.incorporate_cache_hit(nvars);
        switch(conf.cache_time_update) {
          case 0: break;
          case 1: entry(act_id).set_last_used_time(my_time); break;
          case 2: entry(act_id).avg_last_used_time(my_time, 2); break;
          case 3: entry(act_id).avg_last_used_time(my_time, 3); break;
          default: release_assert(false);
        }
        debug_print(COLYEL2 << "Cache hit. cache ID: " << act_id);
        top.include_solution(entry(act_id).model_count());
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
  inline void store_value(const CacheEntryID id, const FF& model_count);

  double calc_cutoff() const;
  bool delete_some_entries();

  // delete entries, keeping the descendants tree consistent
  uint64_t unlink_from_tree(CacheEntryID id);
  uint64_t num_descendants(CacheEntryID id);
  uint64_t num_siblings(CacheEntryID id);

  // test function to ensure consistency of the descendant tree
  inline void test_descendantstree_consistency();
  void debug_mem_data() const;

private:
  uint64_t freeram();
  void consider_table_resize() {
    // NOTE: it's possible to e.g. half the table.size() here, but
    // the performance gain vs mem use is not worth it
    // Good example file to stress: mc2023_track1_138.cnf
    if (entry_base.size() > table.size()) {
      double vm_before;
      auto used_before = mem_used(vm_before);
      rehash_table(2*table.size());
      double vm_after;
      auto used_after = mem_used(vm_after);
      verb_print(2,
        "table resize -- used before: " << used_before/(double)(1024*1024)
        << " vm used before: " << vm_before/(double)(1024*1024)
        << " used after: " << used_after/(double)(1024*1024)
        << " vm used after: " << vm_after/(double)(1024*1024)
        << " total T: " << cpu_time());
    }
  }

  void rehash_table(const uint32_t size) {
    table.clear();
    table.resize(size);
    table.shrink_to_fit();
    std::fill(table.begin(), table.end(), 0);
    assert((table.size() & (table.size() - 1)) == 0 && "Table size must be a power of 2");
    tbl_size_mask = table.size() - 1;

    for (uint32_t id = 2; id < entry_base.size(); id++)
      if (!entry_base[id].is_free()) {
        entry_base[id].set_next_bucket_element(0);
        if(entry_base[id].model_count_found()) {
          uint32_t table_ofs=table_pos(id);
          entry_base[id].set_next_bucket_element(table[table_ofs]);
          table[table_ofs] = id;
        }
    }
  }

  uint32_t table_pos(CacheEntryID id) const {
    return entry(id).get_hashkey() & tbl_size_mask;
  }

  void add_descendant(CacheEntryID compid, CacheEntryID descendantid) {
      assert(descendantid != entry(compid).first_descendant());
      entry(descendantid).set_next_sibling(entry(compid).first_descendant());
      entry(compid).set_first_descendant(descendantid);
    }

  vec<CacheableComp> entry_base;
  vec<CacheEntryID> free_entry_base_slots;

  // the actual hash table
  // by means of which the cache is accessed
  vec<CacheEntryID> table;

  uint32_t tbl_size_mask; // table is always power-of-two size

  DataAndStatistics &stats;
  const CounterConfiguration &conf;
  uint64_t my_time = 0;
  uint32_t num_vars;
};

inline uint64_t CompCache::calc_extra_mem_after_push() const {
  bool at_capacity = entry_base.capacity() == entry_base.size();
  bool at_capacity_table = table.capacity() == table.size();
  uint64_t extra_will_be_added = 0;

  // assume it will be multiplied by 1.5
  if (at_capacity) extra_will_be_added =
    (entry_base.capacity()*sizeof(CacheableComp) + stats.sum_bignum_bytes)/2;
  if (at_capacity_table) extra_will_be_added += (table.capacity()*sizeof(CacheEntryID))/2;
  return extra_will_be_added;
}

inline CacheEntryID CompCache::new_comp(CacheableComp &ccomp, CacheEntryID super_comp_id) {
  uint64_t extra_mem_with_push = calc_extra_mem_after_push();
  while (cache_full(extra_mem_with_push)) {
    verb_print(1, "Cache full. Deleting some entries.");
    if (conf.verb >= 2) debug_mem_data();
    delete_some_entries();
    if (conf.verb >= 2) debug_mem_data();
    extra_mem_with_push = calc_extra_mem_after_push();
  }

  ccomp.set_last_used_time(my_time++);

  CacheEntryID id;
  if (free_entry_base_slots.empty()) {
    bool at_capacity = entry_base.capacity() == entry_base.size();
    if (at_capacity && conf.verb >= 3) {
      double vm_dat;
      auto dat = mem_used(vm_dat);
      verb_print(3,std::setw(40) << "After enlarge entry_base mem use MB: " <<
        (double)(entry_base.capacity()*sizeof(CacheableComp))/(double)(1024*1024));
      verb_print(3,
        "Before entry enlarge Total process MB : " << dat/(double)(1024*1024)
        << " Total process vm MB: " << vm_dat/(double)(1024*1024));

    }
    entry_base.push_back(ccomp);
    if (at_capacity && conf.verb >= 3) {
      double vm_dat;
      double dat = mem_used(vm_dat);
      verb_print(3,std::setw(40) << "After enlarge entry_base mem use MB: " <<
        (double)(entry_base.capacity()*sizeof(CacheableComp))/(double)(1024*1024));
      verb_print(3,
        "After entry enlarge Total process MB  : " << dat/(double)(1024*1024)
        << " Total process vm MB: " << vm_dat/(double)(1024*1024));
    }
    id = entry_base.size() - 1;
  } else {
    id = free_entry_base_slots.back();
    assert(id < entry_base.size());
    assert(entry_base[id].is_free());
    free_entry_base_slots.pop_back();
    entry_base[id] = ccomp;
  }
  compute_size_allocated(); // TODO expensive... and should not be needed
  VERBOSE_DEBUG_DO(if (stats.total_num_cached_comps % 100000 == 99999) debug_mem_data());

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
inline uint64_t CompCache::clean_pollutions_involving(const CacheEntryID id) {
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

inline void CompCache::unlink(CacheEntryID id) {
  uint32_t act_id = table[table_pos(id)];
  if (act_id == id){
    table[table_pos(id)] = entry(act_id).next_bucket_element();
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

inline void CompCache::store_value(const CacheEntryID id, const FF& model_count) {
  consider_table_resize();
  uint32_t table_ofs = table_pos(id);
  // when storing the new model count the size of the model count
  // and hence that of the comp will change
  SLOW_DEBUG_DO(assert(!entry(id).is_free()));
  stats.sum_bignum_bytes -= entry(id).bignum_bytes();
  entry(id).set_model_count(model_count);
  entry(id).set_last_used_time(my_time);
  entry(id).set_next_bucket_element(table[table_ofs]);
  table[table_ofs] = id;
  stats.sum_bignum_bytes += entry(id).bignum_bytes();
}

inline CompCache::CompCache(
    uint32_t _num_vars,
    DataAndStatistics &_stats, const CounterConfiguration &_conf):
  stats(_stats), conf(_conf), num_vars(_num_vars) {}

inline void CompCache::init(Comp &super_comp, void* hash_seed){
  CacheableComp *packed_super_comp;
  vector<uint32_t> tmp(100+super_comp.nVars()+super_comp.num_long_cls());
  packed_super_comp = new CacheableComp(hash_seed,super_comp);
  my_time = 1;

  entry_base.clear();
  auto x = CacheableComp();
  entry_base.push_back(x); // dummy Element
  stats.incorporate_cache_store(x, super_comp.nVars());
  table.clear();
  table.resize(1024*1024);
  std::fill(table.begin(), table.end(), 0);
  tbl_size_mask = table.size() - 1;

  free_entry_base_slots.clear();

  const uint64_t free_ram = freeram();
  uint64_t max_cache_bound = 80 * (free_ram / 100);

  if (stats.max_cache_size_bytes == 0) {
    stats.max_cache_size_bytes = max_cache_bound;
  }

  if (stats.max_cache_size_bytes > free_ram) {
    verb_print(1, "WARNING: Maximum cache size larger than free RAM available");
    verb_print(1, "Free RAM " << std::setprecision(2)
      << (double)free_ram / (1024.0*1024.0) << "MB");
  }
  verb_print(2, "Max cache size (80% free mem-200MB): "
    << stats.max_cache_size_bytes / (1024ULL*1024ULL) << " MB");


  stats.sum_bignum_bytes = 0;
  stats.cache_infra_bytes_mem_usage = 0;
  assert(!cache_full());
  entry_base.push_back(*packed_super_comp);
  stats.incorporate_cache_store(*packed_super_comp, super_comp.nVars());
  delete packed_super_comp;
  super_comp.set_id(1);
  compute_size_allocated();
}

}
