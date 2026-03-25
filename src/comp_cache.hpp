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
#include <algorithm>
#include <numeric>
#include <gmpxx.h>
#include "structures.hpp"
#include "time_mem.hpp"
#include "comp_cache_if.hpp"

#include "comp_types/comp.hpp"
#include "stack.hpp"
namespace GanakInt {

// There is EXACTLY ONE of this
template<typename T>
class CompCache final: public CompCacheIF {
public:
  CompCache(DataAndStatistics &_stats, const CounterConfiguration &_conf) : stats(_stats), conf(_conf) { }
  ~CompCache() override = default;

  void init(Comp &super_comp, uint64_t hash_seed, const BPCSizes& bpc) override;
  [[nodiscard]] uint64_t get_num_entries_used() const override {
    // Skip slots 0 and 1 which are sentinel/reserved entries
    return static_cast<uint64_t>(std::count_if(
      entry_base.begin() + 2, entry_base.end(),
      [](const auto& e) { return !e.is_free(); }));
  }
  uint64_t get_extra_bytes(void* c) const override {
    T* comp = reinterpret_cast<T*>(c);
    return comp->extra_bytes();
  }
  void* create_new_comp(const Comp &comp, uint64_t hash_seed, const BPCSizes& bpc) override {
    return new T(comp, hash_seed, bpc);
  }

  [[nodiscard]] uint64_t get_max_num_entries() const override { return entry_base.size(); }

  // compute the size in bytes of the comp cache from scratch
  // the value is stored in bytes_memory_usage_
  uint64_t compute_size_allocated() override;
  [[nodiscard]] bool cache_full(uint64_t extra_will_be_added = 0) const override {
    // Each slot in free_entry_base_slots represents a T-sized slot in entry_base that
    // is already allocated (paid for in cache_infra_bytes_mem_usage) but currently unused.
    // Subtracting this avoids counting free entry_base slots as live memory.
    return stats.cache_full(free_entry_base_slots.size() * sizeof(T),
        extra_will_be_added);
  }

  void make_entry_deletable(CacheEntryID id) override {
    entry(id).set_deletable();
  }

  [[nodiscard]] bool exists(CacheEntryID id) const override {
    return !entry_base[id].is_free();
  }

  // removes the entry id from the hash table
  // but not from the entry base
  void unlink(CacheEntryID id) override;

  // we delete the Comp with ID id
  // and all its descendants from the cache
  uint64_t clean_pollutions_involving(CacheEntryID id) override;

  // creates a CCacheEntry in the entry base
  // which contains a packed copy of comp
  // returns the id of the entry created
  // stores in the entry the position of
  // comp which is a part of the comp stack
  CacheEntryID add_new_comp(void* comp, CacheEntryID super_comp_id) override;
  void update_entry_time(CacheEntryID id) {
    switch(conf.cache_time_update) {
      case 0: break;
      case 1: entry(id).set_last_used_time(my_time); break;
      case 2: entry(id).avg_last_used_time(my_time, 2); break;
      case 3: entry(id).avg_last_used_time(my_time, 3); break;
      default: release_assert(false);
    }
  }

  bool find_comp_and_incorporate_cnt(StackLevel &top, const uint32_t nvars, const void* c) override {
    const T& comp = *reinterpret_cast<const T*>(c);
    stats.num_cache_look_ups++;
    uint32_t table_ofs = (uint32_t)comp.get_hashkey() & tbl_size_mask;
    CacheEntryID act_id = table[table_ofs];
    if (!act_id) return false;
    while(act_id){ // LINEAR PROBING through collision chain
      if (entry(act_id).equals(comp)) {
        stats.incorporate_cache_hit(nvars);
        update_entry_time(act_id);
        debug_print(COLYEL2 << "Cache hit. cache ID: " << act_id);
        top.include_solution(entry(act_id).model_count());
        return true;
      }
      act_id = entry(act_id).next_bucket_element();
    }
    return false;
  }

  // unchecked erase of an entry from entry_base
  void erase(CacheEntryID id) override {
    stats.incorporate_cache_erase(entry_base[id].extra_bytes());
    entry_base[id].set_free();
    free_entry_base_slots.push_back(id);
    SLOW_DEBUG_DO({
      // Verify sum_extra_bytes matches a from-scratch computation across all live entries.
      // Catches any mismatch between the incremental stat and the true allocation.
      uint64_t actual_sum = std::accumulate(entry_base.begin(), entry_base.end(), 0ULL,
          [](uint64_t s, const auto& e) { return e.is_free() ? s : s + e.extra_bytes(); });
      assert(actual_sum == stats.sum_extra_bytes);
    });
  }

  // store the number in model_count as the model count of CacheEntryID id
  void store_value(const CacheEntryID id, const FF& model_count) override;

  [[nodiscard]] uint64_t calc_cutoff() const override;
  bool delete_some_entries() override;

  // delete entries, keeping the descendants tree consistent
  uint64_t unlink_from_tree(CacheEntryID id) override;
  uint64_t num_direct_children(CacheEntryID id) override;
  uint64_t num_siblings(CacheEntryID id) override;

  // test function to ensure consistency of the descendant tree
  void test_descendantstree_consistency() override;
  void debug_mem_data() const override;
  void free_comp(void* c) override {
    T* comp = reinterpret_cast<T*>(c);
    assert(comp != nullptr);
    delete comp;
  }

private:
  [[nodiscard]] uint64_t calc_extra_mem_after_push() const;
  T& entry(CacheEntryID id) { return entry_base[id]; }
  const T &entry(CacheEntryID id) const { return entry_base[id]; }
  T& entry(const Comp& comp) { return entry(comp.id()); }
  const T &entry(const Comp& comp) const { return entry(comp.id()); }
  void consider_table_resize() {
    // NOTE: only entries with model_count_found() are inserted in the hash table
    // so entry_base.size() is a conservative upper bound on
    // actual load factor. The threshold below (resize when entries exceed table size)
    // is intentionally conservative; the real load factor is always ≤ 1.0.
    if (entry_base.size() > table.size()) {
      double vm_before;
      auto used_before = mem_used(vm_before);
      rehash_table(2*table.size());
      double vm_after;
      auto used_after = mem_used(vm_after);
      verb_print(2,
        "table resize -- used before: " << in_mb(used_before)
        << " vm used before: " << in_mb(vm_before)
        << " used after: " << in_mb(used_after)
        << " vm used after: " << in_mb(vm_after)
        << " total T: " << cpu_time());
    }
  }

  void rehash_table(const uint32_t size) {
    table.clear();
    table.resize(size); // vec::growTo zero-initializes new elements
    table.shrink_to_fit();
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

  [[nodiscard]] uint32_t table_pos(CacheEntryID id) const {
    return (uint32_t)entry(id).get_hashkey() & tbl_size_mask;
  }

  void add_descendant(CacheEntryID compid, CacheEntryID descendantid) {
    assert(descendantid != entry(compid).first_descendant());
    entry(descendantid).set_next_sibling(entry(compid).first_descendant());
    entry(compid).set_first_descendant(descendantid);
  }

  vec<T> entry_base;
  vec<CacheEntryID> free_entry_base_slots;

  // the actual hash table
  // by means of which the cache is accessed
  vec<CacheEntryID> table;

  uint32_t tbl_size_mask = 0; // table is always power-of-two size; 0 until init()

  DataAndStatistics &stats;
  const CounterConfiguration &conf;
  uint64_t my_time = 0;
};

template<typename T>
uint64_t CompCache<T>::calc_extra_mem_after_push() const {
  bool at_capacity = entry_base.capacity() == entry_base.size();
  bool at_capacity_table = table.capacity() == table.size();
  uint64_t extra_will_be_added = 0;

  // entry_base backing store grows by ~50% on realloc; sum_extra_bytes is per-entry
  // heap data whose pointers are merely copied (not reallocated) on entry_base resize,
  // so it must not be included in the growth estimate
  if (at_capacity) extra_will_be_added = entry_base.capacity() * sizeof(T) / 2;
  if (at_capacity_table) extra_will_be_added += (table.capacity()*sizeof(CacheEntryID))/2;
  return extra_will_be_added;
}

template<typename T>
CacheEntryID CompCache<T>::add_new_comp(void* c, CacheEntryID super_comp_id) {
  T& comp = *reinterpret_cast<T*>(c);
  uint64_t extra_mem_with_push = calc_extra_mem_after_push();
  while (cache_full(extra_mem_with_push)) {
    verb_print(1, "Cache full. Deleting some entries.");
    debug_mem_data();
    const bool freed_any = delete_some_entries();
    debug_mem_data();
    extra_mem_with_push = calc_extra_mem_after_push();
    // if nothing was freed (all deletable entries had too many children/siblings),
    // break to avoid an infinite loop; we accept a slightly over-budget allocation
    if (!freed_any) break;
  }

  comp.set_last_used_time(my_time++);

  CacheEntryID id;
  if (free_entry_base_slots.empty()) {
    bool at_capacity = entry_base.capacity() == entry_base.size();
    if (at_capacity && conf.verb >= 3) {
      double vm_dat;
      const auto dat = mem_used(vm_dat);
      verb_print(3, "entry_base cap MB: " << in_mb(entry_base.capacity()*sizeof(comp)));
      verb_print(3,
        "Before entry enlarge Total process MB: " << in_mb(dat)
        << " Total process vm MB: " << in_mb(vm_dat));

    }
    entry_base.emplace_back(std::move(comp));
    if (at_capacity && conf.verb >= 3) {
      double vm_dat;
      const auto dat = mem_used(vm_dat);
      verb_print(3, "entry_base cap MB: " << in_mb(entry_base.capacity()*sizeof(comp)));
      verb_print(3,
        "After entry enlarge Total process MB: " << in_mb(dat)
        << " Total process vm MB: " << in_mb(vm_dat));
    }
    id = entry_base.size() - 1;
    // only recompute infra size when entry_base actually grew (potentially reallocated)
    if (at_capacity) compute_size_allocated();
  } else {
    id = free_entry_base_slots.back();
    assert(id < entry_base.size());
    assert(entry_base[id].is_free());
    free_entry_base_slots.pop_back();
    std::swap(entry_base[id], comp);
  }
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

// Iteratively unlinks & removes id and all its descendants.
// Previously recursive — can stack-overflow on deeply nested descendant trees.
template<typename T>
uint64_t CompCache<T>::clean_pollutions_involving(const CacheEntryID id) {
  // First, detach id from its parent's sibling list
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

  // Collect id and all descendants using an explicit work stack, then erase them all.
  // Using two passes avoids modifying the tree while traversing it.
  vector<CacheEntryID> to_erase;
  vector<CacheEntryID> work;
  work.push_back(id);
  while (!work.empty()) {
    const CacheEntryID cur = work.back();
    work.pop_back();
    CacheEntryID child = entry(cur).first_descendant();
    while (child) {
      work.push_back(child);
      child = entry(child).next_sibling();
    }
    to_erase.push_back(cur);
  }

  VERBOSE_DEBUG_DO(cout << "Pollution-based unlinking from tree. ID: " << id
                        << " total removed (incl. descendants): " << to_erase.size() << endl;);
  for (const CacheEntryID cur : to_erase) {
    unlink(cur);
    erase(cur);
  }
  return to_erase.size();
}

template<typename T>
void CompCache<T>::unlink(CacheEntryID id) {
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

template<typename T>
void CompCache<T>::store_value(const CacheEntryID id, const FF& model_count) {
  consider_table_resize();
  uint32_t table_ofs = table_pos(id);
  // when storing the new model count the size of the model count
  // and hence that of the comp will change
  SLOW_DEBUG_DO(assert(!entry(id).is_free()));
  stats.sum_extra_bytes -= entry(id).extra_bytes();
  entry(id).set_model_count(model_count);
  entry(id).set_last_used_time(my_time);
  entry(id).set_next_bucket_element(table[table_ofs]);
  table[table_ofs] = id;
  stats.sum_extra_bytes += entry(id).extra_bytes();
}

template<typename T>
void CompCache<T>::init(Comp &super_comp, uint64_t hash_seed, const BPCSizes& bpc) {
  // Release mem
  for (auto& e : entry_base) e.set_free();
  entry_base.clear();

  // Setup initial state
  my_time = 1;
  free_entry_base_slots.clear();
  stats.sum_extra_bytes = 0;
  stats.cache_infra_bytes_mem_usage = 0;

  // Size the initial table to fit within ~25% of the cache budget
  // Cap at 1M entries; floor at 256 entries.
  const uint32_t max_tbl = (uint32_t)std::min(
      (uint64_t)1024*1024,
      stats.max_cache_size_bytes / sizeof(CacheEntryID) / 4);
  uint32_t init_tbl = 256;
  while (init_tbl * 2 <= max_tbl) init_tbl *= 2;
  table.clear();
  table.resize(init_tbl);
  tbl_size_mask = table.size() - 1;
  verb_print(2, "Max cache size set: " << std::setprecision(2)
      << (double)stats.max_cache_size_bytes / (double)(1024ULL*1024ULL) << " MB");

  auto x = T();
  entry_base.push_back(x); // dummy Element
  stats.incorporate_cache_store(x.extra_bytes(), super_comp.nVars());

  entry_base.emplace_back(T(super_comp, hash_seed, bpc));
  stats.incorporate_cache_store(entry_base.back().extra_bytes(), super_comp.nVars());
  super_comp.set_id(1);
  compute_size_allocated();
}

template<typename T>
void CompCache<T>::test_descendantstree_consistency() {
  for (uint32_t id = 2; id < entry_base.size(); id++)
    if (!entry_base[id].is_free()) {
      CacheEntryID act_child = entry(id).first_descendant();
      while (act_child) {
        CacheEntryID next_child = entry(act_child).next_sibling();
        assert(entry(act_child).father() == id);
        act_child = next_child;
      }
      CacheEntryID father = entry(id).father();
      CacheEntryID act_sib = entry(father).first_descendant();

      [[maybe_unused]] bool found = false;
      while (act_sib) {
        CacheEntryID next_sib = entry(act_sib).next_sibling();
        if (act_sib == id) found = true;
        act_sib = next_sib;
      }
      assert(found);
    }
}

template<typename T>
uint64_t CompCache<T>::calc_cutoff() const {
  vector<uint64_t> scores;
  // skip index 0 (dummy) and 1 (root formula, never deletable)
  for (auto it = entry_base.begin() + 2; it != entry_base.end(); ++it)
    if (!it->is_free() && it->is_deletable()) scores.push_back(it->last_used_time());
  assert(!scores.empty());
  verb_print(1, "deletable:           " << scores.size());
  const auto mid_idx = scores.size() / 2;
  std::nth_element(scores.begin(), scores.begin() + mid_idx, scores.end());
  return scores[mid_idx];
}

template<typename T>
bool CompCache<T>::delete_some_entries() {
  // Check for deletable entries before doing anything — the cache may be full
  // purely due to infrastructure (e.g. table) with no evictable entries yet.
  const bool has_deletable = std::any_of(entry_base.begin() + 2, entry_base.end(),
    [](const auto& e) { return !e.is_free() && e.is_deletable(); });
  if (!has_deletable) return false;

  const auto start_del_time = cpu_time();
  uint64_t cutoff = (uint64_t)calc_cutoff();
  verb_print(1, "Deleting entries. Num entries: " << entry_base.size());
  verb_print(1, "cache_bytes_memory_usage() in MB: " << in_mb(stats.cache_bytes_memory_usage()));
  verb_print(1, "max_cache_size_bytes in MB: " << in_mb(stats.max_cache_size_bytes));
  verb_print(1, "free entries before: " << free_entry_base_slots.size());

  // note we start at index 2, since index 1 is the whole formula, should always stay here!
  uint64_t tot = 0;
  int64_t num = 0;
  uint64_t desc = 0;
  uint64_t max_desc = 0;
  uint64_t siblings = 0;
  uint64_t max_siblings = 0;
  // NOTE: between here and the rehash_table() call below, the hash table contains
  // stale entries for erased IDs. Do NOT call find_comp_and_incorporate_cnt() in this
  // window. The full rehash at the end restores consistency.
  for (uint32_t id = 2; id < entry_base.size(); id++)
    if (!entry_base[id].is_free() && entry_base[id].is_deletable() &&
        // conf.lru_eviction=1: evict least-recently-used
        // conf.lru_eviction=0: evict most-recently-used
        (conf.lru_eviction ? entry_base[id].last_used_time() <= cutoff
                           : entry_base[id].last_used_time() >= cutoff)) {
      auto d = num_direct_children(id);
      max_desc = std::max(max_desc, d);
      auto s = num_siblings(id);
      max_siblings = std::max(max_siblings, s);
      if (d < 100 && s < 200) {
        tot += unlink_from_tree(id);
        num++;
        desc += d;
        siblings += s;
        erase(id);
      }
    }
  verb_print(1, "max descendants: " << max_desc << " avg descendants: " << safe_div((double)desc, (double)num));
  verb_print(1, "max siblings: " << max_siblings << " avg siblings: " << safe_div((double)siblings, (double)num));
  verb_print(1, "free entries after:  " << free_entry_base_slots.size() << " avg len: " << safe_div((double)tot, (double)num));

  SLOW_DEBUG_DO(test_descendantstree_consistency());
  rehash_table(table.size());

  // Recompute mem usage — must start from 0 to include entry_base[1] (the root
  // formula entry, which has non-zero extra_bytes() from its packed comp data).
  stats.sum_extra_bytes = std::accumulate(entry_base.begin(), entry_base.end(), 0ULL,
      [](uint64_t s, const auto& e) { return e.is_free() ? s : s + e.extra_bytes(); });
  compute_size_allocated();
  verb_print(1, "deletion done. T: " << cpu_time()-start_del_time);
  return num > 0;
}

template<typename T>
uint64_t CompCache<T>::compute_size_allocated() {
  stats.cache_infra_bytes_mem_usage =
      sizeof(CompCache)
      + sizeof(CacheEntryID) * table.capacity()
      + sizeof(CacheEntryID) * free_entry_base_slots.capacity()
      + sizeof(T)* entry_base.capacity();
  return stats.cache_infra_bytes_mem_usage;
}

template<typename T>
void CompCache<T>::debug_mem_data() const {
    if (conf.verb < 2) return;
    cout << std::setw(40) << "c o sizeof (CacheableComp<T> + CacheEntryID) "
         << sizeof(T) << ", "
         << sizeof(CacheEntryID) << endl;
    cout << std::setw(40) << "c o table (size/capa) M " << table.size()/(double)1e6
         << "/" << table.capacity()/(double)1e6 << endl;
    cout << std::setw(40) << "c o entry_base (size/capa) M " << entry_base.size()/(double)1e6
         << "/" << entry_base.capacity()/(double)1e6 << endl;
    cout << std::setw(40) << "c o free_entry_base_slots (size/capa) M "
         << free_entry_base_slots.size()/(double)1e6
         << "/" << free_entry_base_slots.capacity()/(double)1e6 << endl;
    cout << std::setw(40) << "c o table mem use MB: "
         << (double)(table.capacity()*sizeof(CacheEntryID))/(double)(1024*1024) << endl;
    cout << std::setw(40) << "c o entry_base mem use MB: "
         << (double)(entry_base.capacity()*sizeof(T))/(double)(1024*1024) << endl;
    cout << std::setw(40) << "c o free_entry_base_slots mem use MB "
         << (double)(free_entry_base_slots.capacity()*sizeof(CacheEntryID))/(double)(1024*1024)
         << endl;
    uint64_t tot_extra_bytes = std::accumulate(entry_base.begin(), entry_base.end(), 0ULL,
        [](uint64_t s, const auto& e) { return e.is_free() ? s : s + e.extra_bytes(); });
    cout << std::setw(40) << "c o bignum(+packed comp if used) uses MB "
      << (double)tot_extra_bytes/(double)(1024*1024) << endl;

    double vm_dat;
    auto dat = mem_used(vm_dat);
    verb_print(1, "Total process MB : " << (double)dat/(double)(1024*1024)
      << " Total process vm MB: " << vm_dat/(double)(1024*1024));
}

template<typename T>
uint64_t CompCache<T>::num_direct_children(CacheEntryID id) {
  uint64_t ret = 0;
  for (auto c = entry(id).first_descendant(); c; c = entry(c).next_sibling()) ret++;
  return ret;
}

template<typename T>
uint64_t CompCache<T>::num_siblings(CacheEntryID id) {
  uint64_t ret = 0;
  CacheEntryID father = entry(id).father();
  if (entry(father).first_descendant() == id) return 0;
  CacheEntryID act_sibl = entry(father).first_descendant();
  while (act_sibl) {
    ret++;
    CacheEntryID next_sibl = entry(act_sibl).next_sibling();
    if (next_sibl == id) return ret;
    act_sibl = next_sibl;
  }
  release_assert(false);
  return ret;
}

// Used only during cache freeing. Unlinks from descendants tree
template<typename T>
uint64_t CompCache<T>::unlink_from_tree(CacheEntryID id) {
  assert(exists(id));
  // we need a father for this all to work
  assert(entry(id).father());
  assert(exists(entry(id).father()));
  stats.num_cache_dels++;

  // unlink id from the father's siblings list
  uint64_t len = 0;
  CacheEntryID father = entry(id).father();
  if (entry(father).first_descendant() == id) {
    entry(father).set_first_descendant(entry(id).next_sibling());
  } else {
    CacheEntryID act_sibl = entry(father).first_descendant();
    while (act_sibl) {
      len++;
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
  return len;
}

}
