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
#include <atomic>
#include <type_traits>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#include <gmpxx.h>
#include "structures.hpp"
#include "time_mem.hpp"
#include "comp_cache_if.hpp"
#include "cache_table.hpp"

#include "comp_types/comp.hpp"
#include "stack.hpp"
namespace GanakInt {

inline std::atomic<int> num_live_caches{0};
inline std::atomic<int> num_parallel_counts{0}; // OuterCounter::count_with_parallel() is running

// Bytes malloc has handed out, main arena only. 0 if we cannot tell.
inline uint64_t heap_in_use() {
// Nested: without glibc, __GLIBC_PREREQ is undefined and the && does not save the #if
#if defined(__GLIBC__)
#if __GLIBC_PREREQ(2, 33)
  const auto mi = mallinfo2();
  return mi.uordblks + mi.hblkhd;
#else
  return 0;
#endif
#else
  return 0;
#endif
}


// There is EXACTLY ONE of this. Entries are kept in creation order with no
// holes: eviction compacts (and renumbers the comps on the stack), and cache
// pollution removal truncates.
template<typename T>
class CompCache final: public CompCacheIF {
  static_assert(std::is_trivially_copyable_v<T>);
  static constexpr double max_load = 0.85;
  static constexpr double rebuild_load = 0.5;

public:
  CompCache(DataAndStatistics &_stats, const CounterConfiguration &_conf, const CMSat::FieldGen& fg) :
      stats(_stats), conf(_conf), scratch(fg.zero()),
      int_inline(SmallIntCodec::supported && dynamic_cast<const ArjunNS::FGenMpz*>(&fg) != nullptr) {
    num_live_caches++;
  }
  ~CompCache() override {
    for (auto& e : entry_base) e.release();
    probe.release();
    num_live_caches--;
  }

  void init(Comp &super_comp, uint32_t hash_seed, const BPCSizes& bpc,
      const vector<Comp*>* comp_stack) override;
  void set_probe(const Comp &comp, uint32_t hash_seed, const BPCSizes& bpc) override {
    probe.set_comp(comp, hash_seed, bpc);
    table.prefetch(probe.get_full_hashkey());
  }
  bool find_probe_and_incorporate_cnt(StackLevel &top, uint32_t nvars, int* out_node = nullptr) override;
  CacheEntryID add_probe(uint32_t nvars, const Comp& super_comp) override;
  void set_compile_node(CacheEntryID id, int node) override {
    if (id >= compile_nodes.size()) compile_nodes.resize(id + 1, -1);
    compile_nodes[id] = node;
  }
  void make_entry_deletable(CacheEntryID id) override { entry(id).set_deletable(); }
  [[nodiscard]] uint64_t get_num_entries_used() const override { return entry_base.size() - 2; }
  [[nodiscard]] uint64_t get_max_num_entries() const override { return entry_base.size(); }
  [[nodiscard]] bool exists(CacheEntryID id) const override { return id != 0 && id < entry_base.size(); }
  uint64_t remove_from(CacheEntryID from) override;
  void store_value(CacheEntryID id, const FF& model_count) override;
  void debug_mem_data() const override;
#ifdef SLOW_DEBUG
  void check_removal(CacheEntryID from, const vector<CacheEntryID>& roots) const override;
#endif

private:
  T& entry(CacheEntryID id) { return entry_base[id]; }
  const T& entry(CacheEntryID id) const { return entry_base[id]; }
  void update_entry_time(CacheEntryID id) {
    switch(conf.cache_time_update) {
      case 0: break;
      case 1: entry(id).set_last_used_time(my_time); break;
      case 2: entry(id).avg_last_used_time(my_time, 2); break;
      case 3: entry(id).avg_last_used_time(my_time, 3); break;
      default: release_assert(false);
    }
  }
  [[nodiscard]] bool cache_full(uint64_t extra_will_be_added) const {
    // allocated but unused entry_base slots are not live memory
    return stats.cache_full((entry_base.capacity() - entry_base.size()) * sizeof(T), extra_will_be_added);
  }
  [[nodiscard]] uint64_t calc_extra_mem_after_push() const;
  void compute_size_allocated();
  [[nodiscard]] bool calc_cutoff(uint64_t& cutoff) const;
  bool delete_some_entries();
  void rebuild_table();
  void remap_stack_ids();
  void check_freed_mem(uint64_t heap_before, uint64_t accounted_before) const;

  vec<T> entry_base;
  CacheTable table;
  T probe;
  const vector<Comp*>* comp_stack = nullptr;
  // eviction: (old id, new id) of the entries that are not deletable, i.e. on the comp stack
  vector<std::pair<CacheEntryID, CacheEntryID>> moved;
  // d-DNNF only: CacheEntryID -> compiled sub-DAG node id. Lazily sized; -1 = none.
  std::vector<int> compile_nodes;
#ifdef SLOW_DEBUG
  // the comp each entry was found in, with evicted fathers skipped
  vector<CacheEntryID> dbg_father;
#endif

  DataAndStatistics &stats;
  const CounterConfiguration &conf;
  FF scratch; // inline counts are unpacked into this on a hit
  const bool int_inline;
  uint64_t my_time = 0;
};

template<typename T>
void CompCache<T>::init(Comp &super_comp, uint32_t hash_seed, const BPCSizes& bpc,
    const vector<Comp*>* _comp_stack) {
  assert(entry_base.empty());
  comp_stack = _comp_stack;
  my_time = 1;
  stats.sum_extra_bytes = 0;
  stats.cache_infra_bytes_mem_usage = 0;
  table.reset(1ULL << 16, rebuild_load);
  verb_print(2, "Max cache size set: " << std::setprecision(2)
      << (double)stats.max_cache_size_bytes / (double)(1024ULL*1024ULL) << " MB"
      << " inline int counts: " << int_inline);

  entry_base.push_back(T()); // dummy
  stats.incorporate_cache_store(0, super_comp.nVars());

  T root;
  root.set_comp(super_comp, hash_seed, bpc);
  root.reset_meta(my_time);
  entry_base.push_back(root);
  stats.incorporate_cache_store(root.extra_bytes(), super_comp.nVars());
  super_comp.set_id(1);
  SLOW_DEBUG_DO(dbg_father.assign(2, 0));
  compute_size_allocated();
}

template<typename T>
bool CompCache<T>::find_probe_and_incorporate_cnt(StackLevel &top, const uint32_t nvars, int* out_node) {
  stats.num_cache_look_ups++;
  const HashVal& h = probe.get_full_hashkey();
  const CacheEntryID id = table.find(h, [&](const CacheEntryID i) {
    const T& e = entry(i);
    if (e.get_full_hashkey() != h) return false;
    if (!e.content_differs(probe)) return true;
    verb_print(1, "WARNING: hash collision -- cache entry ID " << i
         << " has the same hashkey (0x" << std::hex << h << std::dec
         << ") as the looked-up component but DIFFERENT content."
         << " Under --prob 1 this would have been a wrong cache hit.");
    return false;
  });
  if (!id) return false;
  stats.incorporate_cache_hit(nvars);
  update_entry_time(id);
  debug_print(COLYEL2 << "Cache hit. cache ID: " << id);
  top.include_solution(entry(id).model_count(*scratch));
  if (out_node != nullptr) *out_node = (id < compile_nodes.size()) ? compile_nodes[id] : -1;
  return true;
}

template<typename T>
uint64_t CompCache<T>::calc_extra_mem_after_push() const {
  uint64_t extra = 0;
  // entry_base grows by ~50% on realloc. The packed data and counts of the
  // entries are not copied, so they are not part of the growth
  if (entry_base.capacity() == entry_base.size()) extra += entry_base.capacity() * sizeof(T) / 2;
  if (table.over_load(max_load)) {
    const uint64_t now = table.bytes();
    const uint64_t after = CacheTable::bytes_for(table.num_live() + 1, rebuild_load);
    if (after > now) extra += after - now;
  }
  return extra;
}

template<typename T>
CacheEntryID CompCache<T>::add_probe(const uint32_t nvars, [[maybe_unused]] const Comp& super_comp) {
  uint64_t extra_mem_with_push = calc_extra_mem_after_push();
  while (cache_full(extra_mem_with_push)) {
    verb_print(1, "Cache full. Deleting some entries.");
    debug_mem_data();
    const bool freed_any = delete_some_entries();
    debug_mem_data();
    extra_mem_with_push = calc_extra_mem_after_push();
    // nothing deletable: accept a slightly over-budget allocation
    if (!freed_any) break;
  }

  const bool at_capacity = entry_base.capacity() == entry_base.size();
  probe.reset_meta(my_time++);
  entry_base.push_back(probe);
  probe.forget();
  const CacheEntryID id = entry_base.size() - 1;
  if (at_capacity) compute_size_allocated();
  stats.incorporate_cache_store(entry(id).extra_bytes(), nvars);
#ifdef SLOW_DEBUG
  // super_comp is on the stack, so an eviction above renumbered it
  assert(super_comp.id() < id && exists(super_comp.id()));
  assert(!entry(super_comp.id()).is_deletable());
  dbg_father.push_back(super_comp.id());
  assert(dbg_father.size() == entry_base.size());
#endif
  return id;
}

template<typename T>
uint64_t CompCache<T>::remove_from(const CacheEntryID from) {
  assert(from >= 2 && from <= entry_base.size());
  const uint64_t n = entry_base.size() - from;
  for (CacheEntryID id = entry_base.size() - 1; id >= from; id--) {
    T& e = entry(id);
    if (e.model_count_found()) table.erase(e.get_full_hashkey(), id);
    stats.incorporate_cache_erase(e.extra_bytes());
    e.release();
  }
  entry_base.shrink_(n);
  if (compile_nodes.size() > from) compile_nodes.resize(from);
  SLOW_DEBUG_DO(dbg_father.resize(from));
  VERBOSE_DEBUG_DO(cout << "Pollution removal from ID: " << from << " removed: " << n << endl);
  return n;
}

#ifdef SLOW_DEBUG
// Removing everything from `from` on must remove exactly the roots and their
// descendants: every other removed entry descends from a removed one, and a
// descendant always has a higher id than its father
template<typename T>
void CompCache<T>::check_removal(const CacheEntryID from, const vector<CacheEntryID>& roots) const {
  for (CacheEntryID id = from; id < entry_base.size(); id++) {
    const bool is_root = std::find(roots.begin(), roots.end(), id) != roots.end();
    if (is_root) continue;
    assert(dbg_father[id] >= from && "removing an entry that is not a descendant");
    assert(entry(id).is_deletable() && "removing a non-root entry that is on the stack");
  }
}
#endif

template<typename T>
void CompCache<T>::store_value(const CacheEntryID id, const FF& model_count) {
  assert(exists(id));
  if (table.over_load(max_load)) {
    rebuild_table();
    compute_size_allocated();
  }
  T& e = entry(id);
  stats.sum_extra_bytes -= e.extra_bytes();
  e.set_model_count(model_count, int_inline);
  e.set_last_used_time(my_time);
  stats.sum_extra_bytes += e.extra_bytes();
  stats.cache_inline_cnts += e.count_is_inline();
  if (e.model_count_found()) table.insert(e.get_full_hashkey(), id);
}

template<typename T>
void CompCache<T>::rebuild_table() {
  uint64_t n = 0;
  for (CacheEntryID id = 2; id < entry_base.size(); id++) n += entry(id).model_count_found();
  table.reset(n + 1, rebuild_load);
  // the inserts are random, keep many of them in flight
  constexpr CacheEntryID ahead = 16;
  const CacheEntryID sz = entry_base.size();
  for (CacheEntryID id = 2; id < sz; id++) {
    if (id + ahead < sz) table.prefetch(entry(id + ahead).get_full_hashkey());
    if (entry(id).model_count_found()) table.insert(entry(id).get_full_hashkey(), id);
  }
}

template<typename T>
bool CompCache<T>::calc_cutoff(uint64_t& cutoff) const {
  // The median of a fixed-stride sample: collecting every entry's time took
  // 8 bytes per entry, a transient well above the cache budget.
  constexpr uint64_t max_samples = 1ULL << 18;
  vector<uint64_t> scores;
  uint64_t stride = std::max<uint64_t>(1, entry_base.size() / max_samples);
  while (true) {
    // skip index 0 (dummy) and 1 (root formula, never deletable)
    for (uint64_t id = 2; id < entry_base.size(); id += stride)
      if (entry_base[id].is_deletable()) scores.push_back(entry_base[id].last_used_time());
    if (!scores.empty() || stride == 1) break;
    stride = 1;
  }
  if (scores.empty()) return false;
  verb_print(1, "deletable sampled:   " << scores.size() << " stride: " << stride);
  const auto mid_idx = scores.size() / 2;
  std::nth_element(scores.begin(), scores.begin() + mid_idx, scores.end());
  cutoff = scores[mid_idx];
  return true;
}

template<typename T>
bool CompCache<T>::delete_some_entries() {
  uint64_t cutoff;
  // the cache may be full purely due to infrastructure, with nothing evictable yet
  if (!calc_cutoff(cutoff)) return false;

  const auto start_del_time = cpu_time();
  verb_print(1, "Deleting entries. Num entries: " << entry_base.size());
  verb_print(1, "cache_bytes_memory_usage() in MB: " << in_mb(stats.cache_bytes_memory_usage()));
  verb_print(1, "max_cache_size_bytes in MB: " << in_mb(stats.max_cache_size_bytes));
  const uint64_t heap_before = heap_in_use();
  compute_size_allocated();
  const uint64_t accounted_before = stats.cache_bytes_memory_usage();

#ifdef SLOW_DEBUG
  // old id -> new id, an evicted entry maps to where its father went
  vector<CacheEntryID> remap(entry_base.size());
  remap[0] = 0;
  remap[1] = 1;
#endif
  moved.clear();
  uint64_t num = 0;
  CacheEntryID w = 2;
  for (CacheEntryID r = 2; r < entry_base.size(); r++) {
    T& e = entry(r);
    // conf.lru_eviction=1: evict least-recently-used
    // conf.lru_eviction=0: evict most-recently-used
    if (e.is_deletable() && (conf.lru_eviction ? e.last_used_time() <= cutoff
                                               : e.last_used_time() >= cutoff)) {
      stats.incorporate_cache_erase(e.extra_bytes());
      e.release();
      stats.num_cache_dels++;
      num++;
      SLOW_DEBUG_DO(remap[r] = remap[dbg_father[r]]);
      continue;
    }
    if (!e.is_deletable()) moved.emplace_back(r, w);
    if (w != r) {
      entry_base[w] = e;
      if (w < compile_nodes.size()) compile_nodes[w] = r < compile_nodes.size() ? compile_nodes[r] : -1;
    }
#ifdef SLOW_DEBUG
    remap[r] = w;
    dbg_father[w] = remap[dbg_father[r]];
#endif
    w++;
  }
  entry_base.shrink_(entry_base.size() - w);
  if (compile_nodes.size() > w) compile_nodes.resize(w);
  SLOW_DEBUG_DO(dbg_father.resize(w));
  verb_print(1, "evicted: " << num << " left: " << entry_base.size());

  rebuild_table();
  remap_stack_ids();

  // must start from 0 to include entry_base[1], the root formula
  stats.sum_extra_bytes = std::accumulate(entry_base.begin(), entry_base.end(), 0ULL,
      [](uint64_t s, const auto& e) { return s + e.extra_bytes(); });
  compute_size_allocated();
  check_freed_mem(heap_before, accounted_before);
  verb_print(1, "deletion done. T: " << cpu_time()-start_del_time);
  return num > 0;
}

template<typename T>
void CompCache<T>::remap_stack_ids() {
  for (Comp* c : *comp_stack) {
    if (c == nullptr || c->id() < 2) continue;
    const auto it = std::lower_bound(moved.begin(), moved.end(), std::make_pair(c->id(), CacheEntryID(0)));
    release_assert(it != moved.end() && it->first == c->id() && "comp on the stack has no cache entry");
    c->set_id(it->second);
  }
}

// The cache budget is only as good as extra_bytes(). An eviction frees nothing
// but cache memory, so malloc must agree with what we think we just freed.
template<typename T>
void CompCache<T>::check_freed_mem(const uint64_t heap_before, const uint64_t accounted_before) const {
  // mallinfo2 is process-wide: other threads would pollute it, also those
  // that have no cache yet
  if (num_live_caches != 1 || num_parallel_counts != 0) return;
  const uint64_t heap_after = heap_in_use();
  const uint64_t accounted_after = stats.cache_bytes_memory_usage();
  if (heap_before <= heap_after) return;
  const uint64_t real = heap_before - heap_after;
  const uint64_t accounted = accounted_before > accounted_after ? accounted_before - accounted_after : 0;
  verb_print(2, "eviction freed MB: " << in_mb(real) << " accounted MB: " << in_mb(accounted));
  // page rounding of mmap'd vectors etc. is a few KB of noise
  if (real < 1024*1024 || accounted + accounted/8 >= real) return;
  verb_print(1, "WARNING: cache memory accounting is off, --maxcache will be exceeded."
      << " Eviction freed " << real << " bytes but only " << accounted << " were accounted for");
  assert(false && "extra_bytes()/bytes_used() undercounts the real heap usage");
}

template<typename T>
void CompCache<T>::compute_size_allocated() {
  stats.cache_infra_bytes_mem_usage =
      sizeof(CompCache)
      + table.bytes()
      + sizeof(T) * entry_base.capacity()
      + sizeof(moved[0]) * moved.capacity()
      + sizeof(int) * compile_nodes.capacity();
}

template<typename T>
void CompCache<T>::debug_mem_data() const {
  if (conf.verb < 2) return;
  cout << "c o sizeof(entry) " << sizeof(T) << endl;
  cout << "c o table MB: " << in_mb(table.bytes()) << " live: " << table.num_live() << endl;
  cout << "c o entry_base (size/capa) M " << entry_base.size()/(double)1e6
       << "/" << entry_base.capacity()/(double)1e6
       << " MB: " << in_mb(entry_base.capacity()*sizeof(T)) << endl;
  const uint64_t tot_extra_bytes = std::accumulate(entry_base.begin(), entry_base.end(), 0ULL,
      [](uint64_t s, const auto& e) { return s + e.extra_bytes(); });
  cout << "c o bignum(+packed comp if used) uses MB " << in_mb(tot_extra_bytes) << endl;

  double vm_dat;
  auto dat = mem_used(vm_dat);
  verb_print(1, "Total process MB : " << in_mb(dat) << " Total process vm MB: " << in_mb(vm_dat));
}

}
