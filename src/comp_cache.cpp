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

#include "comp_cache.hpp"
#include "primitive_types.hpp"
#include "time_mem.hpp"
#include <algorithm>
#include <iomanip>

#ifdef __linux__

#include <sys/sysinfo.h>
#include <cstdint>

uint64_t freeram() {

  struct sysinfo info;
  sysinfo(&info);

  return info.freeram *(uint64_t) info.mem_unit;
}

#elif __APPLE__ && __MACH__

#include <sys/types.h>
#include <sys/sysctl.h>

uint64_t freeram() {
  int mib[2];
  int64_t physical_memory;
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  size_t length = sizeof(int64_t);
  sysctl(mib, 2, &physical_memory, &length, nullptr, 0);

  return physical_memory;
}

#else

#endif

#include "stack.hpp"

CompCache::CompCache(
    DataAndStatistics &_stats, const CounterConfiguration &_conf) : stats(_stats), conf(_conf) {}

void CompCache::init(Comp &super_comp, void* hash_seed){
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

  if (stats.maximum_cache_size_bytes_ == 0) {
    stats.maximum_cache_size_bytes_ = max_cache_bound;
  }

  if (stats.maximum_cache_size_bytes_ > free_ram) {
    verb_print(1, "WARNING: Maximum cache size larger than free RAM available");
    verb_print(1, "Free RAM " << std::setprecision(2)
      << (double)free_ram / (1024.0*1024.0) << "MB");
  }
  verb_print(1, "Max cache size (80% free mem-200MB): "
    << stats.maximum_cache_size_bytes_ / (1024ULL*1024ULL) << " MB");


  stats.sum_bignum_bytes = 0;
  stats.cache_infrastructure_bytes_memory_usage_ = 0;
  assert(!cache_full());
  entry_base.push_back(*packed_super_comp);
  stats.incorporate_cache_store(*packed_super_comp, super_comp.nVars());
  delete packed_super_comp;
  super_comp.set_id(1);
  compute_size_allocated();
}

void CompCache::test_descendantstree_consistency() {
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

      bool found = false;
      while (act_sib) {
        CacheEntryID next_sib = entry(act_sib).next_sibling();
        if (act_sib == id) found = true;
        act_sib = next_sib;
      }
      assert(found);
    }
}

double CompCache::calc_cutoff() const {
  vector<uint32_t> scores;
  // TODO: this score is VERY simplistic, we actually don't touch it at all, ever
  //       just create it and that's it. Not bumped with usage(!)
  for (auto it = entry_base.begin() + 1; it != entry_base.end(); it++)
    if (!it->is_free() && it->is_deletable()) {
      scores.push_back(it->last_used_time());
    }
  if (scores.empty()){
    cout<< "c ERROR Memory out!"<<endl;
    exit(-1);
    assert(!scores.empty());
  }
  verb_print(1, "deletable:           " << scores.size())
  sort(scores.begin(), scores.end());
  return scores[scores.size() / 2];
}

bool CompCache::delete_some_entries() {
  const auto start_del_time = cpuTime();
  double cutoff = calc_cutoff();
  verb_print(1, "Deleting entires. Num entries: " << entry_base.size());
  verb_print(1, "cache_bytes_memory_usage() in MB: " << (stats.cache_bytes_memory_usage())/(1024ULL*1024ULL));
  verb_print(1, "maximum_cache_size_bytes_ in MB: " << (stats.maximum_cache_size_bytes_)/(1024ULL*1024ULL));
  verb_print(1, "free entries before: " << free_entry_base_slots.size());

  // note we start at index 2, since index 1 is the whole formula, should always stay here!
  uint64_t tot = 0;
  int64_t num = 0;
  for (uint32_t id = 2; id < entry_base.size(); id++)
    if (!entry_base[id].is_free() &&
        entry_base[id].is_deletable() &&
        entry_base[id].get_dont_delete_before() < my_time &&
        ((!conf.do_cache_reverse_sort && entry_base[id].last_used_time() <= cutoff)
         || (conf.do_cache_reverse_sort && entry_base[id].last_used_time() >= cutoff))) {
      tot += unlink_from_tree(id);
      num++;
      erase(id); // Note: no need to incorporate erase, we recompute bignum bytes below
    }
  verb_print(1, "free entries after:  " << free_entry_base_slots.size() << " avg len: " << (double) tot/(double) num);

  SLOW_DEBUG_DO(test_descendantstree_consistency());
  rehash_table(table.size());

  // Recompute mem usage
  stats.sum_bignum_bytes = 0;
  for (uint32_t id = 2; id < entry_base.size(); id++)
    if (!entry_base[id].is_free()) {
      stats.sum_bignum_bytes += entry_base[id].bignum_bytes();
    }
  compute_size_allocated();

  stats.num_cached_comps_ = entry_base.size();
  verb_print(1, "deletion done. T: " << cpuTime()-start_del_time);
  return true;
}

uint64_t CompCache::compute_size_allocated() {
  stats.cache_infrastructure_bytes_memory_usage_ =
      sizeof(CompCache)
      + sizeof(CacheEntryID) * table.capacity()
      + sizeof(CacheEntryID) * free_entry_base_slots.capacity()
      + sizeof(CacheableComp)* entry_base.capacity();
  return stats.cache_infrastructure_bytes_memory_usage_;
}

void CompCache::debug_mem_data() const {
    cout << std::setw(40) << "c o sizeof (CacheableComp, CacheEntryID) "
         << sizeof(CacheableComp) << ", "
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
         << (double)(entry_base.capacity()*sizeof(CacheableComp))/(double)(1024*1024) << endl;
    cout << std::setw(40) << "c o free_entry_base_slots mem use MB "
         << (double)(free_entry_base_slots.capacity()*sizeof(uint32_t))/(double)(1024*1024)
         << endl;
    uint64_t tot_bignum_bytes = 0;
    for (auto &entry : entry_base)
      if (!entry.is_free()) tot_bignum_bytes += entry.bignum_bytes();
    cout << std::setw(40) << "c o bignum uses MB "
      << tot_bignum_bytes/(double)(1024*1024) << endl;


    double vm_dat;
    auto dat = memUsedTotal(vm_dat);
    verb_print(1, "Total process MB : " << dat/(double)(1024*1024)
      << " Total process vm MB: " << vm_dat/(double)(1024*1024));
}

// Used only during cache freeing. Unlinks from descendants tree
uint64_t CompCache::unlink_from_tree(CacheEntryID id) {
  assert(exists(id));
  // we need a father for this all to work
  assert(entry(id).father());
  assert(exists(entry(id).father()));
  stats.num_cache_dels_++;

  // unlink id from the father's siblings list
  uint64_t len = 0;
  CacheEntryID father = entry(id).father();
  if (entry(father).first_descendant() == id) {
    entry(father).set_first_descendant(entry(id).next_sibling());
  } else {
    CacheEntryID act_sibl = entry(father).first_descendant();
    while (act_sibl) {
    len ++;
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
