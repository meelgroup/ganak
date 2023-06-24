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

#include "comp_cache.h"
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
  sysctl(mib, 2, &physical_memory, &length, NULL, 0);

  return physical_memory;
}

#else

#endif

#include "stack.h"

ComponentCache::ComponentCache(
    DataAndStatistics &_stats, const CounterConfiguration &_conf, const BPCSizes& _sz) :
    stats(_stats), conf(_conf), sz(_sz) {
}

void ComponentCache::init(Component &super_comp, void* randomseedforCLHASH){
  CacheableComponent *packed_super_comp;
  vector<uint32_t> tmp(100+super_comp.nVars()+super_comp.numLongClauses());
  packed_super_comp = new CacheableComponent(randomseedforCLHASH,super_comp, sz, tmp.data());
  my_time_ = 1;

  entry_base_.clear();
  auto x = CacheableComponent();
  entry_base_.push_back(x); // dummy Element
  stats.incorporate_cache_store(x, 0);
  table_.clear();
  table_.resize(1024*1024, 0);
  table_size_mask_ = table_.size() - 1;

  free_entry_base_slots_.clear();

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
  verb_print(1, "c Max cache size (80% free mem-200MB): "
    << stats.maximum_cache_size_bytes_ / (1024ULL*1024ULL) << " MB");

  assert(!cache_full());
  entry_base_.push_back(*packed_super_comp);
  stats.incorporate_cache_store(*packed_super_comp, 0);
  delete packed_super_comp;
  super_comp.set_id(1);
  compute_size_allocated();
}

void ComponentCache::test_descendantstree_consistency() {
  for (uint32_t id = 2; id < entry_base_.size(); id++)
    if (!entry_base_[id].is_free()) {
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

bool ComponentCache::deleteEntries() {
  assert(cache_full());
  vector<double> scores;
  verb_print(1, "Deleting entires. Num entries: " << entry_base_.size());
  verb_print(1, "cache_bytes_memory_usage() in MB: " << (stats.cache_bytes_memory_usage())/(1024ULL*1024ULL));
  verb_print(1, "maximum_cache_size_bytes_ in MB: " << (stats.maximum_cache_size_bytes_)/(1024ULL*1024ULL));
  verb_print(1, "free entries before: " << free_entry_base_slots_.size());
  for (auto it = entry_base_.begin() + 1; it != entry_base_.end(); it++)
    if (!it->is_free() && (it)->isDeletable()) {
      scores.push_back((double) (it)->creation_time());
    }
  if (scores.empty()){
    cout<< "c ERROR Memory out!"<<endl;
    exit(-1);
    assert(!scores.empty());
  }
  verb_print(1, "deletable:           " << scores.size())
  sort(scores.begin(), scores.end());
  double cutoff = scores[scores.size() / 2];
  // first : go through the EntryBase and mark the entries to be deleted as deleted (i.e. EMPTY
  // note we start at index 2,
  // since index 1 is the whole formula,
  // should always stay here!
  for (uint32_t id = 2; id < entry_base_.size(); id++)
    if (!entry_base_[id].is_free() &&
        entry_base_[id].isDeletable() &&
        (double) entry_base_[id].creation_time() <= cutoff) {
      removeFromDescendantsTree(id);
      eraseEntry(id);
    }
  verb_print(1, "free entries after:  " << free_entry_base_slots_.size());

  // then go through the Hash Table and erase all Links to empty entries
  SLOW_DEBUG_DO(test_descendantstree_consistency());
  reHashTable(table_.size());
  stats.sum_bytes_cached_comps_ = 0;

  for (uint32_t id = 2; id < entry_base_.size(); id++)
    if (!entry_base_[id].is_free()) {
      stats.sum_bytes_cached_comps_ += entry_base_[id].SizeInBytes();
    }

  stats.num_cached_comps_ = entry_base_.size();
  compute_size_allocated();
  return true;
}

uint64_t ComponentCache::compute_size_allocated() {
  stats.cache_infrastructure_bytes_memory_usage_ =
      sizeof(ComponentCache)
      + sizeof(CacheEntryID)* table_.capacity()
      + sizeof(CacheableComponent)* entry_base_.capacity()
      + sizeof(CacheEntryID) * free_entry_base_slots_.capacity();
  return stats.cache_infrastructure_bytes_memory_usage_;
}

void ComponentCache::debug_dump_data() {
    cout << "sizeof (CacheableComponent, CacheEntryID) "
         << sizeof(CacheableComponent) << ", "
         << sizeof(CacheEntryID) << endl;
    cout << "table (size/capacity) " << table_.size()
         << "/" << table_.capacity() << endl;
    cout << "entry_base_ (size/capacity) " << entry_base_.size()
             << "/" << entry_base_.capacity() << endl;
    cout << "free_entry_base_slots_ (size/capacity) " << free_entry_base_slots_.size()
             << "/" << free_entry_base_slots_.capacity() << endl;

    cout << "-" << endl;
    cout << std::setw(40) << "table mem use MB: " <<
      (double)(table_.capacity()*sizeof(CacheableComponent))/(double)(1024*1024)
      << endl;
    cout << std::setw(40) << "entry_base_ mem use MB: " <<
      (double)(entry_base_.capacity()*sizeof(CacheEntryID))/(double)(1024*1024)
      << endl;
    cout << std::setw(40) << "free_entry_base_slots_ mem use MB " <<
      (double)(free_entry_base_slots_.capacity()*sizeof(uint32_t))/(double)(1024*1024)
      << endl;

    uint64_t alloc_model_counts = 0;
    for (auto &pentry : entry_base_)
      if (!pentry.is_free()) alloc_model_counts += pentry.alloc_of_model_count();
    cout << "model counts size " << alloc_model_counts << endl;
}
