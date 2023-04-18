/*
 * comp_cache.cpp
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

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

ComponentCache::ComponentCache(DataAndStatistics &statistics, const CounterConfiguration &config, const BPCSizes& _sz) :
		stats(statistics), config_(config), sz(_sz) {
}

void ComponentCache::init(Component &super_comp, void* randomseedforCLHASH){
  CacheableComponent *packed_super_comp;
#ifdef DOPCC
	vector<uint32_t> tmp(100+super_comp.nVars()+super_comp.numLongClauses());
	packed_super_comp = new CacheableComponent(randomseedforCLHASH,super_comp, sz, tmp.data());
	packed_super_comp->finish_hashing(packed_super_comp->SizeInBytes(sz), packed_super_comp->nVars(sz));
#else
	packed_super_comp = new CacheableComponent(super_comp, sz);
#endif
	my_time_ = 1;

	entry_base_.clear();
	entry_base_.push_back(new CacheableComponent()); // dummy Element
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
		cout <<"c WARNING: Maximum cache size larger than free RAM available" << endl;
		cout << "c Free RAM " << std::setprecision(2)
			<< (double)free_ram / (1024.0*1024.0) << "MB" << endl;
	}
	if (config_.verb) cout << "c Max cache size (80% free mem): "
		<< stats.maximum_cache_size_bytes_ / (1024ULL*1024ULL) << " MB" << endl;

	assert(!stats.cache_full());

	if (entry_base_.capacity() == entry_base_.size())
		entry_base_.reserve(2 * entry_base_.size());

	entry_base_.push_back(packed_super_comp);

	stats.incorporate_cache_store(*packed_super_comp, sz);

	super_comp.set_id(1);
	compute_size_used();
}

void ComponentCache::test_descendantstree_consistency() {
	for (uint32_t id = 2; id < entry_base_.size(); id++)
		if (entry_base_[id] != nullptr) {
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
				if (act_sib == id)
					found = true;
				act_sib = next_sib;
			}
			assert(found);
		}
}

#ifndef DOPCC
void ComponentCache::delete_comps_with_vars(const set<uint32_t>& vars) {
	size_t num_deleted = 0;
	size_t orig_num = entry_base_.size();
	// note we start at index 2,
	// since index 1 is the whole formula,
	// should always stay here!
	for (uint32_t id = 2; id < entry_base_.size(); id++)
		if (entry_base_[id] != nullptr && entry_base_[id]->isDeletable()) {
		  DifferencePackedComponent* d = entry_base_[id];
		  if (d->contains_any_var(vars, sz)) {
		    removeFromDescendantsTree(id);
		    eraseEntry(id);
				num_deleted++;
			}
		}
	cout << "c Num deleted: " << num_deleted << " of: " << orig_num
		<< " percent: " << (double)num_deleted/(double)orig_num * 100.0 << "%" << endl;
}
#endif

bool ComponentCache::deleteEntries() {
  assert(stats.cache_full());
	vector<double> scores;
	cout << "c Deleting entires. Num entries: " << entry_base_.size() << endl;
	cout << "c cache_bytes_memory_usage() in MB: " << (stats.cache_bytes_memory_usage())/(1024ULL*1024ULL) << endl;
	cout << "c maximum_cache_size_bytes_ in MB: " << (stats.maximum_cache_size_bytes_)/(1024ULL*1024ULL) << endl;
	for (auto it = entry_base_.begin() + 1; it != entry_base_.end(); it++)
		if (*it != nullptr && (*it)->isDeletable()) {
			scores.push_back((double) (*it)->creation_time());
		}
	if (scores.empty()){
		cout<< "c Memory out!"<<endl;
		exit(-1);
		assert(!scores.empty());
	}
	sort(scores.begin(), scores.end());
	double cutoff = scores[scores.size() / 2];
	// first : go through the EntryBase and mark the entries to be deleted as deleted (i.e. EMPTY
	// note we start at index 2,
	// since index 1 is the whole formula,
	// should always stay here!
	for (uint32_t id = 2; id < entry_base_.size(); id++)
		if (entry_base_[id] != nullptr &&
		    entry_base_[id]->isDeletable() &&
		      (double) entry_base_[id]->creation_time() <= cutoff) {
				removeFromDescendantsTree(id);
				eraseEntry(id);
		}
	// then go through the Hash Table and erase all Links to empty entries


#ifdef DEBUG
	test_descendantstree_consistency();
#endif

	reHashTable(table_.size());
	stats.sum_size_cached_comps_ = 0;
	stats.sum_bytes_cached_comps_ = 0;
	 stats.sys_overhead_sum_bytes_cached_comps_ =0;


	for (uint32_t id = 2; id < entry_base_.size(); id++)
		if (entry_base_[id] != nullptr) {
			stats.sum_size_cached_comps_ += entry_base_[id]->nVars(sz);
			stats.sum_bytes_cached_comps_ += entry_base_[id]->SizeInBytes(sz);
			stats.sys_overhead_sum_bytes_cached_comps_ += entry_base_[id]->sys_overhead_SizeInBytes(sz);
		}

	stats.num_cached_comps_ = entry_base_.size();
	compute_size_used();
	return true;
}


uint64_t ComponentCache::compute_size_used() {
  stats.cache_infrastructure_bytes_memory_usage_ =
      sizeof(ComponentCache)
      + sizeof(CacheEntryID)* table_.capacity()
      + sizeof(CacheableComponent *)* entry_base_.capacity()
      + sizeof(CacheEntryID) * free_entry_base_slots_.capacity();
  return stats.cache_infrastructure_bytes_memory_usage_;
}

void ComponentCache::debug_dump_data() {
    cout << "sizeof (CacheableComponent *, CacheEntryID) "
         << sizeof(CacheableComponent *) << ", "
         << sizeof(CacheEntryID) << endl;
    cout << "table (size/capacity) " << table_.size()
         << "/" << table_.capacity() << endl;
    cout << "entry_base_ (size/capacity) " << entry_base_.size()
             << "/" << entry_base_.capacity() << endl;
    cout << "free_entry_base_slots_ (size/capacity) " << free_entry_base_slots_.size()
             << "/" << free_entry_base_slots_.capacity() << endl;

    uint64_t alloc_model_counts = 0;
    for (auto &pentry : entry_base_)
              if (pentry != nullptr){
                alloc_model_counts += pentry->alloc_of_model_count();
              }
    cout << "model counts size " << alloc_model_counts << endl;
}
