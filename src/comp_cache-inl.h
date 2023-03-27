/*
 * comp_cache-inl.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#pragma once

#include "primitive_types.h"
#include "comp_cache.h"

CacheEntryID ComponentCache::storeAsEntry(CacheableComponent &ccomp, CacheEntryID super_comp_id){
    CacheEntryID id;

    while (stats.cache_full()){
      if (config_.verb) cout << "c Cache full!!" << endl;
      deleteEntries();
    }

    assert(!stats.cache_full());

    ccomp.set_creation_time(my_time_++);

    if (free_entry_base_slots_.empty()) {
        if (entry_base_.capacity() == entry_base_.size()) {
            entry_base_.reserve(2 * entry_base_.size());
            compute_size_used();
        }
        entry_base_.push_back(&ccomp);
        id = entry_base_.size() - 1;
    } else {
        id = free_entry_base_slots_.back();
        assert(id < entry_base_.size());
        assert(entry_base_[id] == nullptr);
        free_entry_base_slots_.pop_back();
        entry_base_[id] = &ccomp;
    }

    entry(id).set_father(super_comp_id);
    add_descendant(super_comp_id, id);

    assert(hasEntry(id));
    assert(hasEntry(super_comp_id));

    stats.incorporate_cache_store(ccomp, config_.do_pcc && ccomp.get_hacked());

  #ifdef DEBUG
      for (uint32_t u = 2; u < entry_base_.size(); u++)
            if (entry_base_[u] != nullptr) {
              assert(entry_base_[u]->father() != id);
              assert(entry_base_[u]->first_descendant() != id);
              assert(entry_base_[u]->next_sibling() != id);
            }
  #endif
    return id;
}

void ComponentCache::cleanPollutionsInvolving(const CacheEntryID id) {
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
    cleanPollutionsInvolving(act_child);
  }
  removeFromHashTable(id);
  eraseEntry(id);
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
  considerCacheResize();
  uint32_t table_ofs = tableEntry(id);
  // when storing the new model count the size of the model count
  // and hence that of the comp will change
  stats.sum_bytes_cached_comps_ -= entry(id).SizeInBytes();
  stats.overall_bytes_comps_stored_ -= entry(id).SizeInBytes();

  stats.sys_overhead_sum_bytes_cached_comps_ -= entry(id).sys_overhead_SizeInBytes();
  stats.sys_overhead_overall_bytes_comps_stored_ -= entry(id).sys_overhead_SizeInBytes();
#ifdef DOPCC
  if (config_.do_pcc)
    entry(id).set_hacked(entry(id).SizeInBytes(), entry(id).nVars());
#endif

  entry(id).set_model_count(model_count,my_time_);
  entry(id).set_creation_time(my_time_);

  entry(id).set_next_bucket_element(table_[table_ofs]);
  table_[table_ofs] = id;

  if (config_.do_pcc){
    stats.sum_bytes_cached_comps_ += entry(id).SizeInBytes_CLHASH();
    stats.overall_bytes_comps_stored_ += entry(id).SizeInBytes_CLHASH();
  }
  else{
    stats.sum_bytes_cached_comps_ += entry(id).SizeInBytes();
    stats.overall_bytes_comps_stored_ += entry(id).SizeInBytes();
  }
  stats.sys_overhead_sum_bytes_cached_comps_ += entry(id).sys_overhead_SizeInBytes();
  stats.sys_overhead_overall_bytes_comps_stored_ += entry(id).sys_overhead_SizeInBytes();
}
