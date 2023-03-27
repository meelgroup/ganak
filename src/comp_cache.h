/*
 * comp_cache.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#ifndef COMPONENT_CACHE_H_
#define COMPONENT_CACHE_H_

#include "statistics.h"
#include "solver_config.h"
#include <gmpxx.h>

#include "comp_types/comp.h"

#include "stack.h"


// There is EXACTLY ONE of this
class ComponentCache {
public:
  ComponentCache(DataAndStatistics &statistics, const SolverConfiguration& config);

  ~ComponentCache() {
   // debug_dump_data();
    for (auto &pentry : entry_base_)
          if (pentry != nullptr)
            delete pentry;
  }

  void init(Component &super_comp, vector <void*>  &randomseedforCLHASH);

  // compute the size in bytes of the comp cache from scratch
  // the value is stored in bytes_memory_usage_
  uint64_t compute_size_used();

  CacheableComponent &entry(CacheEntryID id) {
    assert(entry_base_.size() > id);
    assert(entry_base_[id] != nullptr);
    return *entry_base_[id];
  }

  CacheableComponent &entry(const Component& comp) {
      return entry(comp.id());
  }

  bool hasEntry(CacheEntryID id) {
    assert(entry_base_.size() > id);
    return entry_base_[id];
  }

  // removes the entry id from the hash table
  // but not from the entry base
  inline void removeFromHashTable(CacheEntryID id);

  // we delete the Component with ID id
  // and all its descendants from the cache
  inline void cleanPollutionsInvolving(CacheEntryID id);

  // creates a CCacheEntry in the entry base
  // which contains a packed copy of comp
  // returns the id of the entry created
  // stores in the entry the position of
  // comp which is a part of the comp stack
  inline CacheEntryID storeAsEntry(CacheableComponent &ccomp,
                            CacheEntryID super_comp_id);

  // check quickly if the model count of the comp is cached
  // if so, incorporate it into the model count of top
  // if not, store the packed version of it in the entry_base of the cache
  // bool manageNewComponent(StackLevel &top, CacheableComponent &packed_comp) {
  //   stats.num_cache_look_ups_++;
  //   uint32_t table_ofs =  packed_comp.hashkey() & table_size_mask_;

  //   CacheEntryID act_id = table_[table_ofs];
  //   while(act_id){
  //     if (entry(act_id).equals(packed_comp)) {
  //       stats.incorporate_cache_hit(packed_comp);
  //       top.includeSolution(entry(act_id).model_count());
  //       return true;
  //     }
  //     act_id = entry(act_id).next_bucket_element();
  //   }
  //   return false;
  // }

  bool manageNewComponent(StackLevel &top, CacheableComponent &packed_comp) {
    stats.num_cache_look_ups_++;
    uint64_t *clhash_key;
    uint32_t table_ofs =  packed_comp.hashkey() & table_size_mask_;
    CacheEntryID act_id = table_[table_ofs];
    if (config_.do_pcc) {
#ifdef DOPCC
      if (!act_id) return false;
      clhash_key = packed_comp.compute_clhash();
      while(act_id){
        if (entry(act_id).equals(packed_comp, clhash_key)) {
          stats.incorporate_cache_hit(packed_comp);
          top.includeSolution(entry(act_id).model_count());
          return true;
        }
        act_id = entry(act_id).next_bucket_element();
      }
      return false;
#else
      assert(false);
      exit(-1);
#endif
    } else{
      while(act_id){
        if (entry(act_id).equals(packed_comp)) {
          stats.incorporate_cache_hit(packed_comp);
          top.includeSolution(entry(act_id).model_count());
          return true;
        }
        act_id = entry(act_id).next_bucket_element();
      }
      return false;
    }
  }

  // unchecked erase of an entry from entry_base_
  void eraseEntry(CacheEntryID id) {
    stats.incorporate_cache_erase(*entry_base_[id], config_.do_pcc && entry_base_[id]->get_hacked());
    delete entry_base_[id];
    entry_base_[id] = nullptr;
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
  void reHashTable(uint32_t size){

    table_.clear();
    table_.resize(size,0);
    // we assert that table size is a power of 2
    // otherwise the table_size_mask_ doesn't work
    assert((table_.size() & (table_.size() - 1)) == 0);
    table_size_mask_ = table_.size() - 1;
    for (uint32_t id = 2; id < entry_base_.size(); id++)
      if (entry_base_[id] != nullptr ){
        entry_base_[id]->set_next_bucket_element(0);
       if(entry_base_[id]->modelCountFound()) {
        uint32_t table_ofs=tableEntry(id);
        entry_base_[id]->set_next_bucket_element(table_[table_ofs]);
        table_[table_ofs] = id;
       }
    }
  }

  uint32_t tableEntry(CacheEntryID id){
    return entry(id).hashkey() & table_size_mask_;
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

  vector<CacheableComponent *> entry_base_;
  vector<CacheEntryID> free_entry_base_slots_;

  // the actual hash table
  // by means of which the cache is accessed
  vector<CacheEntryID> table_;

  uint32_t table_size_mask_;

  DataAndStatistics &stats;
  const SolverConfiguration &config_;

  uint64_t my_time_ = 0;
};


#include "comp_cache-inl.h"


#endif /* COMPONENT_CACHE_H_ */
