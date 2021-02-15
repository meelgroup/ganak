/*
 * component_cache.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#ifndef COMPONENT_CACHE_H_
#define COMPONENT_CACHE_H_

#include "component_types/cacheable_component.h"
#include "statistics.h"
#include "solver_config.h"
#include <gmpxx.h>

#include "component_types/component.h"

#include "stack.h"


class ComponentCache {
public:

    ComponentCache(DataAndStatistics &statistics, SolverConfiguration &config);


    ~ComponentCache() {
        // debug_dump_data();
        for (auto &pentry : entry_base_)
            if (pentry != nullptr)
                delete pentry;
    }

    void init(Component &super_comp);

    // compute the size in bytes of the component cache from scratch
    // the value is stored in bytes_memory_usage_
    uint64_t compute_byte_size_infrasture();

    CacheableComponent &entry(CacheEntryID id) {
        assert(entry_base_.size() > id);
        assert(entry_base_[id] != nullptr);
        return *entry_base_[id];
    }

    CacheableComponent &entry(const Component &comp) {
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
    // comp which is a part of the component stack
    inline CacheEntryID storeAsEntry(CacheableComponent &ccomp,
                                     CacheEntryID super_comp_id);

    // check quickly if the model count of the component is cached
    // if so, incorporate it into the model count of top
    // if not, store the packed version of it in the entry_base of the cache
    // bool manageNewComponent(StackLevel &top, CacheableComponent &packed_comp) {
    //   statistics_.num_cache_look_ups_++;
    //   unsigned table_ofs =  packed_comp.hashkey() & table_size_mask_;

    //   CacheEntryID act_id = table_[table_ofs];
    //   while(act_id){
    //     if (entry(act_id).equals(packed_comp)) {
    //       statistics_.incorporate_cache_hit(packed_comp);
    //       top.includeSolution(entry(act_id).model_count());
    //       return true;
    //     }
    //     act_id = entry(act_id).next_bucket_element();
    //   }
    //   return false;
    // }

    /**
     * Manager new component
     * @param top
     * @param packed_comp
     * @return The found cached component, nullptr otherwise.
     */
    CacheableComponent* manageNewComponent(StackLevel &top, CacheableComponent &packed_comp) {
        statistics_.num_cache_look_ups_++;
        unsigned table_ofs = packed_comp.hashkey() & table_size_mask_;
        CacheEntryID act_id = table_[table_ofs];
        if (config_.perform_pcc && !act_id) {
            return nullptr;
        }
        while (act_id) {
            if (entry(act_id).equals(packed_comp)) {
                statistics_.incorporate_cache_hit(packed_comp);
                entry(act_id).increase_cache_hit();
                top.includeSolution(entry(act_id).model_count());
                //std::cout << "found hit of mc: " << entry(act_id).model_count() << endl;
                return &entry(act_id);
            }
            act_id = entry(act_id).next_bucket_element();
        }
        return nullptr;
    }

    // unchecked erase of an entry from entry_base_
    void eraseEntry(CacheEntryID id) {
        statistics_.incorporate_cache_erase(*entry_base_[id]);
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

    vector<CacheableComponent*>& get_entry_base() {
        return entry_base_;
    }

private:

    void considerCacheResize() {
        if (entry_base_.size() > table_.size()) {
            reHashTable(2 * table_.size());
        }
    }

    void reHashTable(unsigned size) {

        table_.clear();
        table_.resize(size, 0);
        // we assert that table size is a power of 2
        // otherwise the table_size_mask_ doesn't work
        assert((table_.size() & (table_.size() - 1)) == 0);
        table_size_mask_ = table_.size() - 1;
        // if (!config_.quiet)
        //   cout << "ts " << table_.size() << " " << table_size_mask_ << endl;
        unsigned collisions = 0;
        for (unsigned id = 2; id < entry_base_.size(); id++)
            if (entry_base_[id] != nullptr) {
                entry_base_[id]->set_next_bucket_element(0);
                if (entry_base_[id]->modelCountFound()) {
                    unsigned table_ofs = tableEntry(id);
                    collisions += (table_[table_ofs] > 0 ? 1 : 0);
                    entry_base_[id]->set_next_bucket_element(table_[table_ofs]);
                    table_[table_ofs] = id;
                }
            }
        // if(!config_.quiet)
        //   cout << "coll " << collisions << endl;
    }

    unsigned tableEntry(CacheEntryID id) {
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

    unsigned table_size_mask_;

    DataAndStatistics &statistics_;
    SolverConfiguration &config_;

    unsigned long my_time_ = 0;
};


#include "component_cache-inl.h"


#endif /* COMPONENT_CACHE_H_ */
