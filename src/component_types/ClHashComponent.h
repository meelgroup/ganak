/*
 * ClHashComponent.h
 *
 *  Created on: May 5, 2020
 *      Author: VincentDerk
 *  Based on changes by smsharma1 to mthurley's component class
 */

#ifndef GANAK_CLHASHCOMPONENT_H
#define GANAK_CLHASHCOMPONENT_H

#include "cacheable_component.h"
#include "component.h"
#include "../clhash/clhash.h"
#include "difference_packed_component.h"
#include "sparsegraph_component.h"
#include <algorithm>

#include <math.h>

class ClHashComponent : public BaseComponent {

public:

    ClHashComponent() = default;

    ClHashComponent(uint64_t* clhashkey, unsigned cl_size, unsigned num_vars) : clhashkey_(clhashkey), cl_size_(cl_size), num_vars_(num_vars) {}

    inline ClHashComponent(vector<void *> &randomseedforCLHASH, DifferencePackedComponent &comp);

    inline ClHashComponent(vector<void *> &randomseedforCLHASH, SparseGraphCacheComponent &comp);

    ~ClHashComponent() override {
        delete[] clhashkey_;
    }

    unsigned num_variables() const final {
        return num_vars_;
    }

    unsigned data_only_byte_size() const override {
        return cl_size_ * sizeof(uint64_t);
    }

    unsigned long raw_data_byte_size() const override {
        return cl_size_ * sizeof(uint64_t) +
         model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t) + 
         variables_.size() * sizeof(VariableIndex) + sizeof(cache_hit);
    }

    // raw data size with the overhead
    // for the supposed 16byte alignment of malloc
    unsigned long sys_overhead_raw_data_byte_size() const override {
        unsigned ds = cl_size_ * sizeof(uint64_t);
        unsigned ms = model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
//      unsigned mask = 0xfffffff8;
//      return (ds & mask) + ((ds & 7)?8:0)
//            +(ms & mask) + ((ms & 7)?8:0);
        unsigned mask = 0xfffffff0;
        return (ds & mask) + ((ds & 15) ? 16 : 0)
               + (ms & mask) + ((ms & 15) ? 16 : 0);
    }

    /**
     * Check equality. Assumes comp is a ClHashComponent
     * @param comp The component to check equality with.
     * @return True iff this and comp are equivalent (same nb of variables, hashkey and clhash).
     */
    bool equals(const CacheableComponent &comp) const override {
        if (num_vars_ != comp.num_variables() || hashkey_ != comp.hashkey())
            return false;

        const auto& cast_comp = static_cast<const ClHashComponent &>(comp);
        if (cl_size_ != cast_comp.cl_size_)
            return false;

        uint64_t *clhashkey2 = cast_comp.clhashkey_;
        for (unsigned int i = 0; i < cl_size_; i++) {
            if (clhashkey2[i] != clhashkey_[i])
                return false;
        }

        return true;
    }

    void clear() override {
        // before deleting the contents of this component,
        // we should make sure that this component is not present in the component stack anymore!
        assert(isDeletable());
        delete[] clhashkey_;
        clhashkey_ = nullptr;
    }

private:
    uint64_t *clhashkey_ = nullptr;
    unsigned cl_size_ = 0; //TODO: Can be static because same for every component?
    unsigned num_vars_ = 0;
};

ClHashComponent::ClHashComponent(vector<void *> &random, DifferencePackedComponent &comp) {
    clhashkey_ = new uint64_t[random.capacity()];
    for (unsigned int i = 0; i < random.capacity(); i++) {
        clhasher h(random[i]);
        clhashkey_[i] = h((void *) comp.data(), sizeof(unsigned) * comp.data_size());
    }
    cl_size_ = random.capacity();
    num_vars_ = comp.num_variables();
    hashkey_ = comp.hashkey();
    variables_ = boost::container::flat_set<VariableIndex>(comp.variables_);
}

ClHashComponent::ClHashComponent(vector<void *> &random, SparseGraphCacheComponent &comp) {
    cl_size_ = random.capacity();
    num_vars_ = comp.num_variables();
    hashkey_ = comp.hashkey();
    variables_ = boost::container::flat_set<VariableIndex>(comp.variables_);

    // Construct clhashkey
    sparsegraph* graph = comp.get_graph();
    size_t data_size = graph->elen + graph->vlen;
    auto *data = new int[data_size];
    size_t data_offset = 0;
    for(size_t i = 0; i < graph->vlen; i++) {
        std::memcpy(data + data_offset, &graph->e[graph->v[i]], sizeof(int) * graph->d[i]);
        //sort
        std::sort(data + data_offset, data + data_offset + graph->d[i]);
        data_offset += graph->d[i];
        data[data_offset++] = -1; // sentinel
    }

    clhashkey_ = new uint64_t[random.capacity()];
    for (unsigned int i = 0; i < random.capacity(); i++) {
        clhasher h(random[i]);
        clhashkey_[i] = h(data, data_size);
    }
    delete[] data;
}

#endif //GANAK_CLHASHCOMPONENT_H
