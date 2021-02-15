/*
 * sparsegraph_component.h
 *
 *  Created in 2020
 *      Author: VincentDerk
 *      Author: smsharma1
 *      Author: tvanbr
 */

#ifndef GANAK_SPARSEGRAPH_COMPONENT_H
#define GANAK_SPARSEGRAPH_COMPONENT_H

#include "cacheable_component.h"
#include "../nauty/nauty.h"
#include "../nauty/naugroup.h"
#include "../nauty/nautinv.h"
#include "../nauty/naututil.h"

class SparseGraphCacheComponent : public BaseComponent {
public:

    SparseGraphCacheComponent(sparsegraph *cg, long hashkey, unsigned num_of_variables) {
        canonical_graph = cg;
        hashkey_ = hashkey;
        num_of_variables_ = num_of_variables;
    }

    ~SparseGraphCacheComponent() override {
        SG_FREE(*canonical_graph);
        free(canonical_graph);
    }

    unsigned num_variables() const final {
        return num_of_variables_;
    }

    unsigned data_only_byte_size() const override {
        //- size_t nde;  /* Number of directed edges (loops contribute only 1) */
        //- size_t *v;   /* Array of indexes into e[*] */
        //- int nv;      /* Number of vertices */
        //- int *d;      /* Array with out-degree of each vertex */
        //- int *e;      /* Array to hold lists of neighbours */
        //sg_weight *w;      /* Not implemented, should be NULL. */
        //size_t vlen,dlen,elen,wlen;  /* Sizes of arrays in units of type */
        return sizeof(int) * (canonical_graph->dlen + canonical_graph->elen + 1) +
               sizeof(size_t) * (canonical_graph->vlen + 5) +
               sizeof(size_t *) + 2 * sizeof(int *) + sizeof(sg_weight *);
    }

    unsigned long raw_data_byte_size() const override {
        return sizeof(*this) +
               model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t) +
               data_only_byte_size()
               + variables_.size() * sizeof(VariableIndex) + sizeof(cache_hit);
    }

    // raw data size with the overhead
    // for the supposed 16byte alignment of malloc
    unsigned long sys_overhead_raw_data_byte_size() const override {
        //TODO: Still correct?
        //unsigned ds = data_size();
        //unsigned ms = model_count().get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
//      unsigned mask = 0xfffffff8;
//      return (ds & mask) + ((ds & 7)?8:0) +(ms & mask) + ((ms & 7)?8:0);
        //unsigned mask = 0xfffffff0;
        //return (ds & mask) + ((ds & 15) ? 16 : 0)
        //       + (ms & mask) + ((ms & 15) ? 16 : 0);
        return raw_data_byte_size();
    }

    bool equals(const CacheableComponent &comp) const final {
        if (num_of_variables_ != comp.num_variables() || hashkey_ != comp.hashkey())
            return false;

        // DEBUG INFO - are there hashkey collisions and how approximate would it be to only use those and not cg?
//            bool same = aresame_sg(canonical_graph, static_cast<const SparseGraphCacheComponent &>(comp).canonical_graph);
//            if (!same) {
//                std::cout << "Collision\n";
//            }
//            return same;
        return aresame_sg(canonical_graph, static_cast<const SparseGraphCacheComponent &>(comp).canonical_graph);
    }

    void clear() override {
        SG_FREE(*canonical_graph);
        free(canonical_graph);
    }

    sparsegraph* get_graph() const {
        return canonical_graph;
    }

private:
    sparsegraph *canonical_graph;
    unsigned num_of_variables_;
};

#endif //GANAK_SPARSEGRAPH_COMPONENT_H
