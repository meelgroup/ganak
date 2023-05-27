/*
 * cacheable_comp.h
 *
 *  Created on: Feb 21, 2013
 *      Author: mthurley
 */

#ifndef CACHEABLE_COMPONENT_H_
#define CACHEABLE_COMPONENT_H_

#include <assert.h>
#include <vector>

#include "primitive_types.h"
#include "difference_packed_comp.h"

class Component;
class ComponentArchetype;

// GenericCacheableComponent Adds Structure to PackedComponent that is
// necessary to store it in the cache
// namely, the descendant tree structure that
// allows for the removal of cache pollutions

class CacheableComponent: public DifferencePackedComponent {
public:
  CacheableComponent() { }
  CacheableComponent(void* randomseedforCLHASH, Component &comp,
      const BPCSizes& sz, uint32_t* tmp_data) :
      DifferencePackedComponent(randomseedforCLHASH, comp, sz, tmp_data) {
  }

  uint32_t SizeInBytes() const {
    return DifferencePackedComponent::raw_data_byte_size();
  }

  // Cache Pollution Management
  void set_father(CacheEntryID f) { father_ = f; }
  CacheEntryID father() const { return father_; }
  void set_next_sibling(CacheEntryID sibling) { next_sibling_ = sibling; }
  CacheEntryID next_sibling() const { return next_sibling_; }
  void set_first_descendant(CacheEntryID descendant) { first_descendant_ = descendant; }
  CacheEntryID first_descendant() const { return first_descendant_; }
  void set_next_bucket_element(CacheEntryID entry) { next_bucket_element_ = entry; }
  CacheEntryID next_bucket_element() const { return next_bucket_element_; }
  bool is_free() const {
    return father_ == std::numeric_limits<uint32_t>::max();
  }
  void set_free() {
    father_ = std::numeric_limits<uint32_t>::max();
    delete model_count_;
    model_count_ = NULL;
  }

private:

  CacheEntryID next_bucket_element_ = 0;

  // father and descendants:
  // each CCacheEntry is a Node in a tree which represents the relationship
  // of the comps stored
  CacheEntryID father_ = 0;
  CacheEntryID first_descendant_ = 0;
  CacheEntryID next_sibling_ = 0;
};

#endif /* CACHEABLE_COMPONENT_H_ */
