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

#include <cassert>
#include <vector>

#include "../primitive_types.hpp"
#include "hashed_comp.hpp"

class Comp;
class CompArchetype;

// Adds Structure to PackedComp that is
// necessary to store it in the cache
// namely, the descendant tree structure that
// allows for the removal of cache pollutions
class CacheableComp: public HashedComp {
public:
  CacheableComp() = default;
  CacheableComp(void* hash_seed, Comp &comp) : HashedComp(hash_seed, comp) { }

  uint32_t bignum_bytes() const { return HashedComp::bignum_bytes(); }

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
    if (father_ == std::numeric_limits<uint32_t>::max()) assert (model_count_ == nullptr);
    return father_ == std::numeric_limits<uint32_t>::max();
  }
  void set_free() {
    father_ = std::numeric_limits<uint32_t>::max();
    delete model_count_;
    model_count_ = nullptr;
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
