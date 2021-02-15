/*
 * base_packed_component.h
 *
 *  Created on: Feb 5, 2013 as cacheable_component.h
 *      Author: mthurley
 *  Changed in 2020:
 *      Author: VincentDerk, smsharma1 and tvanbr
 */

#ifndef BASE_PACKED_COMPONENT_H_
#define BASE_PACKED_COMPONENT_H_

#include <assert.h>
#include <gmpxx.h>
#include <iostream>
#include "../primitive_types.h"
#include <boost/container/flat_set.hpp>

using namespace std;

class BaseComponent
{
public:
  BaseComponent() = default;

  BaseComponent(unsigned creation_time) : creation_time_(creation_time) {}

  virtual ~BaseComponent() = default;

  /**
     * Get the number of variables in this component.
     * @return The number of variables stored in this component
     */
  virtual unsigned num_variables() const = 0;

  virtual unsigned data_only_byte_size() const = 0;

  unsigned creation_time() const
  {
    return creation_time_;
  }

  const mpz_class &model_count() const
  {
    return model_count_;
  }

  unsigned alloc_of_model_count() const
  {
    return sizeof(mpz_class) + model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
  }

  /**
     * Set the creation time of this component
     * @param time The creation time
     */
  void set_creation_time(unsigned time)
  {
    creation_time_ = time;
  }

  /**
     * Set the model count of this component to rn
     * @param rn The modelcount to set
     * @param time The time it was solved at (?)
     */
  void set_model_count(const mpz_class &rn, unsigned time)
  {
    model_count_ = rn;
    length_solution_period_and_flags_ = (time - creation_time_) | (length_solution_period_and_flags_ & 1);
  }

  unsigned hashkey() const
  {
    return hashkey_;
  }

  /**
     * Whether the model count of this component is already found.
     * @return True iff the model count of this was set.
     */
  bool modelCountFound() const
  {
    return (length_solution_period_and_flags_ >> 1);
  }

  virtual bool equals(const BaseComponent &comp) const = 0;

  /**
     * Whether this component is set to be deletable.
     * @return True iff this can be deleted from the cache.
     */
  bool isDeletable() const
  {
    return length_solution_period_and_flags_ & 1;
  }

  /**
     * Set this to be deletable. This must only be called iff this is no longer connected to an active component
     * in the component stack.
     */
  void set_deletable()
  {
    length_solution_period_and_flags_ |= 1;
  }

  virtual void clear() = 0;
  // before deleting the contents of this component,
  // we should make sure that this component is not present in the component stack anymore!
  //        assert(isDeletable());
  //        delete[] data_;
  //        data_ = nullptr;
  //        delete[] clhashkey_;
  //        clhashkey_ = nullptr;
  //    }

  static unsigned _debug_static_val;

  // SIZES

  virtual unsigned long raw_data_byte_size() const = 0;

  virtual unsigned long sys_overhead_raw_data_byte_size() const = 0;

  /**
     * Size of this component (bytes).
     * @return
     */
  unsigned long SizeInBytes() const
  {
    // cout << "Size Divided "<< sizeof(GenericCacheableComponent<T_Component>) <<endl;
    // return sizeof(this);
    // return raw_data_byte_size();
    return sizeof(this)
          + raw_data_byte_size();
  }

  // the 48 = 16*3 in overhead stems from the three parts of the component
  // being dynamically allocated (i.e. the GenericCacheableComponent itself,
  // the data_ and the model_count data
  unsigned long sys_overhead_SizeInBytes() const
  { //TODO: Verify
    return sys_overhead_raw_data_byte_size() + 48 + sizeof(this);
    //return sizeof(GenericCacheableComponent<T_Component>)
    //       + T_Component::sys_overhead_raw_data_byte_size()
    //       +48;
  }

  // BEGIN Cache Pollution Management
  void set_father(CacheEntryID f)
  {
    father_ = f;
  }

  CacheEntryID father() const
  {
    return father_;
  }

  void set_next_sibling(CacheEntryID sibling)
  {
    next_sibling_ = sibling;
  }

  CacheEntryID next_sibling()
  {
    return next_sibling_;
  }

  void set_first_descendant(CacheEntryID descendant)
  {
    first_descendant_ = descendant;
  }

  CacheEntryID first_descendant()
  {
    return first_descendant_;
  }

  void set_next_bucket_element(CacheEntryID entry)
  {
    next_bucket_element_ = entry;
  }

  CacheEntryID next_bucket_element()
  {
    return next_bucket_element_;
  }

  inline unsigned get_cache_hit_count() const
  {
    return cache_hit;
  }

  inline void increase_cache_hit()
  {
    cache_hit++;
  }

  /**
     * Add all variables in the given component to the set of variables in this component.
     * @param comp The component whose variables have to be copied.
     */
  void add_variables_of(BaseComponent &comp)
  {
    for (auto v : comp.variables_)
    {
      variables_.insert(v);
    }
  }

  /**
     * Variables that formed a component isomorphic to this. When config_.use_icsvsads is true,
     * this should contain the variables of this component (the variables when it was stored + the variables of
     * components isomorphic to this, added during cache hits.).
     */
  boost::container::flat_set<VariableIndex> variables_;

protected:
  unsigned hashkey_ = 0;

  unsigned cache_hit = 0;

  mpz_class model_count_;

  unsigned creation_time_ = 1;

  /**
     * Length solution period = length_solution_period_and_flags_ >> 1
     * length_solution_period == 0 means unsolved and the first bit is "delete permitted"
     * Deletion is permitted only after the copy of this component in the stack does not exist anymore.
     */
  unsigned length_solution_period_and_flags_ = 0;

private:
  CacheEntryID next_bucket_element_ = 0;

  // theFather and theDescendants:
  // each CCacheEntry is a Node in a tree which represents the relationship
  // of the components stored
  CacheEntryID father_ = 0;
  CacheEntryID first_descendant_ = 0;
  CacheEntryID next_sibling_ = 0;
};

typedef BaseComponent CacheableComponent;

#endif /* BASE_PACKED_COMPONENT_H_ */
