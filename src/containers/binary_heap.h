/*
 * binary_heap.h
 *
 *  Created on: May 8, 2013
 *      Author: mthurley
 */

#ifndef BINARY_HEAP_H_
#define BINARY_HEAP_H_

#include <assert.h>
#include <vector>


template <class Comparator>
class BinaryHeap {
  // the values corresponding to the usigned keys
  std::vector<unsigned> values_;
  // keys to values, note that a zero key means: not in HEAP
  std::vector<unsigned> keys_;
  // we have a max heap w.r.t. this comparator
  Comparator less_;
public:
  BinaryHeap(const Comparator & comp)
: less_(comp){
  }

  inline void init(unsigned num_values);

  inline unsigned deleteMax();


  void incorporateIncrease(unsigned value) {//*
    assert(keys_[value]);
    swim(keys_[value]);
  }

  bool contains(unsigned value){
    return value < keys_.size() && keys_[value] > 0;
  }

  bool empty(){
    assert(!values_.empty());
    return values_.size() == 1;
  }

  inline void insert(unsigned value);

private:

  bool less(unsigned keyA, unsigned keyB){
    return less_(values_[keyA],values_[keyB]);
  }

  unsigned parentKey(unsigned key){
    return key>>1;
  }

  inline void swim(unsigned key);

  inline void sink(unsigned key);



  /// Test Functions
  inline void test_checkRep();

};


template <class Comparator>
void BinaryHeap<Comparator>::init(unsigned num_values){
  keys_.clear();
  keys_.resize(num_values+1,0);
  values_.clear();
  values_.reserve(num_values+1);
  values_.push_back(0);
  // TODO: add correct initialization
  for (unsigned n = 1; n <= num_values; n++) {
    insert(n);
  }
}

template <class Comparator>
unsigned BinaryHeap<Comparator>::deleteMax(){ //*
  assert(values_.size() >=2);
  unsigned value = values_[1];
  values_[1] = values_.back();
  keys_[values_.back()] = 1;
  keys_[value] = 0;
  values_.pop_back();
  if(!empty())
    sink(1);
  //test_checkRep();
  return value;
}

template <class Comparator>
void BinaryHeap<Comparator>::insert(unsigned value){//*
  assert(!contains(value));
  values_.push_back(value);
  keys_[value] = values_.size() - 1;
  swim(keys_[value]);
  sink(keys_[value]);

  //test_checkRep();
}

template <class Comparator>
void BinaryHeap<Comparator>::swim(unsigned key){
  assert(key < values_.size());
  unsigned value = values_[key];

  unsigned parent_key = parentKey(key);
  while(key > 1 && less_(values_[parent_key],value)){
    keys_[values_[parent_key]] = key;
    values_[key] = values_[parent_key];
    key = parent_key;
    parent_key = parentKey(parent_key);
  }
  values_[key] = value;
  keys_[value] = key;
}

template <class Comparator>
void BinaryHeap<Comparator>::sink(unsigned key){
  assert(key < values_.size());
  unsigned value = values_[key];
  while((key<<1) < values_.size()){
    unsigned j = key<<1;
    if(j+1 < values_.size() && less(j,j+1)) ++j;
    if(!less_(value,values_[j])) break;
    values_[key] = values_[j];
    keys_[values_[j]] = key;
    key=j;
  }
  values_[key] = value;
  keys_[value] = key;
}


////////////////////////////////////////////
//// Testing

template <class Comparator>
void BinaryHeap<Comparator>::test_checkRep(){
  assert(!values_.empty());
  assert(values_[0] == 0);
  for(unsigned k = 1; k < values_.size(); k++){
    unsigned left = k << 1;
    unsigned right = (k << 1) + 1;
    if(left < values_.size()) assert(!less_(values_[k], values_[left]));
    if(right < values_.size()) assert(!less_(values_[k],values_[right]));
  }
}

#endif /* BINARY_HEAP_H_ */
