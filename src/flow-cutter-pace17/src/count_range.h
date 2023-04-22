#ifndef COUNT_RANGE_H
#define COUNT_RANGE_H

#include "range.h"
#include <cassert>
#include <iterator>

struct CountIterator{
	typedef int value_type;
	typedef int difference_type;
	typedef const int* pointer;
	typedef const int& reference;
	typedef std::random_access_iterator_tag iterator_category;

	CountIterator&operator++(){ ++n_; return *this;}
	CountIterator operator++(int) {CountIterator tmp(*this); operator++(); return tmp;}
	CountIterator&operator--(){ --n_; return *this;}
	CountIterator operator--(int) {CountIterator tmp(*this); operator++(); return tmp;}
	int operator*() const {return n_;}

	const int*operator->() const {return &n_;}

	int operator[](int o)const{return n_ + o;}
	CountIterator&operator+=(CountIterator::difference_type o){n_ += o; return *this;}
	CountIterator&operator-=(CountIterator::difference_type o){n_ -= o; return *this;}

	int n_;
};

inline bool operator==(CountIterator l, CountIterator r){return l.n_ == r.n_;}
inline bool operator!=(CountIterator l, CountIterator r){return l.n_ != r.n_;}
inline bool operator< (CountIterator l, CountIterator r){return l.n_ <  r.n_;}
inline bool operator> (CountIterator l, CountIterator r){return l.n_ >  r.n_;}
inline bool operator<=(CountIterator l, CountIterator r){return l.n_ <= r.n_;}
inline bool operator>=(CountIterator l, CountIterator r){return l.n_ >= r.n_;}

inline CountIterator::difference_type operator-(CountIterator l, CountIterator r){return l.n_ - r.n_;}
inline CountIterator operator-(CountIterator l, CountIterator::difference_type r){return {l.n_ - r};}
inline CountIterator operator+(CountIterator l, CountIterator::difference_type r){return {l.n_ + r};}
inline CountIterator operator+(CountIterator::difference_type l, CountIterator r){return {l + r.n_};}

typedef Range<CountIterator> CountRange;

inline CountRange count_range(int n){assert(n >= 0); return {CountIterator{0}, CountIterator{n}}; }
inline CountRange count_range(int begin, int end){assert(begin <= end);return {CountIterator{begin}, CountIterator{end}};}

#endif
