/*
 * containers.h
 *
 *  Created on: Jun 27, 2012
 *      Author: Marc Thurley
 */

#ifndef CONTAINERS_H_
#define CONTAINERS_H_

#include "structures.h"

template<class _T>
class LiteralIndexedVector: protected vector<_T> {

public:
	LiteralIndexedVector(unsigned size = 0) :
			vector<_T>(size * 2) {
	}
	LiteralIndexedVector(unsigned size,
			const typename vector<_T>::value_type& __value) :
			vector<_T>(size * 2, __value) {
	}
	inline _T &operator[](const LiteralID lit) {
		return *(vector<_T>::begin() + lit.raw());
	}

	inline const _T &operator[](const LiteralID &lit) const {
		return *(vector<_T>::begin() + lit.raw());
	}

	inline typename vector<_T>::iterator begin() {
		return vector<_T>::begin() + 2;
	}

	void resize(unsigned _size) {
		vector<_T>::resize(_size * 2);
	}
	void resize(unsigned _size, const typename vector<_T>::value_type& _value) {
		vector<_T>::resize(_size * 2, _value);
	}

	void reserve(unsigned _size) {
		vector<_T>::reserve(_size * 2);
	}

	LiteralID end_lit() {
		return LiteralID(size() / 2, false);
	}

	using vector<_T>::end;
	using vector<_T>::size;
	using vector<_T>::clear;
	using vector<_T>::push_back;
};

#endif /* CONTAINERS_H_ */
