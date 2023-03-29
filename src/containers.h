/*
 * containers.h
 *
 *  Created on: Jun 27, 2012
 *      Author: Marc Thurley
 */

#pragma once

#include "structures.h"

template<class _T>
class LiteralIndexedVector: protected vector<_T> {

public:
	LiteralIndexedVector(const uint32_t sz = 0) :
			vector<_T>(sz * 2) {
	}
	LiteralIndexedVector(const uint32_t sz,
			const typename vector<_T>::value_type& __value) :
			vector<_T>(sz * 2, __value) {
	}
	inline _T &operator[](const Lit lit) {
		return *(vector<_T>::begin() + lit.raw());
	}

	inline const _T &operator[](const Lit &lit) const {
		return *(vector<_T>::begin() + lit.raw());
	}

	inline typename vector<_T>::iterator begin() {
		return vector<_T>::begin() + 2;
	}

	bool empty() const {
		return vector<_T>::empty();
	}

	void resize(uint32_t _size) {
		vector<_T>::resize(_size * 2);
	}
	void resize(uint32_t _size, const typename vector<_T>::value_type& _value) {
		vector<_T>::resize(_size * 2, _value);
	}

	void reserve(uint32_t _size) {
		vector<_T>::reserve(_size * 2);
	}

	Lit end_lit() {
		return Lit(size() / 2, false);
	}

	using vector<_T>::end;
	using vector<_T>::size;
	using vector<_T>::clear;
	using vector<_T>::push_back;
};
