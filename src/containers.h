/******************************************
Copyright (C) 2023 Marc Thurley

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

  inline typename vector<_T>::const_iterator begin() const {
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

  Lit end_lit() const {
    return Lit(size() / 2, false);
  }

  using vector<_T>::end;
  using vector<_T>::size;
  using vector<_T>::clear;
  using vector<_T>::push_back;
};
