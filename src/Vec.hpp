/***************************************************************************************
Copyright (c) 2003-2007, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***************************************************************************************/

#pragma once

#include <cassert>
#include <new>
#include <cstdint>
#include <cstdlib>
#include <cerrno>
#include <limits>
#include <vector>
#include <utility>

using std::numeric_limits;

//=================================================================================================
// Automatically resizable arrays
//
// NOTE! Don't use this vector on datatypes that cannot be re-located in memory (with realloc)

template<class T>
class vec {
public:
    T*  dat;
    T* begin() { return dat; }
    T* end() { return dat + sz; }
    const T* begin() const { return dat; }
    const T* end() const { return dat + sz; }
private:
    uint32_t sz;
    uint32_t cap;

    // Don't allow copying (error prone):
    vec<T>& operator= (vec<T>& ) = delete;
    vec (const vec<T>& other) = delete;

    // Helpers for calculating next capacity:
    static inline uint32_t  imax   (int32_t x, int32_t y)
    {
        int32_t mask = (y - x) >> (sizeof(uint32_t) * 8 - 1);
        return (x & mask) + (y & (~mask));
    }

public:
    vec() : dat(nullptr) , sz(0)   , cap(0)    { }
    explicit vec(uint32_t size) : dat(nullptr), sz(0), cap(0) {
        growTo(size);
    }
    vec(uint32_t size, const T& pad) : dat(nullptr), sz(0), cap(0) {
        growTo(size, pad);
    }
    ~vec() { clear(true); }

    // allow moving, it's easy
    vec (vec<T>&& other) noexcept {
      sz = other.sz;
      cap = other.cap;
      dat = other.dat;
      other.sz = 0;
      other.cap = 0;
      other.dat = nullptr;
    }
    uint32_t size() const { return sz; }
    void shrink   (uint32_t nelems) {
        assert(nelems <= sz);
        for (uint32_t i = 0; i < nelems; i++) {
            sz--, dat[sz].~T();
        }
    }
    void shrink_  (uint32_t nelems) {
        assert(nelems <= sz);
        sz -= nelems;
    }
    uint32_t capacity () const { return cap; }
    void capacity (int32_t min_cap);
    void growTo   (uint32_t size);
    void growTo   (uint32_t size, const T& pad);
    void clear    (bool dealloc = false);
    auto data() const { return dat; }
    auto data() { return dat; }
    vec<T>& operator= (const std::vector<T>& v) {
      clear();
      for (const auto& e : v) push_back(e);
      return *this;
    }
    void push_back(const T& elem) {
        if (sz == cap) capacity(sz + 1);
        dat[sz++] = elem;
    }
    void pop_back() {
        assert(sz > 0);
        sz--, dat[sz].~T();
    }
    // NOTE: it seems possible that overflow can happen in the 'sz+1' expression of 'push()', but
    // in fact it can not since it requires that 'cap' is equal to INT_MAX. This in turn can not
    // happen given the way capacities are calculated (below). Essentially, all capacities are
    // even, but INT_MAX is odd.

    const T& back  () const { return dat[sz - 1]; }
    T&       back  () { return dat[sz - 1]; }
    const T& operator [] (uint32_t index) const { return dat[index]; }
    T&       operator [] (uint32_t index) { return dat[index]; }

    // Duplicatation (preferred instead):
    void copyTo(vec<T>& copy) const {
        copy.clear();
        copy.growTo(sz);
        for (uint32_t i = 0; i < sz; i++) copy[i] = dat[i];
    }
    void moveTo(vec<T>& dest) {
        dest.clear(true);
        dest.dat = dat;
        dest.sz = sz;
        dest.cap = cap;
        dat = nullptr;
        sz = 0;
        cap = 0;
    }
    void swap(vec<T>& dest) {
        std::swap(dest.dat, dat);
        std::swap(dest.sz, sz);
        std::swap(dest.cap, cap);
    }

    void resize(uint32_t s) {
        if (s < sz) shrink(sz - s);
        else growTo(s);
    }
    void insert(uint32_t num) { growTo(sz+num); }
    bool empty() const { return sz == 0; }
    void shrink_to_fit() {
        if (sz == 0) {
          clear(true);
          return;
        }
        for(uint32_t i = sz; i < cap; i++) dat[i].~T();

        T* data2 = (T*)realloc(dat, sz*sizeof(T));
        if (data2 == 0) throw std::bad_alloc();
        dat = data2;
        cap = sz;
     }
};


// Fixes by @Topologist from GitHub. Thank you so much!
template<class T>
void vec<T>::capacity(int32_t min_cap)
{
    if ((int32_t)cap >= min_cap) return;

    // NOTE: grow by approximately 3/2
    uint32_t add = imax((min_cap - (int32_t)cap + 1) & ~1, (((int32_t)cap >> 1) + 2) & ~1);
    if (add > numeric_limits<uint32_t>::max() - cap) throw std::bad_alloc();
    cap += add;

    // This avoids memory fragmentation by many reallocations
    int32_t new_size = 2;
    while (new_size < min_cap) new_size *= 2;
    if ((new_size * 2 / 3) > min_cap) new_size = new_size * 2 / 3;
    cap = new_size;

    if (((dat = (T*)::realloc(dat, cap * sizeof(T))) == nullptr) && errno == ENOMEM)
        throw std::bad_alloc();
    for(uint32_t i = sz; i < cap; i++) new (&dat[i]) T();
}

template<class T>
void vec<T>::growTo(uint32_t size, const T& pad)
{
    if (sz >= size) return;
    capacity(size);
    for (uint32_t i = sz; i < size; i++) new (&dat[i]) T(pad);
    sz = size;
}

template<class T>
void vec<T>::growTo(uint32_t size)
{
    if (sz >= size) return;
    capacity(size);
    for (uint32_t i = sz; i < size; i++) new (&dat[i]) T();
    sz = size;
}

template<class T>
void vec<T>::clear(bool dealloc)
{
    if (dat != nullptr) {
        sz = 0;
        if (dealloc) {
          for(uint32_t i = 0; i < cap; i++) dat[i].~T();
          free(dat);
          dat = nullptr;
          cap = 0;
        }
    }
}
