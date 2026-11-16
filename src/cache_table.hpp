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

#include <bit>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <new>
#include "hashval.hpp"
#ifdef __SSE2__
#include <emmintrin.h>
#endif

namespace GanakInt {

// Open-addressing hash table of 32-bit ids, 12 per 64-byte bucket, each with
// an 8-bit tag from hash2, so a probe is one cache line and a miss rarely reads
// an entry. Linear probing over buckets; a lookup stops at the first bucket
// with an empty slot.
//
// Deletion: a slot is emptied only if its bucket still has an empty slot.
// Slots go empty -> live -> empty/tombstone, and a full bucket never gets an
// empty slot back, so such a bucket was never full and no probe went past it.
// Otherwise it becomes a tombstone, reused by inserts, dropped on rebuild.
class CacheTable {
  struct alignas(64) Bucket {
    uint8_t tag[16]; // only 12 used, the rest stay 0 and are masked off
    uint32_t id[12];
  };
  static_assert(sizeof(Bucket) == 64);
  static constexpr uint32_t slots = 12;
  static constexpr uint32_t slot_mask = (1U << slots) - 1;
  static constexpr uint8_t empty = 0;
  static constexpr uint8_t tomb = 1;

public:
  CacheTable() = default;
  ~CacheTable() { dealloc(); }
  CacheTable(const CacheTable&) = delete;
  CacheTable& operator=(const CacheTable&) = delete;

  // Empties the table and sizes it so that n items fill it to `load`
  void reset(uint64_t n, double load) {
    uint64_t want = (uint64_t)((double)n / (load * slots)) + 1;
    want = std::max<uint64_t>(want, 16);
    if (want != nb) {
      dealloc();
      nb = want;
      bk = static_cast<Bucket*>(::operator new[](nb * sizeof(Bucket), std::align_val_t(64)));
    }
    std::memset(static_cast<void*>(bk), 0, nb * sizeof(Bucket));
    used = 0;
    live = 0;
  }

  [[nodiscard]] bool over_load(double load) const { return (double)(used + 1) > load * (double)(nb * slots); }
  [[nodiscard]] uint64_t num_live() const { return live; }
  [[nodiscard]] uint64_t bytes() const { return nb * sizeof(Bucket); }
  // bytes after a reset(n, load) -- for the memory budget
  [[nodiscard]] static uint64_t bytes_for(uint64_t n, double load) {
    return std::max<uint64_t>((uint64_t)((double)n / (load * slots)) + 1, 16) * sizeof(Bucket);
  }

  void prefetch(const HashVal& h) const { __builtin_prefetch(&bk[home(h)], 1); }

  // eq(id) is called on tag matches only
  template<typename EQ>
  uint32_t find(const HashVal& h, EQ&& eq) const {
    const uint8_t t = tag_of(h);
    for (uint64_t b = home(h);; b = next(b)) {
      const Bucket& x = bk[b];
      uint32_t m = match(x, t);
      while (m) {
        const uint32_t id = x.id[std::countr_zero(m)];
        if (eq(id)) return id;
        m &= m - 1;
      }
      if (match(x, empty)) return 0;
    }
  }

  void insert(const HashVal& h, const uint32_t id) {
    assert(id != 0);
    const uint8_t t = tag_of(h);
    for (uint64_t b = home(h);; b = next(b)) {
      Bucket& x = bk[b];
      const uint32_t m = match(x, empty) | match(x, tomb);
      if (!m) continue;
      const uint32_t k = std::countr_zero(m);
      if (x.tag[k] == empty) used++;
      x.tag[k] = t;
      x.id[k] = id;
      live++;
      return;
    }
  }

  void erase(const HashVal& h, const uint32_t id) {
    const uint8_t t = tag_of(h);
    for (uint64_t b = home(h);; b = next(b)) {
      Bucket& x = bk[b];
      uint32_t m = match(x, t);
      while (m) {
        const uint32_t k = std::countr_zero(m);
        if (x.id[k] == id) {
          if (match(x, empty)) { x.tag[k] = empty; used--; }
          else x.tag[k] = tomb;
          live--;
          return;
        }
        m &= m - 1;
      }
      assert(!match(x, empty) && "erasing an id not in the table");
    }
  }

private:
  Bucket* bk = nullptr;
  uint64_t nb = 0;   // number of buckets, need not be a power of 2
  uint64_t used = 0; // live + tombstones
  uint64_t live = 0;

  void dealloc() {
    if (bk) ::operator delete[](bk, std::align_val_t(64));
    bk = nullptr;
    nb = 0;
  }
  static uint8_t tag_of(const HashVal& h) {
    const uint8_t t = (uint8_t)h.hash2;
    return t < 2 ? t + 2 : t;
  }
  // multiply-shift range reduction, uses the top 32 bits of hash
  [[nodiscard]] uint64_t home(const HashVal& h) const { return ((h.hash >> 32) * nb) >> 32; }
  [[nodiscard]] uint64_t next(uint64_t b) const { return b + 1 == nb ? 0 : b + 1; }
  static uint32_t match(const Bucket& x, const uint8_t t) {
#ifdef __SSE2__
    const __m128i v = _mm_load_si128(reinterpret_cast<const __m128i*>(x.tag));
    return (uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(v, _mm_set1_epi8((char)t))) & slot_mask;
#else
    uint32_t m = 0;
    for (uint32_t k = 0; k < slots; k++) m |= (uint32_t)(x.tag[k] == t) << k;
    return m;
#endif
  }
};

}
