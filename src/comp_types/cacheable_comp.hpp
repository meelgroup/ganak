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
#include <cstdint>
#include <gmp.h>
#include <arjun/arjun.h>
#include "comp.hpp"
#include "hashval.hpp"

namespace GanakInt {

// Non-negative integer counts of at most 128 bits, stored in the entry itself
// instead of as a heap FMpz (two allocations, 64+ bytes).
struct SmallIntCodec {
#if GMP_NUMB_BITS == 64 && GMP_NAIL_BITS == 0
  static constexpr bool supported = true;
#else
  static constexpr bool supported = false;
#endif
  static bool pack(const CMSat::Field& f, uint64_t out[2]) {
    if constexpr (!supported) return false;
    const auto* z = static_cast<const ArjunNS::FMpz&>(f).val.get_mpz_t();
    if (z->_mp_size < 0 || z->_mp_size > 2) return false;
    out[0] = z->_mp_size > 0 ? z->_mp_d[0] : 0;
    out[1] = z->_mp_size > 1 ? z->_mp_d[1] : 0;
    return true;
  }
  static void unpack(const uint64_t in[2], CMSat::Field& f) {
    if constexpr (supported) {
      auto* z = static_cast<ArjunNS::FMpz&>(f).val.get_mpz_t();
      mp_limb_t* p = mpz_limbs_write(z, 2);
      p[0] = in[0];
      p[1] = in[1];
      mpz_limbs_finish(z, in[1] ? 2 : (in[0] ? 1 : 0));
    }
  }
};

// Trivially copyable: the cache moves entries with plain copies when it
// compacts, and frees their heap parts explicitly via release().
template<typename T>
class CacheableComp: private T {
public:
  enum : uint8_t { cnt_none = 0, cnt_zero = 1, cnt_inline = 2, cnt_heap = 3 };

  CacheableComp() = default;
  void set_comp(const Comp &comp, const uint32_t hash_seed, const BPCSizes& bpc) {
    hashkey = T::set_comp(comp, hash_seed, bpc);
  }
  // the packed data now belongs to a copy of this, forget it
  void forget() { T::forget(); }

  [[nodiscard]] uint64_t last_used_time() const { return last_used_time_; }
  void set_last_used_time(uint64_t time) { last_used_time_ = time; }
  void avg_last_used_time(uint64_t time, uint64_t div) {
    assert(time >= last_used_time_);
    last_used_time_ += (time-last_used_time_)/div;
  }
  [[nodiscard]] bool is_deletable() const { return delete_permitted; }
  void set_deletable() { delete_permitted = true; }
  void reset_meta(uint64_t time) {
    last_used_time_ = time;
    delete_permitted = false;
    cnt_state = cnt_none;
  }

  [[nodiscard]] bool model_count_found() const { return cnt_state != cnt_none; }
  void set_model_count(const FF& rn, const bool int_inline) {
    assert(cnt_state == cnt_none);
    if (rn == nullptr) return;
    // zero is a valid model count (UNSAT component), cache it too
    if (rn->is_zero()) { cnt_state = cnt_zero; return; }
    if (int_inline && SmallIntCodec::pack(*rn, cnt.small)) { cnt_state = cnt_inline; return; }
    cnt.heap = rn->dup().release();
    cnt_state = cnt_heap;
  }
  [[nodiscard]] bool count_is_inline() const { return cnt_state == cnt_inline; }
  // nullptr means zero. scratch must be an FMpz, used for inline counts
  const CMSat::Field* model_count(CMSat::Field& scratch) const {
    switch (cnt_state) {
      case cnt_inline: SmallIntCodec::unpack(cnt.small, scratch); return &scratch;
      case cnt_heap: return cnt.heap;
      default: assert(cnt_state == cnt_zero); return nullptr;
    }
  }

  void release() {
    if (cnt_state == cnt_heap) delete cnt.heap;
    cnt_state = cnt_none;
    T::release();
  }

  uint64_t get_hashkey() const  { return hashkey.hash; }
  const HashVal& get_full_hashkey() const  { return hashkey; }
  [[nodiscard]] uint64_t extra_bytes() const {
    return (cnt_state == cnt_heap ? cnt.heap->bytes_used() : 0) + T::comp_bytes();
  }
  // hashkey matches, but the content (if stored) differs
  [[nodiscard]] bool content_differs(const CacheableComp<T>& other) const {
    return !T::equals(other);
  }
  bool equals(const CacheableComp<T>& other) const {
    // full 128-bit key compare (collision needs BOTH words to match, ~2^-128) plus,
    // for the content-comparing variant, the exact-content check in T::equals()
    return hashkey == other.hashkey && T::equals(other);
  }

private:
  HashVal hashkey;
  union {
    CMSat::Field* heap;
    uint64_t small[2];
  } cnt = {nullptr};
  uint64_t last_used_time_:61 = 1; //effectively the score
  uint64_t delete_permitted:1 = false;
  uint64_t cnt_state:2 = cnt_none;
};

}
