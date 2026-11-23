/******************************************
Copyright (C) 2026 Authors of GANAK, see AUTHORS file

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
// Exact bounded component algorithms. Local variable i is bit i in a clause.
#pragma once
#include <algorithm>
#include <array>
#include <bit>
#if defined(__SSE4_2__) && defined(__x86_64__)
#include <nmmintrin.h>
#endif
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#include <tuple>
#include <type_traits>
#include <vector>

namespace GanakInt::frontier {
using Mask = uint64_t;
struct Clause {
  Mask pos = 0, neg = 0;
  Mask vars() const { return pos | neg; }
};
inline Mask bit(unsigned i) { return Mask{1} << i; }

// Clause-index masks for components with more than 64 clauses.
template <unsigned Words> struct WideMask {
  std::array<uint64_t, Words> w{};
  WideMask() = default;
  WideMask(uint64_t low) { w[0] = low; }  // NOLINT: zero and small literals
  WideMask& operator|=(const WideMask& o) { for (unsigned i = 0; i < Words; ++i) w[i] |= o.w[i]; return *this; }
  WideMask& operator&=(const WideMask& o) { for (unsigned i = 0; i < Words; ++i) w[i] &= o.w[i]; return *this; }
  friend WideMask operator|(WideMask a, const WideMask& b) { return a |= b; }
  friend WideMask operator&(WideMask a, const WideMask& b) { return a &= b; }
  friend WideMask operator~(WideMask a) { for (auto& x : a.w) x = ~x; return a; }
  bool operator==(const WideMask&) const = default;
  explicit operator bool() const {
    for (auto x : w) if (x) return true;
    return false;
  }
};
inline int popcount(uint64_t m) { return std::popcount(m); }
template <unsigned Words> int popcount(const WideMask<Words>& m) {
  int total = 0;
  for (auto x : m.w) total += std::popcount(x);
  return total;
}
template <class M> M clause_bit(unsigned j) {
  if constexpr (std::is_same_v<M, uint64_t>) return Mask{1} << j;
  else { M m; m.w[j / 64] = Mask{1} << (j % 64); return m; }
}
inline std::size_t hash_key(uint64_t key) {
  // Clause masks often differ only in their high bits.
#if defined(__SSE4_2__) && defined(__x86_64__)
  return _mm_crc32_u64(0, key);
#else
  Mask hash = (key ^ (key >> 30)) * 0xbf58476d1ce4e5b9ULL;
  hash = (hash ^ (hash >> 27)) * 0x94d049bb133111ebULL;
  return hash ^ (hash >> 31);
#endif
}
template <unsigned Words> std::size_t hash_key(const WideMask<Words>& key) {
  std::size_t h = 0;
  for (auto x : key.w) h = hash_key(x ^ (h * 0x9e3779b97f4a7c15ULL));
  return h;
}

// Keep frontier states contiguous for iteration. The hash index stores entry
// indices plus one, so every clause mask is a valid key (including zero).
template <class Key = Mask> class StateTable {
  using Entry = std::pair<Key, Mask>;
  std::vector<Entry> entries;
  // Ganak caps a layer at a signed-int limit. Four slots per entry keep
  // collision chains short on correlated clause masks.
  std::vector<uint32_t> slots = std::vector<uint32_t>(16);

  std::size_t find_slot(const Key& key) const {
    std::size_t slot = hash_key(key) & (slots.size() - 1);
    while (slots[slot] && !(entries[slots[slot] - 1].first == key))
      slot = (slot + 1) & (slots.size() - 1);
    return slot;
  }

public:
  auto begin() const { return entries.begin(); }
  auto end() const { return entries.end(); }
  std::size_t size() const { return entries.size(); }
  void clear() {
    entries.clear();
    std::fill(slots.begin(), slots.end(), 0);
  }
  void add(const Key& key, Mask count) {
    std::size_t slot = find_slot(key);
    if (slots[slot]) {
      entries[slots[slot] - 1].second += count;
      return;
    }
    if (entries.size() >= slots.size() / 4) {
      slots.assign(slots.size() * 2, 0);
      for (std::size_t i = 0; i < entries.size(); ++i)
        slots[find_slot(entries[i].first)] = i + 1;
      slot = find_slot(key);
    }
    entries.emplace_back(key, count);
    slots[slot] = entries.size();
  }
};

// Bits in a DP key are open clauses whose processed literals are all false.
// With at most 62 variables, even an unconstrained count fits in uint64_t.
// CM holds one bit per clause.
template <class CM>
std::optional<Mask> sweep_impl(unsigned n, const std::vector<Clause>& clauses,
                               unsigned state_limit, uint64_t& peak, bool use_variable_bound) {
  std::array<CM, 62> pos{}, neg{}, others{};
  std::array<Mask, 62> neighbors{};
  for (unsigned j = 0; j < clauses.size(); ++j) {
    auto c = clauses[j];
    if (!c.vars()) return 0;
    const CM cb = clause_bit<CM>(j);
    for (Mask vs = c.vars(); vs; vs &= vs-1) {
      unsigned v = std::countr_zero(vs);
      if (use_variable_bound) neighbors[v] |= c.vars();
      if (c.pos & bit(v)) pos[v] |= cb;
      if (c.neg & bit(v)) neg[v] |= cb;
    }
  }
  CM seen{};
  Mask remaining = bit(n)-1;
  StateTable<CM> counts, next;
  counts.add(CM{}, 1);
  while (remaining) {
    // Prefix/suffix unions give the clauses remaining after each candidate.
    CM suffix{};
    for (unsigned v = n; v-- > 0;) {
      others[v] = suffix;
      if (remaining & bit(v)) suffix |= pos[v] | neg[v];
    }
    // A state is also determined by the assignments of processed variables
    // that share a clause with the suffix. A processed variable leaves this
    // boundary when its only remaining neighbor is the candidate.
    std::array<unsigned, 62> drops{};
    int boundary = 0;
    if (use_variable_bound) {
      for (Mask vs = (bit(n)-1) & ~remaining; vs; vs &= vs-1) {
        Mask future = neighbors[std::countr_zero(vs)] & remaining;
        if (future) ++boundary;
        if (std::has_single_bit(future)) ++drops[std::countr_zero(future)];
      }
    }
    CM prefix{};
    unsigned best = 0;
    std::tuple<int, int, int> best_score{std::numeric_limits<int>::max(), 0, 0};
    for (Mask vs = remaining; vs; vs &= vs-1) {
      unsigned v = std::countr_zero(vs);
      others[v] |= prefix;
      CM touched = pos[v] | neg[v];
      int width = popcount((seen | touched) & others[v]);
      int bound = width;
      if (use_variable_bound) {
        int next_boundary = boundary - int(drops[v]) + bool(neighbors[v] & (remaining ^ bit(v)));
        bound = std::min(width, next_boundary);
      }
      // Closing clauses can reject or merge states immediately. Give two
      // closures the same weight as reducing the frontier bound by one.
      int closed = popcount(touched & ~others[v]);
      auto score = std::tuple{2*bound - closed, width, -popcount(touched & seen)};
      if (score < best_score) { best_score = score; best = v; }
      prefix |= touched;
    }
    CM touched = pos[best] | neg[best];
    CM opened = touched & ~seen, closed = touched & ~others[best];
    next.clear();
    for (const auto& [owed, count] : counts) {
      for (const CM& satisfied : {pos[best], neg[best]}) {
        CM pending = (owed | opened) & ~satisfied;
        if (!(pending & closed)) next.add(pending, count);
      }
      if (next.size() > state_limit) return {};
    }
    peak = std::max<uint64_t>(peak, next.size());
    if (next.size() == 0) return 0;
    std::swap(counts, next);
    seen |= touched;
    remaining &= ~bit(best);
  }
  return counts.begin()->second; // All clauses have closed; only key zero survives.
}

constexpr unsigned max_sweep_clauses = 256;
constexpr unsigned wide_sweep_states = 128; // state cap with more than 64 clauses
inline std::optional<Mask> sweep(unsigned n, const std::vector<Clause>& clauses,
                                 unsigned state_limit, uint64_t& peak, bool use_variable_bound = false) {
  if (n > 62 || !state_limit) return {};
  const auto m = clauses.size();
  if (m <= 64) return sweep_impl<Mask>(n, clauses, state_limit, peak, use_variable_bound);
  if (m <= 128) return sweep_impl<WideMask<2>>(n, clauses, state_limit, peak, use_variable_bound);
  if (m <= max_sweep_clauses) return sweep_impl<WideMask<4>>(n, clauses, state_limit, peak, use_variable_bound);
  return {};
}
}
