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
#include <vector>

namespace GanakInt::frontier {
using Mask = uint64_t;
struct Clause {
  Mask pos = 0, neg = 0;
  Mask vars() const { return pos | neg; }
};
inline Mask bit(unsigned i) { return Mask{1} << i; }

// Keep frontier states contiguous for iteration. The hash index stores entry
// indices plus one, so every 64-bit clause mask is a valid key (including zero).
class StateTable {
  using Entry = std::pair<Mask, Mask>;
  std::vector<Entry> entries;
  // Ganak caps a layer at a signed-int limit. Four slots per entry keep
  // collision chains short on correlated clause masks.
  std::vector<uint32_t> slots = std::vector<uint32_t>(16);

  std::size_t find_slot(Mask key) const {
    // Clause masks often differ only in their high bits.
#if defined(__SSE4_2__) && defined(__x86_64__)
    std::size_t slot = _mm_crc32_u64(0, key) & (slots.size() - 1);
#else
    Mask hash = (key ^ (key >> 30)) * 0xbf58476d1ce4e5b9ULL;
    hash = (hash ^ (hash >> 27)) * 0x94d049bb133111ebULL;
    std::size_t slot = (hash ^ (hash >> 31)) & (slots.size() - 1);
#endif
    while (slots[slot] && entries[slots[slot] - 1].first != key)
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
  void add(Mask key, Mask count) {
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
inline std::optional<Mask> sweep(unsigned n, const std::vector<Clause>& clauses,
                                 unsigned state_limit, uint64_t& peak, uint64_t* visited = nullptr, bool use_variable_bound = false) {
  if (clauses.size() > 64 || n > 62 || !state_limit) return {};
  std::array<Mask, 62> pos{}, neg{}, others, neighbors{};
  for (unsigned j = 0; j < clauses.size(); ++j) {
    auto c = clauses[j];
    if (!c.vars()) return 0;
    for (Mask vs = c.vars(); vs; vs &= vs-1) {
      unsigned v = std::countr_zero(vs);
      if (use_variable_bound) neighbors[v] |= c.vars();
      if (c.pos & bit(v)) pos[v] |= bit(j);
      if (c.neg & bit(v)) neg[v] |= bit(j);
    }
  }
  Mask seen = 0, remaining = bit(n)-1;
  StateTable counts, next;
  counts.add(0, 1);
  while (remaining) {
    // Prefix/suffix unions give the clauses remaining after each candidate.
    Mask suffix = 0;
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
    Mask prefix = 0;
    unsigned best = 0;
    std::tuple<int, int, int> best_score{65, 65, 0};
    for (Mask vs = remaining; vs; vs &= vs-1) {
      unsigned v = std::countr_zero(vs);
      others[v] |= prefix;
      Mask touched = pos[v] | neg[v];
      int width = std::popcount((seen | touched) & others[v]);
      int bound = width;
      if (use_variable_bound) {
        int next_boundary = boundary - int(drops[v]) + bool(neighbors[v] & (remaining ^ bit(v)));
        bound = std::min(width, next_boundary);
      }
      // Closing clauses can reject or merge assignments on this very layer.
      auto score = std::tuple{bound, width, -std::popcount(touched & ~others[v])};
      if (score < best_score) { best_score = score; best = v; }
      prefix |= touched;
    }
    Mask touched = pos[best] | neg[best];
    Mask opened = touched & ~seen, closed = touched & ~others[best];
    next.clear();
    if (visited) *visited += counts.size();
    for (auto [owed, count] : counts) {
      for (Mask satisfied : {pos[best], neg[best]}) {
        Mask pending = (owed | opened) & ~satisfied;
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

// Find one- or two-variable separators. Require at least a quarter of the
// remaining variables outside the largest piece to justify overriding Ganak.
inline int separator(unsigned n, const std::vector<Clause>& clauses) {
  std::vector<Mask> adj(n);
  for (auto c : clauses)
    for (Mask vs = c.vars(); vs; vs &= vs-1) adj[std::countr_zero(vs)] |= c.vars();
  std::vector<unsigned> ranked(n);
  for (unsigned i=0; i<n; ++i) ranked[i]=i;
  std::stable_sort(ranked.begin(), ranked.end(), [&](unsigned a, unsigned b) {
    return std::popcount(adj[a]) > std::popcount(adj[b]);
  });
  auto largest = [&](Mask removed) {
    Mask left = (bit(n)-1) & ~removed;
    int worst = 0;
    while (left) {
      Mask queue = bit(std::countr_zero(left));
      int size = 0;
      while (queue) {
        left &= ~queue;
        size += std::popcount(queue);
        Mask neighbors = 0;
        for (; queue; queue &= queue-1) neighbors |= adj[std::countr_zero(queue)];
        queue = neighbors & left;
      }
      worst = std::max(worst, size);
    }
    return worst;
  };
  int best = -1, bound = int(3*(n-1)/4)+1;
  for (auto v : ranked) {
    int width = largest(bit(v));
    if (width < bound) { bound=width; best=int(v); }
  }
  if (best >= 0) return best;
  bound = int(3*(n-2)/4)+1;
  const unsigned limit = std::min(n, 10u);
  for (unsigned i=0; i<limit; ++i) for (unsigned j=i+1; j<limit; ++j) {
    int width = largest(bit(ranked[i]) | bit(ranked[j]));
    if (width < bound) { bound=width; best=int(ranked[i]); }
  }
  return best;
}
}
