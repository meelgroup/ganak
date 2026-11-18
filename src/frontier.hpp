// Exact bounded component algorithms. Local variable i is bit i in a clause.
#pragma once
#include <algorithm>
#include <bit>
#include <cstdint>
#include <limits>
#include <optional>
#include <unordered_map>
#include <vector>

namespace GanakInt::frontier {
using Mask = uint64_t;
struct Clause {
  Mask pos = 0, neg = 0;
  Mask vars() const { return pos | neg; }
};
inline Mask bit(unsigned i) { return Mask{1} << i; }

// Bits in a DP key are open clauses whose processed literals are all false.
// With at most 62 variables, even an unconstrained count fits in uint64_t.
inline std::optional<Mask> sweep(unsigned n, const std::vector<Clause>& clauses,
                                 unsigned state_limit, uint64_t& peak, uint64_t* visited = nullptr) {
  if (clauses.size() > 64 || n > 62 || !state_limit) return {};
  std::vector<Mask> pos(n), neg(n), others(n);
  for (unsigned j = 0; j < clauses.size(); ++j) {
    auto c = clauses[j];
    if (!c.vars()) return 0;
    for (Mask vs = c.vars(); vs; vs &= vs-1) {
      unsigned v = std::countr_zero(vs);
      if (c.pos & bit(v)) pos[v] |= bit(j);
      if (c.neg & bit(v)) neg[v] |= bit(j);
    }
  }
  Mask seen = 0, remaining = bit(n)-1;
  std::unordered_map<Mask, Mask> counts{{0,1}}, next;
  while (remaining) {
    // Prefix/suffix unions give the clauses remaining after each candidate.
    Mask suffix = 0;
    for (unsigned v = n; v-- > 0;) {
      others[v] = suffix;
      if (remaining & bit(v)) suffix |= pos[v] | neg[v];
    }
    Mask prefix = 0;
    unsigned best = 0;
    int best_width = std::numeric_limits<int>::max();
    for (Mask vs = remaining; vs; vs &= vs-1) {
      unsigned v = std::countr_zero(vs);
      others[v] |= prefix;
      Mask touched = pos[v] | neg[v];
      int width = std::popcount((seen | touched) & others[v]);
      if (width < best_width) { best_width = width; best = v; }
      prefix |= touched;
    }
    Mask touched = pos[best] | neg[best];
    Mask opened = touched & ~seen, closed = touched & ~others[best];
    next.clear();
    if (visited) *visited += counts.size();
    for (auto [owed, count] : counts) {
      for (Mask satisfied : {pos[best], neg[best]}) {
        Mask pending = (owed | opened) & ~satisfied;
        if (!(pending & closed)) next[pending] += count;
      }
      if (next.size() > state_limit) return {};
    }
    peak = std::max<uint64_t>(peak, next.size());
    if (next.empty()) return 0;
    counts.swap(next);
    seen |= touched;
    remaining &= ~bit(best);
  }
  return counts.at(0);
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
