// Exact root counting for disjoint blocks of positive variable pairs linked by
// selector clauses. Each pair permits the first value, the second, or both.
// A block assignment is therefore a subcube of binary points.
#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <compare>
#include <functional>
#include <cstdlib>
#include <map>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>
#include <gmpxx.h>

namespace GanakInt::selector {
using Point = uint32_t;
using Word = uint64_t;
using Scope = std::vector<unsigned>;

struct Block {
  std::vector<Point> forbidden, required;
  auto operator<=>(const Block&) const = default;
};

struct RootStats {
  unsigned blocks = 0, pairs_per_block = 0, rows = 0;
  Word work = 0;
  bool witness_search = false;
};

struct LimitReached {};
struct Budget {
  Word left;
  void spend(Word work) {
    if (work > left) throw LimitReached{};
    left -= work;
  }
};

// At most 128 forbidden points contribute to an intersection.
struct Rows {
  Word low = 0, high = 0;
  bool operator==(const Rows&) const = default;
  void add(unsigned j) {
    if (j < 64) low |= Word{1} << j;
    else high |= Word{1} << (j-64);
  }
  Rows operator&(Rows other) const { return {low & other.low, high & other.high}; }
};
struct RowsHash {
  size_t operator()(Rows r) const {
    Word h = r.low * 0x9e3779b97f4a7c15ULL + std::rotl(r.high, 37) * 0xbf58476d1ce4e5b9ULL;
    return h ^ (h >> 29);
  }
};
using Histogram = std::unordered_map<Rows, Word, RowsHash>;

// For each choice of fixed coordinates on this half, track the forbidden points
// that still agree with those fixed values. Group choices with the same set.
inline Histogram histogram(Rows all, const std::vector<Rows>& columns,
                           unsigned begin, unsigned end, Budget& budget) {
  Histogram counts{{all, 1}}, next;
  for (unsigned j = begin; j < end; ++j) {
    budget.spend(2 * counts.size());
    next = counts;
    for (auto [rows, amount] : counts) next[rows & columns[j]] += amount;
    counts.swap(next);
  }
  return counts;
}

// Inclusion-exclusion on at most seven forbidden subcubes. A choice of fixed
// coordinates C is forbidden by point f exactly when C is a subset of agree[f].
inline Word avoid_few(Point available, const Point* agree, unsigned n) {
  std::array<Point, 128> intersection;
  intersection[0] = available;
  int64_t total = int64_t{1} << std::popcount(available);
  for (unsigned s = 1; s < (1u << n); ++s) {
    intersection[s] = intersection[s & (s-1)] & agree[std::countr_zero(s)];
    int64_t term = int64_t{1} << std::popcount(intersection[s]);
    total += std::popcount(s) % 2 ? -term : term;
  }
  return total;
}

// Positive points have already fixed which coordinates may be constrained and
// their common bit values. Count ways of choosing constraints that exclude all
// forbidden points in this block.
inline Word avoid(Point available, Point values, const std::vector<Point>& forbidden,
                  Budget& budget) {
  budget.spend(forbidden.size() + 1);
  std::vector<Point> agree;
  for (Point f : forbidden) {
    Point mask = available & ~(f ^ values);
    if (mask == available) return 0;
    agree.push_back(mask);
  }
  std::sort(agree.begin(), agree.end(), std::greater<>{});
  agree.erase(std::unique(agree.begin(), agree.end()), agree.end());
  budget.spend(agree.size()*agree.size());
  std::vector<Point> maximal;
  for (Point mask : agree) {
    if (std::none_of(maximal.begin(), maximal.end(), [mask](Point other) {
          return (mask & other) == mask;
        })) maximal.push_back(mask);
  }
  if (maximal.size() <= 7) {
    budget.spend(Word{1} << maximal.size());
    return avoid_few(available, maximal.data(), maximal.size());
  }

  Rows all;
  for (unsigned j = 0; j < maximal.size(); ++j) all.add(j);
  std::vector<Rows> columns;
  budget.spend(std::popcount(available) * maximal.size());
  const unsigned half = std::popcount(available) / 2;
  Point right_mask = available;
  for (Point bits = available; bits; bits &= bits-1) {
    Point coordinate = bits & -bits;
    if (columns.size() < half) right_mask ^= coordinate;
    Rows matching;
    for (unsigned j = 0; j < maximal.size(); ++j)
      if (maximal[j] & coordinate) matching.add(j);
    columns.push_back(matching);
  }
  auto left = histogram(all, columns, 0, half, budget);
  auto right = histogram(all, columns, half, columns.size(), budget);
  Word total = 0;
  for (auto [x, multiplicity] : left) {
    unsigned n = std::popcount(x.low) + std::popcount(x.high);
    Word valid = 0;
    if (n <= 7) {
      budget.spend(Word{1} << n);
      std::array<Point, 7> masks;
      unsigned k = 0;
      for (Word v = x.low; v; v &= v-1) masks[k++] = maximal[std::countr_zero(v)];
      for (Word v = x.high; v; v &= v-1) masks[k++] = maximal[64+std::countr_zero(v)];
      valid = avoid_few(right_mask, masks.data(), n);
    } else {
      budget.spend(right.size());
      for (auto [y, amount] : right)
        if (!(x.low & y.low) && !(x.high & y.high)) valid += amount;
    }
    total += multiplicity * valid;
  }
  return total;
}

inline Word unconstrained_count(unsigned width, Point full, const std::vector<Point>& forbidden,
                                Budget& budget) {
  Word count = 1;
  for (unsigned j = 0; j < width; ++j) count *= 3;
  // Partition the inadmissible subcubes by the first forbidden point they contain.
  std::vector<Point> earlier;
  for (Point f : forbidden) {
    count -= avoid(full, f, earlier, budget);
    earlier.push_back(f);
  }
  return count;
}

// A[S] counts admissible subcubes containing all points whose selectors are
// true in S. A block has at most 32 pairs and 20 selector rows: even after the
// zeta transform each entry is <= 3^32 + (2^20-1)*2^32, well within 64 bits.
inline std::vector<Word> block_table(unsigned width, const Block& block, Budget& budget) {
  const unsigned size = 1u << block.required.size();
  budget.spend(Word(size) * (1 + block.required.size()));
  Point full = width == 32 ? ~Point{0} : (Point{1} << width)-1;
  std::vector<Word> table(size);
  table[0] = unconstrained_count(width, full, block.forbidden, budget);
  std::vector<Point> ones(size, full), zeros(size, full);
  std::unordered_map<Word, Word> cache;
  for (unsigned s = 1; s < size; ++s) {
    Point p = block.required[std::countr_zero(s)];
    unsigned previous = s & (s-1);
    ones[s] = ones[previous] & p;
    zeros[s] = zeros[previous] & ~p;
    auto [it, fresh] = cache.try_emplace((Word(ones[s]) << 32) | zeros[s], 0);
    if (fresh) it->second = avoid(ones[s] | zeros[s], ones[s], block.forbidden, budget);
    table[s] = it->second;
  }
  for (unsigned step = 1; step < size; step *= 2)
    for (unsigned base = 0; base < size; base += 2*step)
      for (unsigned j = 0; j < step; ++j) table[base+step+j] += table[base+j];
  return table;
}

}
