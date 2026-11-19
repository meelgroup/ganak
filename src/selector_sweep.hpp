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

// A[S] counts admissible subcubes containing all points whose selectors are
// true in S. A block has at most 32 pairs and 20 selector rows: even after the
// zeta transform each entry is <= 3^32 + (2^20-1)*2^32, well within 64 bits.
inline std::vector<Word> block_table(unsigned width, const Block& block, Budget& budget) {
  const unsigned size = 1u << block.required.size();
  budget.spend(Word(size) * (1 + block.required.size()));
  Point full = width == 32 ? ~Point{0} : (Point{1} << width)-1;
  std::vector<Word> table(size);
  table[0] = 1;
  for (unsigned j = 0; j < width; ++j) table[0] *= 3;
  // Partition the inadmissible subcubes by the first forbidden point they contain.
  std::vector<Point> earlier;
  for (Point f : block.forbidden) {
    table[0] -= avoid(full, f, earlier, budget);
    earlier.push_back(f);
  }
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

class Encoding {
  unsigned vars;
  std::vector<Scope> positive, forbidden, targets;

public:
  explicit Encoding(unsigned n) : vars(n), targets(n+1) {}
  // 8 blocks of 32 pairs plus 20 selectors per block. Larger roots decline
  // before collecting any clauses. Every accepted clause has at most 32 literals.
  static bool fits(unsigned n, size_t clauses) { return n <= 672 && clauses <= 6800; }

  bool add(const std::vector<int>& clause) {
    if (clause.empty() || clause.size() > 32) return false;
    Scope pos, neg;
    for (int lit : clause) {
      if (!lit || lit < -int(vars) || lit > int(vars)) return false;
      (lit > 0 ? pos : neg).push_back(std::abs(lit));
    }
    if (neg.empty()) positive.push_back(std::move(pos));
    else if (pos.empty()) forbidden.push_back(std::move(neg));
    else if (pos.size() == 1 && neg.size() == 1) targets[neg[0]].push_back(pos[0]);
    else return false;
    return true;
  }

  std::optional<mpz_class> solve(RootStats* stats = nullptr, Word work_limit = 250000000) const {
    // Mixed binary clauses identify the selector variables. All other variables
    // must occur in one disjoint positive pair. Every selector must occur in
    // exactly one all-selector clause, which requires at least one true selector.
    std::vector<unsigned> pair(vars+1), side(vars+1), appearances(vars+1);
    std::vector<Scope> rows;
    unsigned pairs = 0;
    for (const auto& clause : positive) {
      unsigned selectors = 0;
      for (unsigned v : clause) selectors += !targets[v].empty();
      if (selectors == clause.size()) {
        rows.push_back(clause);
        for (unsigned v : clause) if (++appearances[v] > 1) return {};
      } else {
        if (selectors || clause.size() != 2) return {};
        unsigned a = clause[0], b = clause[1];
        if (a == b || pair[a] || pair[b]) return {};
        pair[a] = pair[b] = ++pairs;
        side[b] = 1;
      }
    }
    if (rows.empty() || rows.size() > 20) return {};
    const unsigned blocks = rows.front().size();
    if (blocks < 2 || blocks > 8) return {};
    for (unsigned v = 1; v <= vars; ++v)
      if (targets[v].empty() ? !pair[v] : appearances[v] != 1) return {};

    auto pattern = [&](const Scope& literals, Scope& scope, Point& value) {
      if (literals.size() < 2 || literals.size() > 32) return false;
      Scope ranked;
      for (unsigned v : literals) {
        if (!pair[v]) return false;
        ranked.push_back(2*pair[v] + side[v]);
      }
      std::sort(ranked.begin(), ranked.end());
      for (unsigned j = 0; j < ranked.size(); ++j) {
        unsigned p = ranked[j] / 2;
        if (j && p == scope.back()) return false;
        scope.push_back(p);
        if (ranked[j] % 2) value |= Point{1} << j;
      }
      return true;
    };
    std::map<Scope, unsigned> scopes;
    std::vector<bool> used(pairs+1);
    std::vector<unsigned> block_of(vars+1);
    std::vector<Point> point_of(vars+1);
    std::vector<Block> factors(blocks);
    for (unsigned v = 1; v <= vars; ++v) {
      if (targets[v].empty()) continue;
      Scope scope;
      if (!pattern(targets[v], scope, point_of[v])) return {};
      auto [it, fresh] = scopes.try_emplace(scope, scopes.size());
      if (fresh) {
        if (scopes.size() > blocks) return {};
        for (unsigned p : scope) {
          if (used[p]) return {};
          used[p] = true;
        }
      }
      block_of[v] = it->second;
    }
    if (scopes.size() != blocks || std::count(used.begin(), used.end(), true) != pairs) return {};
    const unsigned width = scopes.begin()->first.size();
    for (const auto& entry : scopes) if (entry.first.size() != width) return {};
    for (const auto& row : rows) {
      if (row.size() != blocks) return {};
      unsigned seen = 0;
      for (unsigned v : row) {
        unsigned b = block_of[v];
        if (seen & (1u << b)) return {};
        seen |= 1u << b;
        factors[b].required.push_back(point_of[v]);
      }
    }
    for (const auto& clause : forbidden) {
      Scope scope;
      Point value = 0;
      if (!pattern(clause, scope, value)) return {};
      auto it = scopes.find(scope);
      if (it == scopes.end()) return {};
      factors[it->second].forbidden.push_back(value);
    }
    for (auto& f : factors) {
      // Pair orientation and order are arbitrary. Align columns using selector
      // rows before comparing whole blocks and sharing their transformed tables.
      const Point anchor = f.required.front();
      std::vector<std::pair<Word, unsigned>> columns;
      for (unsigned j = 0; j < width; ++j) {
        Word signature = 0;
        for (unsigned r = 0; r < rows.size(); ++r)
          signature |= Word(((f.required[r] ^ anchor) >> j) & 1) << r;
        columns.emplace_back(signature, j);
      }
      std::sort(columns.begin(), columns.end());
      auto align = [&](Point p) {
        Point result = 0;
        for (unsigned j = 0; j < width; ++j)
          result |= (((p ^ anchor) >> columns[j].second) & 1) << j;
        return result;
      };
      for (auto& p : f.required) p = align(p);
      for (auto& p : f.forbidden) p = align(p);
      auto& fs = f.forbidden;
      std::sort(fs.begin(), fs.end());
      fs.erase(std::unique(fs.begin(), fs.end()), fs.end());
      if (fs.size() > 128) return {};
    }

    Budget budget{work_limit};
    try {
      std::map<Block, unsigned> ids;
      std::vector<std::vector<Word>> tables;
      std::vector<unsigned> layout;
      for (const auto& f : factors) {
        auto [it, fresh] = ids.try_emplace(f, tables.size());
        if (fresh) tables.push_back(block_table(width, f, budget));
        layout.push_back(it->second);
      }
      // OR-convolution at the full selector mask by subset inversion. Selectors
      // are optional witnesses, so all nonempty choices in each row are counted.
      budget.spend(tables.front().size() * layout.size());
      mpz_class answer = 0, term, value;
      for (unsigned s = 0; s < tables.front().size(); ++s) {
        term = 1;
        for (unsigned b : layout) {
          const Word entry = tables[b][s];
          mpz_import(value.get_mpz_t(), 1, 1, sizeof(entry), 0, 0, &entry);
          term *= value;
        }
        if ((rows.size()-std::popcount(s)) % 2) answer -= term;
        else answer += term;
      }
      if (stats) *stats = {blocks, width, unsigned(rows.size()), work_limit-budget.left};
      return answer;
    } catch (const LimitReached&) {
      return {};
    }
  }
};
}
