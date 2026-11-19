// Exact counting when a dense table over selector rows would exceed the work
// budget. Each nonempty witness choice imposes a hull on the selected blocks.
#pragma once
#include "selector_tables.hpp"

namespace GanakInt::selector {

class HullTable {
  struct Hull {
    Point ones, zeros;
    std::array<int, 32> next;
    std::optional<Word> ways;
  };
  static constexpr int unknown = -2, impossible = -1;
  unsigned width;
  Point full;
  const Block block;
  Budget& budget;
  std::vector<Hull> hulls;
  std::unordered_map<Word, int> ids;

  int intern(Point ones, Point zeros) {
    Word key = (Word(ones) << 32) | zeros;
    if (auto it = ids.find(key); it != ids.end()) return it->second;
    // Bound the memory for rejected and feasible hulls separately.
    if (ids.size() >= 131072) throw LimitReached{};
    budget.spend(1 + block.forbidden.size());
    // These are the coordinates shared by every required point. Any forbidden
    // point agreeing on all of them would lie even in the smallest possible cube.
    for (Point f : block.forbidden)
      if (!((~f & ones) | (f & zeros))) {
        ids.emplace(key, impossible);
        return impossible;
      }
    if (hulls.size() >= 32768) throw LimitReached{};
    int id = int(hulls.size());
    Hull h{ones, zeros, {}, {}};
    h.next.fill(unknown);
    hulls.push_back(h);
    ids.emplace(key, id);
    return id;
  }
public:
  HullTable(unsigned w, const Block& f, Budget& work)
      : width(w), full(w == 32 ? ~Point{0} : (Point{1} << w)-1), block(f), budget(work) {
    intern(full, full); // empty selection, no required point yet
  }

  int add_row(int id, unsigned row) {
    int result = hulls[id].next[row];
    if (result == unknown) {
      Point p = block.required[row];
      result = intern(hulls[id].ones & p, hulls[id].zeros & ~p);
      hulls[id].next[row] = result;
    }
    return result;
  }

  Word count(int id) {
    auto& h = hulls[id];
    if (!h.ways)
      h.ways = id ? avoid(h.ones | h.zeros, h.ones, block.forbidden, budget)
                  : unconstrained_count(width, full, block.forbidden, budget);
    return *h.ways;
  }
};

class WitnessSearch {
  using Ids = std::array<int, 3>;
  struct Key {
    Point remaining;
    Ids ids;
    bool operator==(const Key&) const = default;
  };
  struct Hash {
    size_t operator()(const Key& key) const {
      Word hash = key.remaining;
      for (int id : key.ids) hash = (hash ^ unsigned(id)) * 0x9e3779b97f4a7c15ULL;
      return hash ^ (hash >> 31);
    }
  };
  Budget& budget;
  unsigned blocks;
  unsigned calls = 0;
  std::map<Block, HullTable> distinct;
  std::array<HullTable*, 3> tables{};
  std::unordered_map<Key, mpz_class, Hash> memo;

  mpz_class visit(Point remaining, Ids ids) {
    if (++calls > 650000) throw LimitReached{};
    budget.spend(1);
    for (unsigned i = 0; i < blocks; ++i)
      for (unsigned j = i+1; j < blocks; ++j)
        if (tables[i] == tables[j] && ids[i] > ids[j]) std::swap(ids[i], ids[j]);
    Key key{remaining, ids};
    if (auto it = memo.find(key); it != memo.end()) return it->second;

    std::pair<int, int> best{4, 4};
    unsigned selected = 0, legal = 0, same = 0;
    Ids next{};
    mpz_class factor = 1;
    budget.spend(std::popcount(remaining) * blocks);
    for (Point bits = remaining; bits; bits &= bits-1) {
      unsigned r = std::countr_zero(bits), choices = 0, existing = 0;
      Ids steps{};
      for (unsigned b = 0; b < blocks; ++b) {
        steps[b] = tables[b]->add_row(ids[b], r);
        if (steps[b] >= 0) choices |= 1u << b;
        if (steps[b] == ids[b]) existing |= 1u << b;
      }
      if (!choices) return 0;
      if (choices == existing) {
        // Future hulls can only grow. The same witnesses stay available, and
        // a block already excluded by a forbidden point cannot become feasible.
        remaining ^= Point{1} << r;
        factor *= (1u << std::popcount(existing))-1;
        continue;
      }
      // Minimize the feasible blocks first, then the blocks whose hull changes.
      auto score = std::pair{std::popcount(choices), std::popcount(choices ^ existing)};
      if (score < best) {
        best = score; selected = r; legal = choices; same = existing; next = steps;
      }
    }

    mpz_class result = 0;
    if (!remaining) {
      result = 1;
      mpz_class value;
      for (unsigned b = 0; b < blocks; ++b) {
        Word ways = tables[b]->count(ids[b]);
        mpz_import(value.get_mpz_t(), 1, 1, sizeof(ways), 0, 0, &ways);
        result *= value;
      }
    } else {
      remaining ^= Point{1} << selected;
      unsigned multiplicity = 1u << std::popcount(same);
      if (same) result += (multiplicity-1) * visit(remaining, ids);
      unsigned changes = legal ^ same;
      for (unsigned subset = changes; subset; subset = (subset-1) & changes) {
        Ids child = ids;
        for (unsigned b = 0; b < blocks; ++b)
          if (subset & (1u << b)) child[b] = next[b];
        result += multiplicity * visit(remaining, child);
      }
    }
    result *= factor;
    if (memo.size() < 120000) memo.emplace(key, result);
    return result;
  }

public:
  WitnessSearch(unsigned width, const std::vector<Block>& factors, Budget& work)
      : budget(work), blocks(factors.size()) {
    if (blocks < 2 || blocks > tables.size()) throw LimitReached{};
    for (unsigned b = 0; b < blocks; ++b) {
      const auto& f = factors[b];
      tables[b] = &distinct.try_emplace(f, width, f, budget).first->second;
    }
  }

  mpz_class count(unsigned rows) {
    if (!rows || rows > 32) throw LimitReached{};
    return visit(rows == 32 ? ~Point{0} : (Point{1} << rows)-1, {});
  }
};
}
