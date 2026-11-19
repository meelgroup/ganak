#include "src/selector_sweep.hpp"
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>

using GanakInt::selector::Budget;
using GanakInt::selector::Encoding;
using GanakInt::selector::Point;
using GanakInt::selector::Word;

static void check_block(std::mt19937& rng) {
  // Many distinct equal-sized agreement masks exercise both words of the
  // forbidden-point set and the meet-in-the-middle histogram.
  constexpr unsigned width = 12;
  constexpr Point full = (1u << width)-1;
  for (unsigned trial = 0; trial < 80; ++trial) {
    std::vector<Point> forbidden;
    Point base = rng() & full;
    for (Point d = 1; d <= full; ++d)
      if (std::popcount(d) == 6 && rng()%6 == 0 && forbidden.size() < 120)
        forbidden.push_back(base ^ d);
    Point available = trial%2 ? full : rng() & full;
    Word expected = 0;
    for (Point fixed = available;; fixed = (fixed-1) & available) {
      if (std::none_of(forbidden.begin(), forbidden.end(), [=](Point f) {
            return ((base ^ f) & fixed) == 0;
          })) ++expected;
      if (!fixed) break;
    }
    Budget budget{100000000};
    if (GanakInt::selector::avoid(available, base, forbidden, budget) != expected)
      throw std::runtime_error("incorrect forbidden-point count");
  }
}

static void check_encoding(std::mt19937& rng) {
  for (unsigned trial = 0; trial < 100; ++trial) {
    constexpr unsigned width = 2, blocks = 2, features = 2*width*blocks;
    unsigned rows = 1 + rng()%3, vars = features + rows*blocks;
    std::vector<std::vector<int>> cnf;
    for (unsigned v = 1; v <= features; v += 2)
      cnf.push_back({int(v), int(v+1)});
    for (unsigned b = 0; b < blocks; ++b) {
      for (unsigned p = 0; p < (1u << width); ++p) {
        if (rng()%3) continue;
        std::vector<int> forbidden;
        for (unsigned j = 0; j < width; ++j)
          forbidden.push_back(-int(1+2*(b*width+j)+((p >> j) & 1)));
        cnf.push_back(std::move(forbidden));
      }
    }
    for (unsigned row = 0; row < rows; ++row) {
      std::vector<int> required;
      for (unsigned b = 0; b < blocks; ++b) {
        int z = features + row*blocks + b + 1;
        required.push_back(z);
        for (unsigned j = 0; j < width; ++j)
          cnf.push_back({-z, int(1+2*(b*width+j)+rng()%2)});
      }
      cnf.push_back(std::move(required));
    }
    std::vector<int> perm(vars);
    std::iota(perm.begin(), perm.end(), 1);
    std::shuffle(perm.begin(), perm.end(), rng);
    std::shuffle(cnf.begin(), cnf.end(), rng);
    Encoding encoding(vars);
    for (auto& clause : cnf) {
      for (int& lit : clause) lit = (lit < 0 ? -1 : 1) * perm[std::abs(lit)-1];
      std::shuffle(clause.begin(), clause.end(), rng);
      if (!encoding.add(clause)) throw std::runtime_error("rejected valid clause");
    }
    unsigned expected = 0;
    for (unsigned x = 0; x < (1u << vars); ++x) {
      bool good = true;
      for (const auto& clause : cnf) {
        bool satisfied = false;
        for (int lit : clause)
          satisfied |= bool(x & (1u << (std::abs(lit)-1))) == (lit > 0);
        if (!satisfied) { good = false; break; }
      }
      expected += good;
    }
    auto actual = encoding.solve();
    if (!actual || *actual != expected || encoding.solve(nullptr, 0))
      throw std::runtime_error("incorrect selector count or budget fallback");
    // An additional binary implication between features invalidates the encoding.
    if (encoding.add({perm[0], -perm[2]}) && encoding.solve())
      throw std::runtime_error("accepted an unsupported cross-block clause");
  }
}

// Count by enumerating the 3^width actual subcubes in each block, then
// counting nonempty choices of witnesses for each row independently.
static mpz_class enumerate_witnesses(unsigned width,
                                    const std::vector<GanakInt::selector::Block>& blocks) {
  std::vector<std::vector<Point>> masks(blocks.size());
  const Point full = (1u << width)-1;
  for (unsigned b = 0; b < blocks.size(); ++b) {
    const auto& block = blocks[b];
    for (Point fixed = 0; fixed <= full; ++fixed) {
      for (Point value = fixed;; value = (value-1) & fixed) {
        bool legal = std::none_of(block.forbidden.begin(), block.forbidden.end(),
            [&](Point p) { return (p & fixed) == value; });
        if (legal) {
          Point mask = 0;
          for (unsigned r = 0; r < block.required.size(); ++r)
            if ((block.required[r] & fixed) == value) mask |= Point{1} << r;
          masks[b].push_back(mask);
        }
        if (!value) break;
      }
    }
  }
  if (blocks.size() == 2) masks.push_back({0});
  mpz_class sum = 0;
  for (Point a : masks[0]) for (Point b : masks[1]) for (Point c : masks[2]) {
    mpz_class ways = 1;
    for (unsigned r = 0; r < blocks[0].required.size() && ways != 0; ++r) {
      unsigned count = ((a >> r) & 1) + ((b >> r) & 1) + ((c >> r) & 1);
      ways *= (1u << count)-1;
    }
    sum += ways;
  }
  return sum;
}

static void check_search(std::mt19937& rng) {
  using GanakInt::selector::Block;
  using GanakInt::selector::WitnessSearch;
  for (unsigned trial = 0; trial < 120; ++trial) {
    const unsigned width = trial%10 == 0 ? 3 : 2;
    const unsigned rows = 1 + trial%32;
    std::vector<Block> blocks(2 + trial%2);
    for (unsigned b = 0; b < blocks.size(); ++b) {
      if (b && trial%3 == 0) { blocks[b] = blocks.front(); continue; }
      for (unsigned p = 0; p < (1u << width); ++p)
        if (rng()%4 == 0) blocks[b].forbidden.push_back(p);
      for (unsigned r = 0; r < rows; ++r)
        blocks[b].required.push_back(rng() % (1u << width));
    }
    Budget budget{250000000};
    WitnessSearch search(width, blocks, budget);
    if (search.count(rows) != enumerate_witnesses(width, blocks))
      throw std::runtime_error("incorrect hull search count");
  }

  // The normal recognizer must choose the search for a 32-row CNF.
  constexpr unsigned rows = 32;
  Encoding encoding(8 + 2*rows);
  for (int v = 1; v <= 8; v += 2) encoding.add({v, v+1});
  std::vector<Block> blocks(2);
  for (unsigned r = 0; r < rows; ++r) {
    int z = int(9 + 2*r);
    encoding.add({z, z+1});
    for (unsigned b = 0; b < 2; ++b) {
      Point p = (r+b) % 4;
      blocks[b].required.push_back(p);
      for (unsigned j = 0; j < 2; ++j)
        encoding.add({-(z+int(b)), int(1 + 4*b + 2*j + ((p >> j) & 1))});
    }
  }
  auto actual = encoding.solve();
  if (!actual || *actual != enumerate_witnesses(2, blocks) || encoding.solve(nullptr, 1))
    throw std::runtime_error("incorrect 32-row dispatch or fallback");
}

int main() {
  std::mt19937 rng(20261119);
  check_block(rng);
  check_encoding(rng);
  check_search(rng);
  std::cout << "Selector count checks passed\n";
}
