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

int main() {
  std::mt19937 rng(20261119);
  check_block(rng);
  check_encoding(rng);
  std::cout << "Selector count checks passed\n";
}
