#include "src/frontier_tables.hpp"
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

using namespace GanakInt::frontier;

static void check(std::optional<Mask> result, Mask expected, const char* method,
                  unsigned trial, bool required = true) {
  if ((!result && required) || (result && *result != expected))
    throw std::runtime_error(std::string(method) + " mismatch on trial " + std::to_string(trial));
}

static Mask brute(unsigned n, const std::vector<Clause>& clauses) {
  Mask total = 0;
  for (Mask assignment = 0; assignment < bit(n); ++assignment) {
    bool valid = true;
    for (auto c : clauses) {
      if (!(assignment & c.pos) && !(~assignment & c.neg)) { valid = false; break; }
    }
    total += valid;
  }
  return total;
}

int main() {
  std::mt19937_64 rng(20261118);
  for (unsigned trial = 0; trial < 1000; ++trial) {
    unsigned n = 1 + rng()%12, m = rng()%129;
    std::vector<Clause> clauses(m);
    for (auto& c : clauses) {
      unsigned arity = 1 + rng()%6;
      for (unsigned j = 0; j < arity; ++j)
        (rng()%2 ? c.pos : c.neg) |= bit(rng()%n);
    }
    Mask expected = brute(n, clauses);
    check(word_count(n, clauses, 8192), expected, "word", trial);
    check(bucket_count(n, clauses, 8192), expected, "bucket", trial);
    check(bucket_count(n, clauses, 2), expected, "bounded bucket", trial, false);
    if (m <= 64) {
      for (bool variable_bound : {false, true}) {
        uint64_t peak = 0;
        check(sweep(n, clauses, 8192, peak, nullptr, variable_bound), expected, "sparse", trial);
        check(sweep(n, clauses, 2, peak, nullptr, variable_bound), expected, "bounded sparse", trial, false);
      }
    }
  }

  // Nonconsecutive original variable numbers exercise message projection
  // when the elimination order is different from variable-number order.
  std::vector<Clause> path;
  for (unsigned i = 1; i < 62; ++i)
    path.push_back({bit((17*i)%62) | bit((17*(i-1))%62), 0});
  Mask prev = 1, curr = 2;
  for (unsigned i = 2; i <= 62; ++i) { Mask next = prev+curr; prev=curr; curr=next; }
  check(bucket_count(62, path, 8192), curr, "path", 0);
  check(bucket_count(62, {}, 8192), bit(62), "free", 0);
  check(bucket_count(62, {Clause{}}, 8192), 0, "empty clause", 0);
  check(word_count(16, {{bit(15), 0}, {0, bit(0)}}, 8192), bit(14), "word blocks", 0);
  if (bucket_count(30, {{bit(30)-1, 0}}, 8192))
    throw std::runtime_error("wide clique exceeded the table limit");
  std::cout << "Frontier count checks passed\n";
}
