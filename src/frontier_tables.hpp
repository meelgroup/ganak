// Exact bounded tables for unweighted residuals with local variable masks.
#pragma once
#include "frontier.hpp"
#if defined(__BMI2__) && defined(__x86_64__)
#include <immintrin.h>
#endif

namespace GanakInt::frontier {

// One word represents the low six variables, one bit per assignment. Only
// blocks in which all the higher literals are false need a clause update.
inline std::optional<Mask> word_count(unsigned n, const std::vector<Clause>& clauses,
                                      unsigned state_limit) {
  if (n > 18) return {};
  const unsigned blocks = 1u << (n > 6 ? n-6 : 0);
  if (blocks > state_limit) return {};
  if (n > 12) {
    unsigned work = blocks;
    for (auto c : clauses) {
      work += blocks >> std::popcount(c.vars() >> 6);
      if (work > 4096) return {};
    }
  }
  constexpr Mask truth[] = {
    0xaaaaaaaaaaaaaaaaULL, 0xccccccccccccccccULL, 0xf0f0f0f0f0f0f0f0ULL,
    0xff00ff00ff00ff00ULL, 0xffff0000ffff0000ULL, 0xffffffff00000000ULL
  };
  std::array<Mask, 4096> allowed;
  std::fill_n(allowed.begin(), blocks, ~Mask{0});
  for (auto c : clauses) {
    if (c.pos & c.neg) continue;
    Mask low = 0;
    for (Mask vs = c.pos & 63; vs; vs &= vs-1) low |= truth[std::countr_zero(vs)];
    for (Mask vs = c.neg & 63; vs; vs &= vs-1) low |= ~truth[std::countr_zero(vs)];
    const unsigned fixed = c.neg >> 6;
    const unsigned free = (blocks-1) & ~(c.vars() >> 6);
    for (unsigned subset = free;; subset = (subset-1) & free) {
      allowed[fixed | subset] &= low;
      if (!subset) break;
    }
  }
  Mask total = 0;
  for (unsigned i = 0; i < blocks; ++i) total += std::popcount(allowed[i]);
  // With fewer than six variables each assignment appears repeatedly.
  return n < 6 ? total >> (6-n) : total;
}

// Pack selected bits in order, with the first selected bit in position zero.
inline Mask extract_bits(Mask value, Mask selected) {
#if defined(__BMI2__) && defined(__x86_64__)
  return _pext_u64(value, selected);
#else
  Mask result = 0, dst = 1;
  for (; selected; selected &= selected-1, dst <<= 1)
    if (value & (selected & -selected)) result |= dst;
  return result;
#endif
}

struct BucketMessage {
  Mask scope;
  std::vector<Mask> counts;
};

// Eliminate one variable per bucket. A message counts assignments to already
// eliminated variables for each assignment to its remaining neighbors. Each
// eliminated variable belongs to just one incoming subtree, so products and
// partial sums fit in 64 bits for at most 62 original variables.
inline std::optional<Mask> bucket_count(unsigned n, const std::vector<Clause>& clauses,
                                        unsigned state_limit) {
  if (n > 62 || !state_limit) return {};
  constexpr unsigned max_bucket_width = 15;
  const unsigned max_width = std::min<unsigned>(max_bucket_width, std::bit_width(state_limit)-1);
  std::array<Mask, 62> adjacency{}, later{}, rank_bit{};
  for (auto c : clauses) {
    Mask vars = c.vars();
    if (!vars) return 0;
    if (c.pos & c.neg) continue;
    // A clause is a clique in the primal graph, independent of order.
    if (unsigned(std::popcount(vars)) > max_width+1) return {};
    for (Mask vs = vars; vs; vs &= vs-1) adjacency[std::countr_zero(vs)] |= vars;
  }
  Mask remaining = bit(n)-1;
  uint64_t work = 0;
  for (unsigned step = 0; step < n; ++step) {
    unsigned best = 0, degree = 64;
    for (Mask vs = remaining; vs; vs &= vs-1) {
      unsigned v = std::countr_zero(vs);
      unsigned d = std::popcount(adjacency[v] & remaining & ~bit(v));
      if (d < degree) { best = v; degree = d; }
    }
    if (degree > max_width) return {};
    work += bit(degree+1);
    if (work > uint64_t(state_limit)*16) return {};
    remaining ^= bit(best);
    rank_bit[best] = bit(step);
    Mask neighbors = later[step] = adjacency[best] & remaining;
    // Summing out best connects all of its remaining neighbors.
    for (Mask vs = neighbors; vs; vs &= vs-1)
      adjacency[std::countr_zero(vs)] |= neighbors;
  }
  auto permute = [&](Mask vars) {
    Mask mask = 0;
    for (; vars; vars &= vars-1) mask |= rank_bit[std::countr_zero(vars)];
    return mask;
  };
  // In rank order the current variable is always bit zero of a packed table.
  std::array<std::vector<Clause>, 62> local;
  std::array<std::vector<BucketMessage>, 62> incoming;
  for (auto c : clauses) {
    if (c.pos & c.neg) continue;
    Clause ranked{permute(c.pos), permute(c.neg)};
    local[std::countr_zero(ranked.vars())].push_back(ranked);
  }
  Mask answer = 1;
  for (unsigned step = 0; step < n; ++step) {
    Mask rest = permute(later[step]), bag = rest | bit(step);
    const unsigned size = 1u << std::popcount(bag);
    std::vector<Mask> table(size, 1);
    for (auto c : local[step]) {
      unsigned fixed = extract_bits(c.neg, bag);
      unsigned free = (size-1) ^ extract_bits(c.vars(), bag);
      for (unsigned subset = free;; subset = (subset-1) & free) {
        table[fixed | subset] = 0;
        if (!subset) break;
      }
    }
    for (const auto& msg : incoming[step]) {
      Mask selector = extract_bits(msg.scope, bag);
#if defined(__BMI2__) && defined(__x86_64__)
      for (unsigned i = 0; i < size; ++i)
        table[i] *= msg.counts[extract_bits(i, selector)];
#else
      // Incrementing i flips its low bits through the first carry. In message
      // coordinates these are precisely the low selected bits of that prefix.
      std::array<unsigned, max_bucket_width + 1> carry_mask;
      unsigned selected = 0;
      for (unsigned j = 0; (1u << j) < size; ++j) {
        selected += bool(selector & bit(j));
        carry_mask[j] = (1u << selected) - 1;
      }
      unsigned index = 0;
      table[0] *= msg.counts[0];
      for (unsigned i = 1; i < size; ++i) {
        index ^= carry_mask[std::countr_zero(i)];
        table[i] *= msg.counts[index];
      }
#endif
    }
    incoming[step].clear();
    for (unsigned i = 0; i < size/2; ++i) table[i] = table[2*i]+table[2*i+1];
    table.resize(size/2);
    if (rest) incoming[std::countr_zero(rest)].push_back({rest, std::move(table)});
    else answer *= table.front();
    if (!answer) break;
  }
  return answer;
}
}
