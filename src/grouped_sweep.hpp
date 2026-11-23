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
// Dense variable-frontier sweep over truth tables of equal-scope clauses.
#pragma once
#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <type_traits>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <vector>
#include <gmpxx.h>

namespace GanakInt::grouped {
using Scope = std::vector<unsigned>;
struct Factor {
  Scope vars;
  std::vector<unsigned> forbidden;
};
class Factors {
  // Each clause forbids exactly one assignment to its distinct variables.
  std::map<Scope, unsigned> ids;

public:
  std::vector<Factor> factors;
  // Input literals use signed, one-based DIMACS notation.
  bool add(const std::vector<int> &lits) {
    std::map<unsigned, bool> signs;
    for (int l : lits) {
      unsigned v = std::abs(l);
      auto [it, fresh] = signs.emplace(v, l < 0);
      if (!fresh && it->second != (l < 0))
        return true; // tautology
    }
    if (signs.size() > 18)
      return false;
    Scope scope;
    unsigned forbidden = 0;
    for (auto [v, sign] : signs) {
      if (sign)
        forbidden |= 1u << scope.size();
      scope.push_back(v);
    }
    auto [it, fresh] = ids.emplace(scope, factors.size());
    if (fresh) {
      if (factors.size() >= 768)
        return false;
      factors.push_back({scope, {}});
    }
    factors[it->second].forbidden.push_back(forbidden);
    return true;
  }
};
struct Plan {
  std::vector<unsigned> order;
  uint64_t work = 0;
};
// Choose the cheapest of eight tie orders by total table entries touched.
// Peak width bounds scratch memory independently of the total-work limit.
inline std::optional<Plan> plan(const std::vector<Factor> &factors, unsigned n) {
  std::mt19937 rng(2024);
  std::optional<Plan> best;
  std::vector<unsigned> degree(n + 1);
  for (auto &f : factors)
    for (auto v : f.vars)
      ++degree[v];
  for (int start = 0; start < 8; ++start) {
    Plan trial;
    auto rem = degree;
    std::vector<bool> live(n + 1, false);
    std::vector<unsigned> todo(factors.size());
    std::iota(todo.begin(), todo.end(), 0);
    std::shuffle(todo.begin(), todo.end(), rng);
    int width = 0;
    while (!todo.empty()) {
      auto score = [&](unsigned id) {
        int merged = width, closes = 0;
        for (auto v : factors[id].vars) {
          merged += !live[v];
          closes += rem[v] == 1;
        }
        return std::pair{merged - closes, merged};
      };
      auto chosen = todo.begin();
      auto rank = score(*chosen);
      for (auto it = chosen + 1; it != todo.end(); ++it) {
        auto other = score(*it);
        if (other < rank) {
          chosen = it;
          rank = other;
        }
      }
      if (rank.second > (n > 63 ? 18 : 20))
        break;
      trial.work += uint64_t{1} << rank.second;
      if (trial.work > 8000000 || (best && trial.work >= best->work))
        break;
      for (auto v : factors[*chosen].vars)
        live[v] = --rem[v] > 0;
      width = rank.first;
      trial.order.push_back(*chosen);
      todo.erase(chosen);
    }
    if (todo.empty())
      best = std::move(trial);
  }
  return best;
}
// Incremental projection of an assignment index onto a subset of its bits.
inline std::vector<unsigned> projection(const Scope &from, const Scope &to) {
  std::vector<unsigned> toggles(from.size()), result(1u << from.size());
  for (unsigned i = 0; i < to.size(); ++i)
    toggles[std::lower_bound(from.begin(), from.end(), to[i]) - from.begin()] = 1u << i;
  for (unsigned a = 1; a < result.size(); ++a)
    result[a] = result[a & (a - 1)] | toggles[std::countr_zero(a)];
  return result;
}
// Sum out variables at their final occurrence. The next table has one
// exact count per assignment of the still-live variables.
template <class Count> mpz_class execute(const std::vector<Factor> &fs, const Plan &p, unsigned n) {
  std::vector<unsigned> rem(n + 1);
  for (auto &f : fs)
    for (auto v : f.vars)
      ++rem[v];
  unsigned free = std::count(rem.begin() + 1, rem.end(), 0);
  Scope live;
  std::vector<Count> count(1, 1);
  for (auto id : p.order) {
    auto &f = fs[id];
    Scope merged, next;
    std::set_union(live.begin(), live.end(), f.vars.begin(), f.vars.end(), std::back_inserter(merged));
    for (auto v : f.vars)
      --rem[v];
    for (auto v : merged)
      if (rem[v])
        next.push_back(v);
    auto old = projection(merged, live), factor = projection(merged, f.vars),
         dest = projection(merged, next);
    std::vector<bool> allowed(1u << f.vars.size(), true);
    for (auto a : f.forbidden)
      allowed[a] = false;
    std::vector<Count> out(1u << next.size());
    for (unsigned i = 0; i < old.size(); ++i)
      if (allowed[factor[i]])
        out[dest[i]] += count[old[i]];
    live.swap(next);
    count.swap(out);
  }
  mpz_class result;
  if constexpr (std::is_same_v<Count, uint64_t>)
    mpz_import(result.get_mpz_t(), 1, 1, sizeof(Count), 0, 0, count.data());
  else
    result = count[0];
  return result << free;
}
inline std::optional<mpz_class> solve(const Factors &groups, unsigned n) {
  auto order = plan(groups.factors, n);
  if (!order)
    return {};
  if (n <= 63)
    return execute<uint64_t>(groups.factors, *order, n);
  return execute<mpz_class>(groups.factors, *order, n);
}
} // namespace GanakInt::grouped
