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

// Counts proper k-colorings of a graph recognized as the WHOLE input CNF
// (by Christian Reitwiessner, PR #129). The decoders only accept a CNF whose
// models are in bijection with the proper colorings; anything else declines.
#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <optional>
#include <vector>
#include <gmpxx.h>

namespace GanakInt::color {

using Graph = std::vector<std::vector<uint32_t>>;
struct Coloring { Graph edges; uint32_t colors; };
constexpr uint32_t max_vertices = 200;

// Bit-plane encoding: vertex v is the binary code of vars v, v+m, v+2m, ...
// Clauses of size `bits` on one vertex forbid a code (must be the same set
// for every vertex); clauses of size 2*bits forbid equal codes on an edge.
template<class Clauses>
std::optional<Coloring> decode_bits(uint32_t n, const Clauses& clauses, uint32_t bits) {
  if (n % bits || n / bits > max_vertices || n / bits < 2) return {};
  const uint32_t m = n / bits, colors = 1u << bits;
  std::vector<uint32_t> forbidden(m*m), unavailable(m);
  for (const auto& cl : clauses) {
    if (cl.size() == bits) {
      uint32_t vertex = cl[0].var()%m, planes = 0, code = 0;
      for (const auto& l : cl) {
        uint32_t plane = l.var()/m;
        if (l.var()%m != vertex || (planes & (1u<<plane))) return {};
        planes |= 1u<<plane;
        if (l.sign()) code |= 1u<<plane;
      }
      unavailable[vertex] |= 1u<<code;
      continue;
    }
    if (cl.size() != 2*bits) return {};
    uint32_t u = cl[0].var() % m, v = m;
    uint32_t seen[2]{}, code[2]{};
    for (const auto& l : cl) {
      uint32_t vertex = l.var() % m, plane = l.var() / m;
      if (vertex != u && v == m) v = vertex;
      if (vertex != u && vertex != v) return {};
      uint32_t side = vertex != u;
      if (seen[side] & (1u << plane)) return {};
      seen[side] |= 1u << plane;
      if (l.sign()) code[side] |= 1u << plane;
    }
    if (v == m || code[0] != code[1]) return {};
    forbidden[std::min(u,v)*m + std::max(u,v)] |= 1u << code[0];
  }
  // A per-vertex palette would break the color symmetry the counters rely on
  if (!std::all_of(unavailable.begin(), unavailable.end(),
        [&](uint32_t mask){ return mask == unavailable.front(); })) return {};
  const uint32_t palette = ((1u<<colors)-1) & ~unavailable.front();
  Graph edges(m);
  for (uint32_t u = 0; u < m; u++) for (uint32_t v = u+1; v < m; v++) {
    uint32_t mask = forbidden[u*m+v] & palette;
    if (!mask) continue;
    if (mask != palette) return {};
    edges[u].push_back(v);
    edges[v].push_back(u);
  }
  return Coloring{std::move(edges), (uint32_t)std::popcount(palette)};
}

// One-hot encoding: the positive clauses (all of size k) partition the vars
// into vertices, every pair inside a vertex is excluded, and an edge excludes
// every color at the same clause position on both endpoints.
template<class Clauses>
std::optional<Coloring> decode_onehot(uint32_t n, const Clauses& clauses) {
  std::vector<int> block(n, -1);
  std::vector<uint32_t> color(n);
  uint32_t k = 0, m = 0;
  for (const auto& cl : clauses) {
    if (cl.size() == 2 && cl[0].sign() && cl[1].sign()) continue;
    if (cl.size() < 3 || cl.size() > 8 || (k && k != cl.size()) || m >= max_vertices) return {};
    k = cl.size();
    for (uint32_t j = 0; j < k; j++) {
      const auto l = cl[j];
      if (l.sign() || block[l.var()] != -1) return {};
      block[l.var()] = m;
      color[l.var()] = j;
    }
    m++;
  }
  if (!k || m*k != n) return {};
  std::vector<uint64_t> within(m);
  std::vector<uint32_t> cross(m*m);
  uint64_t required = 0;
  for (uint32_t j = 0; j < k; j++) for (uint32_t i = 0; i < j; i++)
    required |= uint64_t{1}<<(j*8+i);
  for (const auto& cl : clauses) {
    if (cl.size() != 2) continue;
    uint32_t a = cl[0].var(), b = cl[1].var();
    uint32_t u = block[a], v = block[b], i = color[a], j = color[b];
    if (u == v) {
      if (i == j) return {};
      within[u] |= uint64_t{1}<<(8*std::max(i,j)+std::min(i,j));
    } else {
      if (i != j) return {};
      cross[std::min(u,v)*m+std::max(u,v)] |= 1u<<i;
    }
  }
  Graph edges(m);
  for (uint32_t u = 0; u < m; u++) {
    if (within[u] != required) return {};
    for (uint32_t v = u+1; v < m; v++) {
      if (!cross[u*m+v]) continue;
      if (cross[u*m+v] != (1u<<k)-1) return {};
      edges[u].push_back(v);
      edges[v].push_back(u);
    }
  }
  return Coloring{std::move(edges), k};
}

template<class Clauses>
std::optional<Coloring> decode(uint32_t n, const Clauses& clauses) {
  if (n < 3 || n > 1600 || clauses.empty() || clauses.size() > 64000) return {};
  std::optional<Coloring> g;
  const uint32_t width = clauses.front().size();
  if (width == 2 || width == 3) g = decode_bits(n, clauses, width);
  if (width == 4 || width == 6) g = decode_bits(n, clauses, width/2);
  if (!g) g = decode_onehot(n, clauses);
  return g;
}

// DSATUR-style search; colors never used so far are interchangeable, so
// only the first unused one is tried and weighted by the number unused.
class Search {
  const Coloring& g;
  std::vector<uint32_t> forbidden;
  std::vector<bool> assigned;
  uint64_t remaining_nodes;
public:
  Search(const Coloring& graph, uint64_t budget) :
    g(graph), forbidden(g.edges.size()), assigned(g.edges.size()), remaining_nodes(budget) {}
  // A decline leaves the object in a dirty state: throw it away
  std::optional<mpz_class> count(uint32_t left, uint32_t used) {
    if (!remaining_nodes) return {};
    remaining_nodes--;
    if (!left) return mpz_class(1);
    uint32_t best = g.edges.size(), fewest = g.colors+1;
    for (uint32_t v = 0; v < g.edges.size(); v++) if (!assigned[v]) {
      uint32_t choices = g.colors-std::popcount(forbidden[v]);
      if (!choices) return mpz_class(0);
      if (choices < fewest || (choices == fewest && g.edges[v].size() > g.edges[best].size())) {
        best = v;
        fewest = choices;
      }
    }
    assigned[best] = true;
    mpz_class total = 0;
    for (uint32_t i = 0; i < std::min(used+1, g.colors); i++) {
      const uint32_t bit = 1u<<i;
      if (forbidden[best]&bit) continue;
      uint32_t changed[max_vertices], size = 0;
      for (auto v : g.edges[best]) if (!assigned[v] && !(forbidden[v]&bit)) {
        changed[size++] = v;
        forbidden[v] |= bit;
      }
      auto value = count(left-1, used+(i == used));
      if (!value) return {};
      for (uint32_t j = 0; j < size; j++) forbidden[changed[j]] ^= bit;
      total += *value * (i == used ? g.colors-used : 1u);
    }
    assigned[best] = false;
    return total;
  }
};

template<size_t N> struct Mask {
  std::array<uint64_t, N> words{};
  auto operator<=>(const Mask&) const = default;
  static Mask singleton(uint32_t v) {
    Mask r;
    r.words[v/64] = uint64_t(1) << (v%64);
    return r;
  }
  explicit operator bool() const {
    for (auto w : words) if (w) return true;
    return false;
  }
  uint32_t count() const {
    uint32_t sum = 0;
    for (auto w : words) sum += std::popcount(w);
    return sum;
  }
  uint32_t first() const {
    for (uint32_t i = 0; i < N; i++) if (words[i]) return 64*i+std::countr_zero(words[i]);
    return uint32_t(-1);
  }
  Mask operator|(const Mask& b) const { Mask c; for (size_t i = 0; i < N; i++) c.words[i] = words[i]|b.words[i]; return c; }
  Mask operator&(const Mask& b) const { Mask c; for (size_t i = 0; i < N; i++) c.words[i] = words[i]&b.words[i]; return c; }
  Mask without(const Mask& b) const { Mask c; for (size_t i = 0; i < N; i++) c.words[i] = words[i]&~b.words[i]; return c; }
};

inline uint64_t mix(uint64_t x) {
  x ^= x>>30; x *= 0xbf58476d1ce4e5b9ULL; x ^= x>>27; x *= 0x94d049bb133111ebULL;
  return x^(x>>31);
}

// State d[c] = vertices that may still take color c; colored vertices are in
// none. The count only depends on the multiset of the d[c]s (colors are
// symmetric), so sorted d is the cache key, and equal d[c]s branch once.
template<size_t N, uint32_t K> class DomainSearch {
  static_assert(K == 4 || K == 5);
  using M = Mask<N>;
  using State = std::array<M, K>;
  struct Entry { State key{}; mpz_class value; bool active = false; };
  std::vector<M> neighbors;
  std::vector<Entry> cache;
  uint64_t limit;
  uint64_t nodes = 0;
  bool exhausted = false;

  static M vertices(const State& d) {
    M all = (d[0]|d[1]) | (d[2]|d[3]);
    if constexpr (K == 5) all = all|d[4];
    return all;
  }
  // Only neighbors sharing a possible color constrain one another
  M component(const State& d, M active) const {
    M reached, front = M::singleton(active.first());
    while (front) {
      reached = reached|front;
      M next;
      while (front) {
        uint32_t v = front.first(); M bit = M::singleton(v);
        front = front.without(bit);
        M colors;
        for (const M& domain : d) if (domain&bit) colors = colors|domain;
        next = next|(neighbors[v]&colors);
      }
      front = next.without(reached);
    }
    return reached;
  }
  mpz_class search(State d, uint32_t depth) {
    if (nodes == limit) { exhausted = true; return 0; }
    nodes++;
    M active = vertices(d), choice;
    while (active) {
      M a = d[0]&d[1], b = d[2]&d[3], x = d[0]|d[1], y = d[2]|d[3];
      M two = a|b|(x&y), three = (a&y)|(b&x), four = a&b, five;
      if constexpr (K == 5) {
        five = four&d[4];
        four = four|(three&d[4]);
        three = three|(two&d[4]);
        two = two|((x|y)&d[4]);
      }
      M forced = active.without(two);
      if (!forced) {
        choice = active.without(three);
        if (!choice) choice = active.without(four);
        if constexpr (K == 5) { if (!choice) choice = active.without(five); }
        if (!choice) choice = active;
        break;
      }
      uint32_t v = forced.first(); M bit = M::singleton(v);
      for (M& mask : d) if (mask&bit) { mask = mask.without(neighbors[v]|bit); break; }
      active = active.without(bit);
      if (vertices(d) != active) return 0;
    }
    if (!active) return 1;
    std::sort(d.begin(), d.end());
    const uint32_t remaining = active.count();
    const bool memoize = remaining >= 5 && remaining <= 52;
    size_t slot = 0;
    if (memoize) {
      uint64_t hash = 0;
      for (const M& domain : d) for (uint64_t w : domain.words) hash = mix(hash^w);
      slot = hash&(cache.size()-1);
      if (cache[slot].active && cache[slot].key == d) return cache[slot].value;
    }
    if (remaining <= 36 && depth%4 == 0) {
      M part = component(d, active);
      if (part != active) {
        State left = d, right = d;
        for (uint32_t c = 0; c < K; c++) { left[c] = left[c]&part; right[c] = right[c].without(part); }
        mpz_class result = search(left, depth+1);
        if (!exhausted && result != 0) result *= search(right, depth+1);
        return result;
      }
    }
    uint32_t v = choice.first(); int best = -1;
    while (choice) {
      uint32_t u = choice.first();
      choice = choice.without(M::singleton(u));
      int degree = (neighbors[u]&active).count();
      if (degree > best) { v = u; best = degree; }
    }
    M bit = M::singleton(v);
    State base = d;
    for (M& m : base) m = m.without(bit);
    if (best == 0) {
      uint32_t num = 0;
      for (const M& mask : d) num += bool(mask&bit);
      return search(base, depth+1)*num;
    }
    mpz_class total = 0;
    for (uint32_t c = 0; c < K;) {
      uint32_t next = c+1;
      while (next < K && d[next] == d[c]) next++;
      if (d[c]&bit) {
        auto child = base;
        child[c] = child[c].without(neighbors[v]);
        if (vertices(child) == active.without(bit)) total += search(child, depth+1)*(next-c);
        if (exhausted) return 0;
      }
      c = next;
    }
    if (memoize) cache[slot] = {d, total, true};
    return total;
  }

public:
  DomainSearch(const Graph& g, uint64_t budget) :
      neighbors(g.size()), cache(size_t(1) << std::min<uint32_t>(N == 2 ? 19 : 18, g.size())), limit(budget) {
    for (uint32_t v = 0; v < g.size(); v++)
      for (uint32_t u : g[v]) neighbors[v] = neighbors[v]|M::singleton(u);
  }
  std::optional<mpz_class> solve() {
    M all;
    for (uint32_t i = 0; i < neighbors.size(); i++) all = all|M::singleton(i);
    State initial; initial.fill(all);
    auto value = search(initial, 0);
    if (exhausted) return {};
    return value;
  }
};

// domain_budget 0 = unlimited
inline std::optional<mpz_class> count(const Coloring& g, uint64_t search_budget, uint64_t domain_budget) {
  if (auto v = Search(g, search_budget).count(g.edges.size(), 0)) return v;
  const uint64_t lim = domain_budget ? domain_budget : UINT64_MAX;
  const size_t n = g.edges.size();
  if (g.colors == 4 && n <= 128) return DomainSearch<2, 4>(g.edges, lim).solve();
  if (g.colors == 4 && n <= 192) return DomainSearch<3, 4>(g.edges, lim).solve();
  if (g.colors == 5 && n <= 128) return DomainSearch<2, 5>(g.edges, lim).solve();
  if (g.colors == 5 && n <= 192) return DomainSearch<3, 5>(g.edges, lim).solve();
  return {};
}

}
