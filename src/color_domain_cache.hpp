// Four masks describe the residual CSP: mask c contains the vertices that
// may still take color c. Removed vertices have already been counted.
// Sorting masks gives a cache key invariant under color permutations.
#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <optional>
#include <vector>
#include <gmpxx.h>

namespace GanakInt::colorcache {
template<size_t N> struct Mask {
  std::array<uint64_t,N> words{};
  auto operator<=>(const Mask&) const = default;
  static Mask singleton(unsigned v) {
    Mask result;
    result.words[v/64] = uint64_t(1) << (v%64);
    return result;
  }
  explicit operator bool() const {
    for (auto w: words) if (w) return true;
    return false;
  }
  unsigned count() const {
    unsigned sum=0;
    for (auto w: words) sum+=std::popcount(w);
    return sum;
  }
  unsigned first() const {
    for (unsigned i=0; i<N; ++i) if (words[i]) return 64*i+std::countr_zero(words[i]);
    return unsigned(-1);
  }
  Mask operator|(const Mask& b) const { Mask c; for(size_t i=0;i<N;++i)c.words[i]=words[i]|b.words[i]; return c; }
  Mask operator&(const Mask& b) const { Mask c; for(size_t i=0;i<N;++i)c.words[i]=words[i]&b.words[i]; return c; }
  Mask without(const Mask& b) const { Mask c; for(size_t i=0;i<N;++i)c.words[i]=words[i]&~b.words[i]; return c; }
};
inline uint64_t mix(uint64_t x) {
  x^=x>>30; x*=0xbf58476d1ce4e5b9ULL; x^=x>>27; x*=0x94d049bb133111ebULL;
  return x^(x>>31);
}

template<size_t N> class Solver {
  using M=Mask<N>;
  using State=std::array<M,4>;
  struct Entry { State key{}; mpz_class value; bool active=false; };
  std::vector<M> neighbors;
  std::vector<Entry> cache;
  uint64_t limit;
  static M vertices(const State& d) { return (d[0]|d[1]) | (d[2]|d[3]); }
  // Only neighbors sharing a possible color can constrain one another.
  M component(const State& d, M active) const {
    M reached, front=M::singleton(active.first());
    while (front) {
      reached=reached|front;
      M next;
      while (front) {
        unsigned v=front.first(); M bit=M::singleton(v);
        front=front.without(bit);
        M colors;
        for (const M& domain : d) if (domain&bit) colors=colors|domain;
        next=next|(neighbors[v]&colors);
      }
      front=next.without(reached);
    }
    return reached;
  }
  mpz_class search(State d, unsigned depth) {
    if (++nodes>limit) { exhausted=true; return 0; }
    M active=vertices(d), choice;
    while (active) {
      M a=d[0]&d[1], b=d[2]&d[3], x=d[0]|d[1], y=d[2]|d[3];
      M forced=active.without(a|b|(x&y));
      if (!forced) {
        choice=active.without((a&y)|(b&x));
        if (!choice) choice=active.without(a&b);
        if (!choice) choice=active;
        break;
      }
      unsigned v=forced.first(); M bit=M::singleton(v);
      for (M& mask : d) if (mask&bit) { mask=mask.without(neighbors[v]|bit); break; }
      active=active.without(bit);
      if (vertices(d)!=active) return 0;
    }
    if (!active) return 1;
    std::sort(d.begin(),d.end());
    unsigned remaining=active.count();
    const bool memoize=remaining>=5 && remaining<=52;
    size_t slot=0;
    if (memoize) {
      uint64_t hash=0;
      for (const M& domain : d) for (uint64_t w : domain.words) hash=mix(hash^w);
      slot=hash&(cache.size()-1);
      if (cache[slot].active && cache[slot].key==d) { ++hits; return cache[slot].value; }
    }
    if (remaining<=36 && depth%4==0) {
      M part=component(d,active);
      if (part!=active) {
        State left=d,right=d;
        for (unsigned c=0;c<4;++c) { left[c]=left[c]&part; right[c]=right[c].without(part); }
        mpz_class result=search(left,depth+1);
        if (!exhausted && result!=0) result*=search(right,depth+1);
        return result;
      }
    }
    unsigned v=choice.first(); int best=-1;
    while (choice) {
      unsigned u=choice.first();
      choice=choice.without(M::singleton(u));
      int degree=(neighbors[u]&active).count();
      if (degree>best) { v=u; best=degree; }
    }
    M bit=M::singleton(v);
    State base=d;
    for (M& m : base) m=m.without(bit);
    if (best==0) {
      unsigned num=0;
      for (const M& mask : d) num+=bool(mask&bit);
      mpz_class result=search(base,depth+1);
      return result*num;
    }
    mpz_class total=0;
    for (unsigned c=0;c<4;) {
      unsigned next=c+1;
      while (next<4 && d[next]==d[c]) ++next;
      if (d[c]&bit) {
        auto child=base;
        child[c]=child[c].without(neighbors[v]);
        if (vertices(child)==active.without(bit)) total+=search(child,depth+1)*(next-c);
        if (exhausted) return 0;
      }
      c=next;
    }
    if (memoize) cache[slot]={d,total,true};
    return total;
  }
 public:
  bool exhausted=false;
  uint64_t nodes=0,hits=0;
  Solver(const std::vector<std::vector<unsigned>>& graph,uint64_t budget):
      neighbors(graph.size()),cache(size_t(1)<<std::min<unsigned>(N==2?19:18,graph.size())),limit(budget) {
    for (unsigned v=0;v<graph.size();++v)
      for (unsigned u:graph[v]) neighbors[v]=neighbors[v]|M::singleton(u);
  }
  std::optional<mpz_class> solve() {
    M all;
    for (unsigned i=0;i<neighbors.size();++i) all=all|M::singleton(i);
    auto value=search({all,all,all,all},0);
    if (exhausted) return {};
    return value;
  }
};
}
