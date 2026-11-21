// Certified ordinary graph-coloring encodings. Each graph coloring extends to
// exactly one input assignment. Declining at any stage leaves the CNF intact.
#pragma once
#include "color_domain_cache.hpp"
#include <algorithm>
#include <bit>
#include <cstdint>
#include <optional>
#include <vector>
#include <gmpxx.h>

namespace GanakInt::color_sweep {
using Graph = std::vector<std::vector<unsigned>>;
struct Coloring { Graph edges; unsigned colors; };

// Bit-plane encoding: vertex v uses v, v+m, ... . An edge forbids all
// assignments where the two endpoints carry the same binary code.
template<class Clauses>
std::optional<Coloring> decode_bits(unsigned n, const Clauses& clauses, unsigned bits) {
    if (n % bits || n / bits > 200 || n / bits < 2) return {};
    const unsigned m = n / bits, colors = 1u << bits;
    std::vector<unsigned> forbidden(m*m), unavailable(m);
    for (const auto& cl : clauses) {
        if (cl.size() == bits) {
            unsigned vertex=cl[0].var()%m, planes=0, code=0;
            for (auto l:cl) {
                unsigned plane=l.var()/m;
                if (l.var()%m!=vertex || (planes & (1u<<plane))) return {};
                planes |= 1u<<plane;
                if (l.sign()) code |= 1u<<plane;
            }
            unavailable[vertex] |= 1u<<code;
            continue;
        }
        if (cl.size() != 2*bits) return {};
        unsigned u = cl[0].var() % m, v = m;
        unsigned seen[2]{}, code[2]{};
        for (auto l : cl) {
            unsigned vertex = l.var() % m, plane = l.var() / m;
            if (vertex != u && v == m) v = vertex;
            if (vertex != u && vertex != v) return {};
            unsigned side = vertex != u;
            if (seen[side] & (1u << plane)) return {};
            seen[side] |= 1u << plane;
            if (l.sign()) code[side] |= 1u << plane;
        }
        if (v == m || code[0] != code[1]) return {};
        forbidden[std::min(u,v)*m + std::max(u,v)] |= 1u << code[0];
    }
    // A global palette restriction preserves color permutation symmetry.
    if (!std::all_of(unavailable.begin(),unavailable.end(),
                     [&](unsigned mask){return mask==unavailable.front();})) return {};
    unsigned palette=((1u<<colors)-1) & ~unavailable.front();
    Graph edges(m);
    for (unsigned u=0; u<m; ++u) for (unsigned v=u+1; v<m; ++v) {
        unsigned mask=forbidden[u*m+v] & palette;
        if (!mask) continue;
        if (mask != palette) return {};
        edges[u].push_back(v);
        edges[v].push_back(u);
    }
    return Coloring{std::move(edges), unsigned(std::popcount(palette))};
}

// Positive clauses partition the variables. Require every within-block
// pair exclusion and every color of every observed across-block inequality.
template<class Clauses>
std::optional<Coloring> decode_onehot(unsigned n, const Clauses& clauses) {
    std::vector<int> block(n,-1);
    std::vector<unsigned> color(n);
    unsigned k=0,m=0;
    for (const auto& cl : clauses) {
        if (cl.size()==2 && cl[0].sign() && cl[1].sign()) continue;
        if (cl.size()<3 || cl.size()>8 || (k && k!=cl.size()) || m>=200) return {};
        k=cl.size();
        for (unsigned j=0; j<k; ++j) {
            auto l=cl[j];
            if (l.sign() || block[l.var()]!=-1) return {};
            block[l.var()]=m;
            color[l.var()]=j;
        }
        ++m;
    }
    if (!k || m*k!=n) return {};
    std::vector<uint64_t> within(m);
    std::vector<unsigned> cross(m*m);
    uint64_t required=0;
    for (unsigned j=0; j<k; ++j) for (unsigned i=0;i<j;++i)
        required |= uint64_t{1}<<(j*8+i);
    for (const auto& cl : clauses) {
        if (cl.size()!=2) continue;
        unsigned a=cl[0].var(), b=cl[1].var();
        unsigned u=block[a],v=block[b],i=color[a],j=color[b];
        if (u==v) {
            if (i==j) return {};
            within[u] |= uint64_t{1}<<(8*std::max(i,j)+std::min(i,j));
        } else {
            if (i!=j) return {};
            cross[std::min(u,v)*m+std::max(u,v)] |= 1u<<i;
        }
    }
    Graph edges(m);
    for (unsigned u=0; u<m; ++u) {
        if (within[u]!=required) return {};
        for (unsigned v=u+1;v<m;++v) {
            if (!cross[u*m+v]) continue;
            if (cross[u*m+v]!=(1u<<k)-1) return {};
            edges[u].push_back(v);
            edges[v].push_back(u);
        }
    }
    return Coloring{std::move(edges),k};
}

class Search {
    const Coloring& g;
    std::vector<unsigned> forbidden;
    std::vector<bool> assigned;
    unsigned remaining_nodes=250000;
public:
    explicit Search(const Coloring& graph):g(graph), forbidden(g.edges.size()),assigned(g.edges.size()){}
    std::optional<mpz_class> count(unsigned left, unsigned used) {
        if (!remaining_nodes) return {};
        --remaining_nodes;
        if (!left) return mpz_class(1);
        unsigned best=g.edges.size(), fewest=g.colors+1;
        for (unsigned v=0; v<g.edges.size(); ++v) if (!assigned[v]) {
            unsigned choices=g.colors-std::popcount(forbidden[v]);
            if (!choices) return mpz_class(0);
            if (choices<fewest || (choices==fewest && g.edges[v].size()>g.edges[best].size())) {
                best=v;
                fewest=choices;
            }
        }
        assigned[best]=true;
        mpz_class total=0;
        for (unsigned i=0; i<std::min(used+1,g.colors); ++i) {
            const unsigned bit=1u<<i;
            if (forbidden[best]&bit) continue;
            unsigned changed[200], size=0;
            for (auto v:g.edges[best]) if (!assigned[v] && !(forbidden[v]&bit)) {
                changed[size++]=v;
                forbidden[v]|=bit;
            }
            auto value=count(left-1,used+(i==used));
            // A decline discards the entire Search object.
            if (!value) return {};
            for (unsigned j=0;j<size;++j) forbidden[changed[j]]^=bit;
            // Unused labels are interchangeable under every remaining constraint.
            total += *value * (i==used ? g.colors-used : 1u);
        }
        assigned[best]=false;
        return total;
    }
};

template<class Clauses>
std::optional<mpz_class> solve(unsigned n,const Clauses& clauses) {
    if (n<3 || n>1600 || clauses.empty() || clauses.size()>64000) return {};
    std::optional<Coloring> graph;
    unsigned width=clauses.front().size();
    if (width==2 || width==3) graph=decode_bits(n,clauses,width);
    if (width==4 || width==6) graph=decode_bits(n,clauses,width/2);
    if (!graph) graph=decode_onehot(n,clauses);
    if (!graph) return {};
    auto value=Search(*graph).count(graph->edges.size(),0);
    if (value || graph->colors!=4 || graph->edges.size()>128) return value;
    // Store normalized residual color domains only after the cheap search
    // declines. A capped search discards partial counts and falls through.
    return colorcache::Solver<2>(graph->edges,180000000).solve();
}
}
