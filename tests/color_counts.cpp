#include "src/color_sweep.hpp"
#include "src/lit.hpp"
#include <random>
#include <stdexcept>

using namespace GanakInt;
using CNF = std::vector<std::vector<Lit>>;

static void check(unsigned n, const CNF& cnf, bool must_decode) {
    auto answer = color_sweep::solve(n, cnf);
    if (!answer && !must_decode) return;
    unsigned total = 0;
    for (unsigned assignment=0; assignment<(1u<<n); ++assignment) {
        bool satisfied=true;
        for (const auto& clause:cnf) {
            bool hit=false;
            for (auto l:clause) hit |= ((assignment>>l.var())&1u) != l.sign();
            if (!hit) { satisfied=false; break; }
        }
        total+=satisfied;
    }
    if (!answer || *answer!=total) throw std::runtime_error("incorrect color count");
}

int main() {
    std::mt19937 rng(2898948);
    for (unsigned t=0;t<160;++t) {
        bool onehot=t%2;
        unsigned bits=(t%3==1 ? 3 : 2), k=onehot ? 3 : (1u<<bits), m=onehot ? 4 : 3;
        unsigned n=m*(onehot ? k : bits);
        CNF cnf;
        unsigned palette=onehot ? (1u<<k)-1 : rng()%(1u<<k);
        for (unsigned v=0;v<m;++v) {
            if (onehot) {
                std::vector<Lit> choices;
                for (unsigned a=0;a<k;++a) {
                    choices.emplace_back(v*k+a,false);
                    for (unsigned b=0;b<a;++b)
                        cnf.push_back({Lit(v*k+a,true),Lit(v*k+b,true)});
                }
                cnf.push_back(choices);
            } else for (unsigned a=0;a<k;++a) if (!(palette&(1u<<a))) {
                std::vector<Lit> cl;
                for (unsigned j=0;j<bits;++j) cl.emplace_back(v+m*j,(a>>j)&1);
                cnf.push_back(cl);
            }
            for (unsigned u=0;u<v;++u) if ((u==0 && v==1) || rng()%2) {
                for (unsigned a=0;a<k;++a) if (palette&(1u<<a)) {
                    if (onehot) cnf.push_back({Lit(u*k+a,true),Lit(v*k+a,true)});
                    else {
                        std::vector<Lit> cl;
                        for (unsigned w:{u,v}) for (unsigned j=0;j<bits;++j)
                            cl.emplace_back(w+m*j,(a>>j)&1);
                        cnf.push_back(cl);
                    }
                }
            }
        }
        std::shuffle(cnf.begin(),cnf.end(),rng);
        check(n,cnf,true);
        // Incomplete edge gadgets and vertex-dependent palettes must fall
        // through, unless the altered formula still has a valid certificate.
        cnf.erase(cnf.begin()+rng()%cnf.size());
        check(n,cnf,false);
        cnf.push_back({Lit(rng()%n,false)});
        check(n,cnf,false);
    }
    // Verify the cached path separately, including isolated vertices,
    // disconnected components, conflicts and repeated canonical states.
    for (unsigned t=0;t<100;++t) {
        color_sweep::Graph g(7);
        for(unsigned v=0;v<7;++v) for(unsigned u=0;u<v;++u) if(rng()%3==0) {
            g[u].push_back(v);g[v].push_back(u);
        }
        unsigned ref=0;
        for(unsigned a=0;a<(1u<<14);++a) {
            bool ok=true;
            for(unsigned v=0;v<7;++v) for(auto u:g[v])
                if (((a>>(2*u))&3)==((a>>(2*v))&3)) ok=false;
            ref+=ok;
        }
        if (colorcache::Solver<2>(g,100000).solve()!=mpz_class(ref))
            throw std::runtime_error("incorrect cached color count");
    }
    for (unsigned n : {63u,64u,65u,127u,128u}) {
        color_sweep::Graph tree(n);
        for (unsigned v=1;v<n;++v) {
            unsigned parent = n<=65 ? v-1 : 0;
            tree[parent].push_back(v);
            tree[v].push_back(parent);
        }
        mpz_class expected;
        mpz_ui_pow_ui(expected.get_mpz_t(),3,n-1);
        expected*=4;
        if (colorcache::Solver<2>(tree,100000).solve()!=expected)
            throw std::runtime_error("incorrect tree count at mask boundary");
        if (colorcache::Solver<2>(tree,1).solve().has_value())
            throw std::runtime_error("partial search must decline");
    }
}
