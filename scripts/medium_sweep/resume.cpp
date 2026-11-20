#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
#include <gmpxx.h>

using Mask = uint64_t;
using Domains = std::vector<Mask>;
struct Edge { unsigned to; std::vector<Mask> support; };
struct Exhausted { mpz_class done=0; std::vector<std::pair<unsigned,unsigned>> path; };
struct Solver {
    std::vector<std::vector<Edge>> edges;
    uint64_t nodes=0,limit=UINT64_MAX;
    bool mac, adaptive=false;
    std::vector<uint64_t> weighted_degree;
    bool propagate(Domains& d,std::vector<unsigned> queue) {
        std::vector<bool> pending(d.size());
        for (auto v:queue) pending[v]=true;
        for(size_t i=0;i<queue.size();++i){
            unsigned v=queue[i]; pending[v]=false;
            for(const Edge& e:edges[v]){
                Mask allowed=0, target=d[e.to];
                for(Mask a=d[v];a && (allowed&target)!=target;a&=a-1)
                    allowed|=e.support[std::countr_zero(a)];
                Mask next=target&allowed;
                if(!next){
                    if(adaptive){++weighted_degree[v];++weighted_degree[e.to];}
                    return false;
                }
                if(next!=target){
                    d[e.to]=next;
                    if(!pending[e.to] && (mac || std::has_single_bit(next))){
                        queue.push_back(e.to);pending[e.to]=true;
                    }
                }
            }
        }
        return true;
    }
    mpz_class count(const Domains& d){
        if(nodes==limit)throw Exhausted{};
        ++nodes;
        unsigned best=d.size();int smallest=65;size_t degree=0;
        for(unsigned i=0;i<d.size();++i){
            int choices=std::popcount(d[i]);
            if(choices<=1)continue;
            uint64_t weight=adaptive?weighted_degree[i]:1;
            uint64_t old_weight=(adaptive && best<d.size())?weighted_degree[best]:1;
            if(best==d.size() || choices*old_weight<smallest*weight ||
               (choices*old_weight==smallest*weight && edges[i].size()>degree)){
                best=i;smallest=choices;degree=edges[i].size();
            }
        }
        if(best==d.size())return 1;
        mpz_class result=0;
        for(Mask a=d[best];a;a&=a-1){
            Domains next=d;next[best]=a&-a;
            if(propagate(next,{best})){
                try {result+=count(next);}
                catch(Exhausted& e){
                    e.done+=result;
                    e.path.emplace_back(best,std::countr_zero(a));
                    throw;
                }
            }
        }
        return result;
    }
};
int main(int argc,char**argv){
    unsigned n,m;if((argc!=2 && argc!=3) || !(std::cin>>n>>m) || !n || n>128)return 2;
    Solver s; if(argc==3)s.limit=std::stoull(argv[2]); s.mac=std::string(argv[1])=="mac";
    s.adaptive=std::string(argv[1])=="adaptive";
    s.edges.resize(n);Domains full(n);
    for(auto& mask:full){unsigned w;if(!(std::cin>>w) || w<2 || w>63)return 2;mask=(Mask(1)<<w)-1;}
    std::vector<std::vector<int>> slot(n,std::vector<int>(n,-1));
    auto forbid=[&](unsigned v,unsigned x,unsigned u,unsigned y){
        if(slot[v][u]==-1){
            slot[v][u]=s.edges[v].size();
            s.edges[v].push_back({u,std::vector<Mask>(std::popcount(full[v]),full[u])});
        }
        s.edges[v][slot[v][u]].support[x]&=~(Mask(1)<<y);
    };
    for(unsigned i=0;i<m;++i){
        unsigned v,x,u,y;if(!(std::cin>>v>>x>>u>>y) || v>=n || u>=n || v==u ||
            x>=unsigned(std::popcount(full[v])) || y>=unsigned(std::popcount(full[u])))return 2;
        forbid(v,x,u,y);forbid(u,y,v,x);
    }
    for(const auto& es:s.edges)s.weighted_degree.push_back(1+es.size());
    std::vector<unsigned> q;
    if(s.mac)for(unsigned i=0;i<n;++i)q.push_back(i);
    mpz_class result=0;
    unsigned jobs;
    if(std::cin>>jobs){
        while(jobs--){
            Domains current(n);std::vector<unsigned> changed;
            for(unsigned i=0;i<n;++i){
                if(!(std::cin>>current[i]) || !current[i] || (current[i]&~full[i]))return 2;
                if(current[i]!=full[i] || s.mac)changed.push_back(i);
            }
            if(s.propagate(current,changed))result+=s.count(current);
        }
        std::cout<<"exact "<<result<<" "<<s.nodes<<"\n";return 0;
    }
    try{if(s.propagate(full,q))result=s.count(full);}
    catch(const Exhausted& e){
        std::cout<<"partial "<<e.done<<" "<<s.nodes;
        for(auto it=e.path.rbegin();it!=e.path.rend();++it)std::cout<<" "<<it->first<<" "<<it->second;
        std::cout<<"\n";return 0;
    }
    std::cout<<"exact "<<result<<" "<<s.nodes<<"\n";
}
