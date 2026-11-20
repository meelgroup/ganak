#include "cadical.hpp"
#include <iostream>
#include <vector>
// Reads total variables, clause count, then zero-terminated clauses.
int main(){
    int n,m;
    if(!(std::cin>>n>>m)||n<0||m<0)return 2;
    CaDiCaL::Solver solver;
    solver.set("quiet",1);
    for(int i=0;i<m;++i){
        int x;
        do {if(!(std::cin>>x)||x>n||x< -n)return 2;solver.add(x);} while(x);
    }
    for(int count=0;count<=32;++count){
        int status=solver.solve();
        if(status==20){std::cout<<count<<"\n";return 0;}
        if(status!=10||count==32)break;
        std::vector<int> model;
        for(int v=1;v<=n;++v)model.push_back(solver.val(v)>0?v:-v);
        for(int x:model)solver.add(-x);
        solver.add(0);
    }
    std::cout<<"decline\n";
}
