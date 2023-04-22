#ifndef TREE_H
#define TREE_H

#include "array_id_func.h"
#include "tiny_id_func.h"
#include <cstdlib>
#include <utility>

template<class Neighbors>
bool is_tree(const Neighbors&neighbors){
	int node_count = neighbors.preimage_count();
	if(node_count <= 1)
		return true;

	int reachable_node_count = 0;

	ArrayIDFunc<int>parent(node_count);
	parent.fill(-1);
	ArrayIDFunc<int>stack(node_count);
	int stack_end = 0;
	stack[stack_end++] = 0;
	while(stack_end != 0){
		++reachable_node_count;
		auto x = stack[--stack_end];
		for(auto y:neighbors(x)){
			if(parent(x) != y){
				if(parent(y) == -1){
					parent[y] = x;
					stack[stack_end++] = y;
				}else{
					return false;
				}
			}
		}
	}
	return reachable_node_count == node_count;
}

#endif

