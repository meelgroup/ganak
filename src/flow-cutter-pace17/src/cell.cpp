#include "cell.h"

int get_treewidth_of_multilevel_partition(const std::vector<Cell>&p){
	int tw = 0;
	for(auto&x:p)
		if(tw < x.bag_size())
			tw = x.bag_size();
	return tw;
}

int get_node_count_of_multilevel_partition(const std::vector<Cell>&p){
	int node_count = 0;
	for(auto&x:p){
		for(auto&y:x.separator_node_list){
			if(node_count <= y)
				node_count = y+1;
		}
	}
	return node_count;
}
