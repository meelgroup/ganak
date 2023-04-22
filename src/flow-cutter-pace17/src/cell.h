#ifndef CELL_H
#define CELL_H

#include <vector>

struct Cell{
	std::vector<int>separator_node_list;
	std::vector<int>boundary_node_list;

	int parent_cell;

	int bag_size()const{
		return separator_node_list.size() + boundary_node_list.size();
	}
};

inline
bool operator<(const Cell&l, const Cell&r){
	return l.bag_size() < r.bag_size();
}

int get_treewidth_of_multilevel_partition(const std::vector<Cell>&p);
int get_node_count_of_multilevel_partition(const std::vector<Cell>&p);

#endif

