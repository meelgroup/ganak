#include "tree_decomposition.hpp"
#include "io_helper.hpp"
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include "heap.hpp"
#include "tiny_id_func.hpp"
#include "contraction_graph.hpp"
#include "id_multi_func.hpp"
using namespace std;

void print_tree_decompostion_of_order(std::ostream&out, ArrayIDIDFunc tail, ArrayIDIDFunc head, const ArrayIDIDFunc&order){
	const int node_count = tail.image_count();

	auto inv_order = inverse_permutation(order);
	tail = chain(tail, inv_order);
	head = chain(head, inv_order);

	vector<vector<int>>nodes_in_bag;
	ArrayIDFunc<vector<int>>bags_of_node(node_count);

	auto is_left_subset_of_right = [](const std::vector<int>&l, const std::vector<int>&r){
		auto i = l.begin(), j = r.begin();

		for(;;){
			if(i == l.end())
				return true;
			if(j == r.end())
				return false;

			if(*i < *j)
				return false;
			if(*i == *j)
				++i;
			++j;
		}
	};

	auto compute_intersection_size = [](const std::vector<int>&l, const std::vector<int>&r){
		auto i = l.begin(), j = r.begin();
		int n = 0;
		for(;;){
			if(i == l.end() || j == r.end())
				return n;
			if(*i < *j)
				++i;
			else if(*i > *j)
				++j;
			else{
				++i;
				++j;
				++n;
			}
		}
	};

	auto on_new_potential_maximal_clique = [&](int lowest_node_in_clique, std::vector<int>clique){
		for(auto b:bags_of_node(lowest_node_in_clique))
			if(is_left_subset_of_right(clique, nodes_in_bag[b]))
				return;
		int bag_id = nodes_in_bag.size();
		for(auto x:clique)
			bags_of_node[x].push_back(bag_id);
		nodes_in_bag.push_back(std::move(clique));
	};


	{
		BitIDFunc is_root(node_count);
		is_root.fill(true);
		std::vector<int>upper_neighborhood_of_z;
		int z = -1;
		compute_chordal_supergraph(
			tail, head,
			[&](int x, int y){
				is_root.set(x, false);
				if(z != -1 && z != x){
					upper_neighborhood_of_z.push_back(z);
					sort(upper_neighborhood_of_z.begin(), upper_neighborhood_of_z.end());
					on_new_potential_maximal_clique(z, std::move(upper_neighborhood_of_z));
					upper_neighborhood_of_z.clear();
				}
				z = x;
				upper_neighborhood_of_z.push_back(y);
			}
		);
		if(z != -1){
			upper_neighborhood_of_z.push_back(z);
			sort(upper_neighborhood_of_z.begin(), upper_neighborhood_of_z.end());
			on_new_potential_maximal_clique(z, std::move(upper_neighborhood_of_z));
		}

		for(int x=0; x<node_count; ++x){
			if(is_root(x)){
				on_new_potential_maximal_clique(x, {x});
			}
		}
	}



	int bag_count = nodes_in_bag.size();

	int maximum_bag_size = 0;
	for(auto&b:nodes_in_bag)
		if((int)b.size() > maximum_bag_size)
			maximum_bag_size = b.size();

	out << "s td "<< nodes_in_bag.size() << ' ' << maximum_bag_size << ' ' << node_count << '\n';
	for(int i=0; i<bag_count; ++i){
		out << "b "<<(i+1);
		for(auto x:nodes_in_bag[i])
			out << ' ' << (order(x)+1);
		out << '\n';

	}

	{
		auto output_backbone_edge = [&](int b, int p){
			out << (b+1) << ' ' << (p+1) << '\n';
		};

		std::vector<int>tail2, head2, weight;

		for(int b=0; b<bag_count; ++b){
			vector<int>neighbor_bags;
			for(auto x:nodes_in_bag[b]){
				vector<int>tmp;
				std::set_union(
					bags_of_node[x].begin(), bags_of_node[x].end(),
					neighbor_bags.begin(), neighbor_bags.end(),
					std::back_inserter(tmp)
				);
				neighbor_bags.swap(tmp);
			}
			for(auto p:neighbor_bags){
				if(p != b){
					tail2.push_back(b);
					head2.push_back(p);
					weight.push_back(compute_intersection_size(nodes_in_bag[b], nodes_in_bag[p]));
				}
			}
		}

		int arc_count = tail2.size();

		auto out_arc = invert_id_id_func(
			id_id_func(
				arc_count, bag_count,
				[&](unsigned a){return tail2[a];}
			)
		);

		BitIDFunc in_tree(bag_count);
		in_tree.fill(false);
		max_id_heap<int>q(arc_count);

		for(int b=0; b<bag_count; ++b){
			if(!in_tree(b)){
				if(b != 0)
					output_backbone_edge(0, b);
				in_tree.set(b, true);
				for(int a:out_arc(b))
					q.push(a, weight[a]);
				while(!q.empty()){
					int xy = q.pop();

					int x = tail2[xy];
					int y = head2[xy];

					assert(in_tree(x));

					if(!in_tree(y)){
						output_backbone_edge(x, y);
						in_tree.set(y, true);
						for(int yz:out_arc(y)){
							assert(!q.contains(yz));
							q.push(yz, weight[yz]);
						}
					}
				}

			}
		}
	}
}

void print_tree_decompostion_of_multilevel_partition(std::ostream&out,
		const ArrayIDIDFunc&/*tail*/, const ArrayIDIDFunc&/*head*/,
		const ArrayIDIDFunc&to_input_node_id, const std::vector<Cell>&cell_list){
	int bag_count = cell_list.size();

	out << "s td "<< bag_count << ' ' << get_treewidth_of_multilevel_partition(cell_list) << ' ' << get_node_count_of_multilevel_partition(cell_list) << '\n';

	for(int i=0; i<bag_count; ++i){
		out << "b "<<(i+1);
		for(auto&x:cell_list[i].separator_node_list)
			out << ' ' << (to_input_node_id(x)+1);
		for(auto&x:cell_list[i].boundary_node_list)
			out << ' ' << (to_input_node_id(x)+1);
		out << '\n';
	}

	for(int i=0; i<bag_count; ++i){
		if(cell_list[i].parent_cell != -1)
			out << (i+1) << ' ' << (cell_list[i].parent_cell+1) << '\n';
	}
}

