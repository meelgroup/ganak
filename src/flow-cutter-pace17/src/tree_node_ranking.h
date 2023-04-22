#ifndef TREE_NODE_RANKING_H
#define TREE_NODE_RANKING_H
#include "array_id_func.h"
#include "count_range.h"
#include "id_sort.h"
#include <vector>

//! input graph must be a symmetric graph
//! The result is guaranteed to be optimal for trees. For non-trees there are no guarantees.
template<class Neighbors>
ArrayIDIDFunc compute_tree_node_ranking(const Neighbors&neighbors){
	const int node_count = neighbors.preimage_count();

	ArrayIDFunc<int>
		parent(node_count),
		child_first_order(node_count);

	// Root tree at node 0
	{
		ArrayIDFunc<int>stack(node_count);
		int stack_end = 0;
		int next_order_pos = node_count;
		
		stack[stack_end++] = 0;
		parent.fill(-1);		
		parent[0] = -2;
		

		while(stack_end != 0){
			auto x = stack[--stack_end];
			child_first_order[--next_order_pos] = x;
			for(auto y:neighbors(x)){
				if(parent(y) == -1){
					assert(y != 0);
					stack[stack_end++] = y;
					parent[y] = x;
				}
			}
		}
		parent[0] = -1;
		assert(next_order_pos == 0);
	}


	ArrayIDIDFunc level(node_count, 1);

	// Compute node levels
	{
		struct SubTreeInfo{
			std::vector<int>critical_list;
			int size;
		};
		ArrayIDFunc<std::vector<SubTreeInfo>>node_children_info(node_count);

		BitIDFunc crit(node_count);
		crit.fill(false);

		for(int i=0; i<node_count; ++i){
			auto x = child_first_order[i];

			auto&&children_info = node_children_info[x];

			int children_count = children_info.size();

			SubTreeInfo tree_info;

			if(children_count == 0){
				level[x] = 0;
				tree_info = {{1}, 1};
			}else if(children_count == 1){
				auto&child_info = children_info[0];

				int l = 1;
				while(!child_info.critical_list.empty() && child_info.critical_list.back() == l){
					child_info.critical_list.pop_back();
					++l;
				}
				tree_info = std::move(child_info);
				tree_info.critical_list.push_back(l);
				++tree_info.size;
				if(level.image_count() < l)
					level.set_image_count(l);
				level[x] = l-1;
			}else{

				tree_info.size = 1;
				for(int i=0; i<children_count; ++i)
					tree_info.size += children_info[i].size;

				// Move largest sub tree to the position 0
				{
					int largest_sub_tree = 0;
					int largest_sub_tree_size = children_info[0].size;
					for(int i=1; i<children_count; ++i)
						if(largest_sub_tree_size < children_info[i].size){
							largest_sub_tree = i;
							largest_sub_tree_size = children_info[i].size;
						}

					if(largest_sub_tree != 0)
						std::swap(children_info[0], children_info[largest_sub_tree]);
				}

				// Move second largest sub tree to the position 1
				{
					int second_largest_sub_tree = 1;
					int second_largest_sub_tree_size = children_info[1].size;
					for(int i=2; i<children_count; ++i)
						if(second_largest_sub_tree_size < children_info[i].size){
							second_largest_sub_tree = i;
							second_largest_sub_tree_size = children_info[i].size;
						}

					if(second_largest_sub_tree != 1)
						std::swap(children_info[1], children_info[second_largest_sub_tree]);
				}

				// variable names come from http://www.sciencedirect.com/science/article/pii/0020019089901610
				int max = 0;
				int p = 0;
				int q = 0;
				
				for(int i = 1; i<children_count; ++i)
					for(auto t:children_info[i].critical_list)
						if(!crit(t)){
							crit.set(t, true);
							if(t > max)
								max = t;
						}else{
							if(t > p)
								p = t;
						}

				auto&&first_child_critical_list = children_info[0].critical_list;

				for(int i = first_child_critical_list.size()-1; i>=0; --i){
					int t = first_child_critical_list[i];

					if(t > max)
						break;
					if(p >= t)
						continue;

					if(!crit(t)){
						crit.set(t, true);
					}else{
						p = t;
					}
				}

				for(int i=0; i<=p; ++i)
					crit.set(i, false);
				for(int i=p+1; i<=max; ++i){
					if(!crit(i)){
						if(q == 0){
							q = i;
							tree_info.critical_list = {q};
						}
					}else{
						crit.set(i, false);
						if(q != 0)
							tree_info.critical_list.push_back(i);
					}
				}

				if(q == 0){

					while(!first_child_critical_list.empty() && first_child_critical_list.back() <= max)
						first_child_critical_list.pop_back();
					q = max+1;
					while(!first_child_critical_list.empty() && first_child_critical_list.back() == q){
						first_child_critical_list.pop_back();
						++q;
					}

					tree_info.critical_list = std::move(first_child_critical_list);
					tree_info.critical_list.push_back(q);
				}else{
					while(!first_child_critical_list.empty() && first_child_critical_list.back() <= max)
						first_child_critical_list.pop_back();
					assert(std::is_sorted(tree_info.critical_list.begin(), tree_info.critical_list.end()));
					for(int i=tree_info.critical_list.size()-1; i>=0; --i)
						first_child_critical_list.push_back(tree_info.critical_list[i]);
					tree_info.critical_list = std::move(first_child_critical_list);
				}

				assert(q > 0);

				if(level.image_count() < q)
					level.set_image_count(q);
				level[x] = q-1;
			}

			assert(std::is_sorted(tree_info.critical_list.begin(), tree_info.critical_list.end(), std::greater<int>()));

			if(parent(x) != -1)
				node_children_info[parent(x)].push_back(std::move(tree_info));
		}
	}

	return level; // NVRO
}

#endif

