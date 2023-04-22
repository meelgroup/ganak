#ifndef CONNECTED_COMPONENTS_H
#define CONNECTED_COMPONENTS_H

#include "tiny_id_func.h"
#include "array_id_func.h"
#include "union_find.h"
#include "filter.h"
#include "chain.h"
#include "back_arc.h"

template<class Tail, class Head>
ArrayIDIDFunc compute_connected_components(const Tail&tail, const Head&head){
	const int node_count = tail.image_count();
	const int arc_count = tail.preimage_count();

	UnionFind uf(node_count);
	for(int i=0; i<arc_count; ++i)
		uf.unite(tail(i), head(i));
	
	auto is_representative = id_func(
		node_count, 
		[&](int i){
			return uf.is_representative(i);
		}
	);
	
	return chain(uf, compute_keep_function(is_representative, count_true(is_representative)));
}

template<class Tail, class Head>
bool is_connected(const Tail&tail, const Head&head){
	const int node_count = tail.image_count();
	const int arc_count = tail.preimage_count();

	if(node_count == 0){
		return true;
	} else {
		UnionFind uf(node_count);
		for(int i=0; i<arc_count; ++i)
			uf.unite(tail(i), head(i));
	
		return  uf.component_size(uf(0)) == node_count;
	}
}

template<
	class OutArc, class Head, 
	class OnRootFirstVisit,
	class OnRootLastVisit,
	class OnTreeUpArcVisit,
	class OnTreeDownArcVisit,
	class OnNonTreeArcVisit
>
void symmetric_depth_first_search(
	const OutArc&out_arc, 
	const Head&head,
	const OnRootFirstVisit&on_root_first_visit, 
	const OnRootLastVisit&on_root_last_visit,
	const OnTreeUpArcVisit&on_tree_down_arc_visit,
	const OnTreeDownArcVisit&on_tree_up_arc_visit,
	const OnNonTreeArcVisit&on_non_tree_arc_visit
){
	const int arc_count = head.preimage_count();
	const int node_count = out_arc.preimage_count();

	(void)arc_count;
	(void)node_count;

	ArrayIDFunc<int> dfs_stack(node_count);  
	int dfs_stack_end = 0;

	ArrayIDFunc<int> parent_arc(node_count);
	parent_arc.fill(-1);

	ArrayIDFunc<int> parent_node(node_count);
	parent_node.fill(-1);
		
	typedef typename std::decay<decltype(out_arc(0).begin())>::type Iter;
	ArrayIDFunc<Iter>next_out(node_count);
	for(int i=0; i<node_count; ++i)
		next_out[i] = std::begin(out_arc(i));

	for(int r=0; r<node_count; ++r){
		if(parent_node[r] == -1){
			dfs_stack_end = 0;
			parent_arc[r] = -2;
			parent_node[r] = -2;

			on_root_first_visit(r);

			int x = r; // the current node
			for(;;){
				if(next_out[x] == std::end(out_arc(x))){
					if(parent_arc[x] == -2)
						break;
					assert(0 <= parent_arc[x] && parent_arc[x] < arc_count);
					assert(0 <= parent_node[x] && parent_node[x] < node_count);
					on_tree_up_arc_visit(x, parent_arc[x], parent_node[x]);
					assert(dfs_stack_end != 0);
					x = dfs_stack[--dfs_stack_end];
				}else{
					int xy = *next_out[x]++;
					int y = head(xy);
					if(y == parent_node[x]){
						parent_arc[x] = xy;
					}else{
						if(parent_node[y] == -1){
							dfs_stack[dfs_stack_end++] = x;
							parent_node[y] = x;
							on_tree_down_arc_visit(x, xy, y);
							x = y;
						}else{
							on_non_tree_arc_visit(x, xy, y);
						}
					}
				}
				
			}

			on_root_last_visit(r);
		}
	}

}

template<
	class OutArc, class Head, class BackArc
>
ArrayIDIDFunc compute_biconnected_components(
	const OutArc&out_arc, const Head&head, const BackArc&back_arc
){
	const int node_count = out_arc.preimage_count();
	const int arc_count = head.preimage_count();

	(void)arc_count;
	(void)node_count;

	ArrayIDFunc<int> arc_stack(arc_count);
	int arc_stack_end = 0;

	ArrayIDIDFunc arc_component(arc_count, 0);
	arc_component.fill(-1);

	ArrayIDFunc<int> depth(node_count);
	ArrayIDFunc<int> min_succ_depth(node_count);

	auto min_to = [](int&x, int y){
		if(y < x)
			x = y;
	};

	auto on_first_root_visit = [&](int x){
		depth[x] = 0;
	};

	auto on_last_root_visit = [&](int x){

	};

	auto on_tree_down_arc_visit = [&](int x, int xy, int y){
		arc_stack[arc_stack_end++] = xy;

		min_succ_depth[y] = std::numeric_limits<int>::max();
		depth[y] = depth[x]+1;
	};

	auto on_tree_up_arc_visit = [&](int x, int xy, int y){
		arc_stack[arc_stack_end++] = xy;

		min_to(min_succ_depth[y], min_succ_depth[x]);
		min_to(min_succ_depth[y], depth[x]);

		if(min_succ_depth[x] >= depth[y]){
			const int new_component_id = arc_component.image_count();
			arc_component.set_image_count(arc_component.image_count() + 1);

			while(arc_stack_end != 0){
				int ab = arc_stack[--arc_stack_end];
				int ba = back_arc(ab);
				if(arc_component[ba] == -1){
					assert(arc_component[ab] == -1);
					arc_component[ab] = new_component_id;
					arc_component[ba] = new_component_id;
				}
				if(ba == xy)
					break;
			}
		}
	};

	auto on_non_tree_arc_visit = [&](int x, int xy, int y){
		arc_stack[arc_stack_end++] = xy;

		min_to(min_succ_depth[x], depth[y]);
	};

	symmetric_depth_first_search(
		out_arc, head,
		on_first_root_visit, on_last_root_visit,
		on_tree_down_arc_visit, on_tree_up_arc_visit,
		on_non_tree_arc_visit
	);

	#ifndef NDEBUG
	for(int i=0; i<arc_count; ++i)
		assert(arc_component[i] != -1);
	#endif

	return arc_component; // NVRO
}

template<class Tail, class Head>
bool is_biconnected(const Tail&tail, const Head&head){
	return compute_biconnected_components(invert_id_id_func(tail), head, compute_back_arc_permutation(tail, head)).image_count() <= 1;
}

#endif

