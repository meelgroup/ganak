#ifndef CONTRACTION_GRAPH_H
#define CONTRACTION_GRAPH_H

#include "array_id_func.h"
#include "tiny_id_func.h"
#include "min_max.h"
#include "multi_arc.h"
#include <cassert>
#include <algorithm>

class EdgeContractionGraph{
public:
	void rewire_arcs_from_second_to_first(int u, int v){
		union_find_parent[v] = u;
		std::swap(next_adjacency_in_ring[u], next_adjacency_in_ring[v]);
	}

	template<class F>
	void forall_nodes_in_last_computed_neighborhood(const F&f){
		for(int i=0; i<neighborhood_size; ++i)
			f(neighborhood(i));
	}

	void compute_neighborhood_of(int v){
		for(int i=0; i<neighborhood_size; ++i)
			in_neighborhood.set(neighborhood(i), false);
		neighborhood_size = 0;

		if(union_find_parent[v] == v){
			const int initial_adjacency = v;
			int current_adjacency = v;
			do{
				// Iterate over the adjacency
				{
					int arc_in_begin = out_arc_begin[current_adjacency];
					int arc_in_end = out_arc_end[current_adjacency];
					
					int arc_out_begin = out_arc_begin[current_adjacency];
					
					while(arc_in_begin != arc_in_end){
						// Compress union find path
						{
							int x = arc_head[arc_in_begin];
							while(union_find_parent[x] != x)
								x = union_find_parent[x];
							int y = arc_head[arc_in_begin];
							while(union_find_parent[y] != y){
								int z = union_find_parent[y];
								union_find_parent[y] = x;
								y = z;
							}
							
						}

						// Replace arc head by representative from union find
						arc_head[arc_in_begin] = union_find_parent[arc_head[arc_in_begin]];

						assert(union_find_parent[arc_head[arc_in_begin]] == arc_head[arc_in_begin]);

						// Only keep the nodes that are not the heads of loops or multi arcs
						if(!in_neighborhood(arc_head[arc_in_begin]) && arc_head[arc_in_begin] != v){
							arc_head[arc_out_begin] = arc_head[arc_in_begin];
							++arc_out_begin;
							in_neighborhood.set(arc_head[arc_in_begin], true);
							neighborhood[neighborhood_size++] = arc_head[arc_in_begin];
						}

						++arc_in_begin;
					}

					out_arc_end[current_adjacency] = arc_out_begin;
				}
				
				// Goto next non-empty adjacency in the ring, and rewire the ring pointer to skip them in future
				int next_adjacency = next_adjacency_in_ring[current_adjacency];
				while(out_arc_begin[next_adjacency] == out_arc_end[next_adjacency] && next_adjacency != initial_adjacency)
					next_adjacency = next_adjacency_in_ring[next_adjacency];
				next_adjacency_in_ring[current_adjacency] = next_adjacency;
				current_adjacency = next_adjacency;
			}while(current_adjacency != initial_adjacency);
		}
	}

	template<class Tail, class Head>
	EdgeContractionGraph(const Tail&tail, const Head&head):
		next_adjacency_in_ring(tail.image_count()),
		union_find_parent(tail.image_count()),
		out_arc_begin(tail.image_count()),
		out_arc_end(tail.image_count()),
		arc_head(tail.preimage_count()),
		in_neighborhood(tail.image_count()),
		neighborhood(tail.image_count()),
		neighborhood_size(0)
	{
		assert(is_symmetric(tail, head));
		for(int i=0; i<tail.image_count(); ++i){
			next_adjacency_in_ring.set(i, i);
			union_find_parent.set(i, i);
		}

		in_neighborhood.fill(false);

		out_arc_end.fill(0);
		for(int i=0; i<tail.preimage_count(); ++i){
			int t = tail(i);
			out_arc_end.set(t, out_arc_end(t)+1);
		}
		if(tail.image_count() != 0){
			out_arc_begin.set(0, 0);
			for(int i=1; i<tail.image_count(); ++i){
				out_arc_begin.set(i, out_arc_end(i-1));
				out_arc_end.set(i, out_arc_end(i) + out_arc_begin(i));
			}
			assert(out_arc_end(tail.image_count()-1) == tail.preimage_count());
		}
		for(int i=0; i<tail.preimage_count(); ++i){
			int t = tail(i);
			arc_head.set(out_arc_begin(t), head(i));
			out_arc_begin.set(t, out_arc_begin(t)+1);
		}
		for(int i=0; i<tail.preimage_count(); ++i){
			int t = tail(i);
			out_arc_begin.set(t, out_arc_begin(t)-1);
		}
	}

private:
	ArrayIDFunc<int> next_adjacency_in_ring;
	ArrayIDFunc<int> union_find_parent;
	ArrayIDFunc<int> out_arc_begin;
	ArrayIDFunc<int> out_arc_end;
	ArrayIDFunc<int> arc_head;
	BitIDFunc in_neighborhood;
	ArrayIDFunc<int> neighborhood;
	int neighborhood_size;
};

class NodeContractionGraph{
public:
	template<class Tail, class Head>
	NodeContractionGraph(const Tail&tail, const Head&head):
		g(tail, head), is_virtual(tail.image_count()){
		assert(is_symmetric(tail, head));
		is_virtual.fill(false);
	}

	template<class F>
	void forall_neighbors_then_contract_node(int v, const F&callback){
		g.compute_neighborhood_of(v);
		g.forall_nodes_in_last_computed_neighborhood(
			[&](int u){
				if(is_virtual(u))
					g.rewire_arcs_from_second_to_first(v, u);
			}
		);
		is_virtual.set(v, true);
		g.compute_neighborhood_of(v);
		g.forall_nodes_in_last_computed_neighborhood(callback);
	}

private:
	EdgeContractionGraph g;
	BitIDFunc is_virtual;
};

template<class Tail, class Head, class OnNewArc>
int compute_chordal_supergraph(const Tail&tail, const Head&head, const OnNewArc&on_new_arc){
	assert(is_symmetric(tail, head));

	NodeContractionGraph g(tail, head);
	int max_upward_degree = 0;
	for(int x=0; x<tail.image_count()-1; ++x){
		int upward_degree = 0;
		g.forall_neighbors_then_contract_node(
			x, 
			[&](int y){
				on_new_arc(x, y); 
				++upward_degree;
			}
		);
		max_to(max_upward_degree, upward_degree);
	}
	return max_upward_degree;
}

/*
class SimpleNodeContractionAdjGraph{
public:
	template<class ArcList>
	SimpleNodeContractionAdjGraph(int node_count, const ArcList&arc_list):
		g(node_count, arc_list), is_virtual(node_count, false), in_neighborhood(node_count, false){
	}

	int node_count()const{
		return g.node_count();
	}

	int was_node_contracted(int v)const{
		assert(0 <= v && v < node_count() && "v is out of bounds");
		return g.was_node_contracted(v) || is_virtual[v];
	}

	void contract_node(int v){
		assert(0 <= v && v < node_count() && "v is out of bounds");
		assert(!was_node_contracted(v));

		is_virtual[v] = true;
		for(auto h:g.out(v))
			if(is_virtual[h.target])
				g.contract_arc(v, h.target);
	}

	std::vector<ArcHead> out_followed_by_contract_node(int v){
		contract_node(v);
		std::vector<ArcHead> r = g.out(v);
		if(r.size() == 1)
			g.contract_arc(r[0].target, v);
		return std::move(r);
	}

	std::vector<ArcHead>out(int v)const{
		assert(0 <= v && v < node_count() && "v is out of bounds");
		if(was_node_contracted(v))
			return {};
		else{
			#ifdef EXPENSIVE_CONTRACTION_GRAPH_TESTS
			for(auto x:in_neighborhood)
				assert(!x);
			#endif
			std::vector<ArcHead>neighborhood;		
			for(auto u:g.out(v)){
				if(!is_virtual[u.target]){
					if(!in_neighborhood[u.target]){
						in_neighborhood[u.target] = true;
						neighborhood.push_back(u);
					}
				}else{
					for(auto w:g.out(u.target)){
						assert(!is_virtual[w.target]);
						if(w.target != v){
							if(!in_neighborhood[w.target]){
								in_neighborhood[w.target] = true;
								neighborhood.push_back(w);
							}
						}
					}
				}
			}
			for(auto x:neighborhood)
				in_neighborhood[x.target] = false;
			#ifdef EXPENSIVE_CONTRACTION_GRAPH_TESTS
			for(auto x:in_neighborhood)
				assert(!x);
			#endif
			return std::move(neighborhood);
		}
	}

private:
	SimpleEdgeContractionAdjGraph g;
	std::vector<bool>is_virtual;
	mutable std::vector<bool>in_neighborhood;
};

template<class NodeList, class ArcList>
std::vector<Arc>compute_simple_contraction_hierarchy(const NodeList&node_list, const ArcList&arc_list){
	SimpleNodeContractionAdjGraph g(range_size(node_list), arc_list);
	std::vector<Arc>shortcuts;
	for(int v:node_list){
		for(auto u:g.out_followed_by_contract_node(v))
			shortcuts.push_back({v, u.target});
	}
	
	arc_sort(shortcuts.begin(), shortcuts.end());

	return std::move(shortcuts);
}*/

#endif
