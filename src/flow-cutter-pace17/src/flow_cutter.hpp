#ifndef FLOW_CUTTER_H
#define FLOW_CUTTER_H

#include "tiny_id_func.hpp"
#include "array_id_func.hpp"
#include "id_func.hpp"
#include "min_max.hpp"
#include <vector>
#include <algorithm>
#include <sstream>
#include <random>
#include <memory>

#include "flow_cutter_config.hpp"

namespace flow_cutter{

	template<class Tail, class Head, class BackArc, class Capacity, class OutArc>
	struct Graph{
		Graph(
			Tail _tail,
			Head _head,
			BackArc _back_arc,
			Capacity _capacity,
			OutArc _out_arc
		):
			tail(std::move(_tail)),
			head(std::move(_head)),
			back_arc(std::move(_back_arc)),
			capacity(std::move(_capacity)),
			out_arc(std::move(_out_arc)){}

		Tail tail;
		Head head;
		BackArc back_arc;
		Capacity capacity;
		OutArc out_arc;

		int node_count()const{
			return tail.image_count();
		}

		int arc_count()const{
			return tail.preimage_count();
		}
	};

	struct TemporaryData{
		TemporaryData(){}
		explicit TemporaryData(int node_count):
			node_space(node_count){}
		ArrayIDFunc<int>node_space;
	};

	template<class Tail, class Head, class BackArc, class Capacity, class OutArc>
	Graph<Tail, Head, BackArc, Capacity, OutArc>
		make_graph(
			Tail tail, Head head, BackArc back_arc,
			Capacity capacity, OutArc out_arc
		){
		return {std::move(tail), std::move(head), std::move(back_arc), std::move(capacity), std::move(out_arc)};
	}

	template<class Tail, class Head, class BackArc, class OutArc>
	Graph<
		Tail, Head, BackArc,
		ConstIntIDFunc<1>,
		OutArc
	>
		make_graph(
			Tail tail, Head head,
			BackArc back_arc, OutArc out_arc
		){
		return {
			std::move(tail), std::move(head), std::move(back_arc),
			ConstIntIDFunc<1>(tail.preimage_count()),
			std::move(out_arc)
		};
	}

	class PseudoDepthFirstSearch{
	public:
		template<class Graph, class WasNodeSeen, class SeeNode, class ShouldFollowArc, class OnNewArc>
		void operator()(
			const Graph&graph, TemporaryData&tmp, int source_node,
			const WasNodeSeen&was_node_seen, const SeeNode&see_node,
			const ShouldFollowArc&should_follow_arc, const OnNewArc&on_new_arc
		)const{
			int stack_end = 1;
			auto&stack = tmp.node_space;
			stack[0] = source_node;

			while(stack_end != 0){
				int x = stack[--stack_end];
				for(auto xy : graph.out_arc(x)){
					on_new_arc(xy);
					int y = graph.head(xy);
					if(!was_node_seen(y)){
						if(should_follow_arc(xy)){
							if(!see_node(y))
								return;
							stack[stack_end++] = y;
						}
					}
				}
			}
		}
	};

	class BreadthFirstSearch{
	public:
		template<class Graph, class WasNodeSeen, class SeeNode, class ShouldFollowArc, class OnNewArc>
		void operator()(
			const Graph&graph, TemporaryData&tmp, int source_node,
			const WasNodeSeen&was_node_seen, const SeeNode&see_node,
			const ShouldFollowArc&should_follow_arc, const OnNewArc&on_new_arc
		)const{
			int queue_begin = 0, queue_end = 1;
			auto&queue = tmp.node_space;
			queue[0] = source_node;
			while(queue_begin != queue_end){
				int x = queue[queue_begin++];
				for(auto xy : graph.out_arc(x)){
					on_new_arc(xy);
					int y = graph.head(xy);
					if(!was_node_seen(y)){
						if(should_follow_arc(xy)){
							if(!see_node(y))
								return;
							queue[queue_end++] = y;
						}
					}
				}
			}
		}
	};

	struct UnitFlow{
		UnitFlow(){}
		explicit UnitFlow(int preimage_count):flow(preimage_count){}

		void clear(){
			flow.fill(1);
		}

		int preimage_count()const{
			return flow.preimage_count();
		}

		template<class Graph>
		void increase(const Graph&graph, int a){
			auto f = flow(a);
#ifdef SLOW_DEBUG
			assert((f == 0 || f == 1) && "Flow is already maximum; can not be increased");
			assert(flow(graph.back_arc(a)) == 2-f && "Back arc has invalid flow");
#endif
			++f;
			flow.set(a, f);
			flow.set(graph.back_arc(a), 2-f);
		}

		template<class Graph>
		void decrease(const Graph&graph, int a){
			auto f = flow(a);
#ifdef SLOW_DEBUG
			assert((f == 1 || f == 2) && "Flow is already minimum; can not be decreased");
			assert(flow(graph.back_arc(a)) == 2-f && "Back arc has invalid flow");
#endif
			--f;
			flow.set(a, f);
			flow.set(graph.back_arc(a), 2-f);
		}

		int operator()(int a)const{
			return static_cast<int>(flow(a))-1;
		}

		void swap(UnitFlow&o){
			flow.swap(o.flow);
		}

		TinyIntIDFunc<2>flow;
	};

	class BasicNodeSet{
	public:
		template<class Graph>
		explicit BasicNodeSet(const Graph&graph):
			node_count_inside_(0),
			inside_flag(graph.node_count()),
			extra_node(-1){}

		void clear(){
			node_count_inside_ = 0;
			inside_flag.fill(false);
		}

		bool can_grow()const{
			return extra_node != -1;
		}

		template<class Graph, class SearchAlgorithm, class OnNewNode, class ShouldFollowArc, class OnNewArc>
		void grow(
			const Graph&graph,
			TemporaryData&tmp,
			const SearchAlgorithm&search_algo,
			const OnNewNode&on_new_node, // on_new_node(x) is called for every node x. If it returns false then the search is stopped, if it returns true it continues
			const ShouldFollowArc&should_follow_arc, // is called for a subset of arcs and must say whether the arc sould be followed
			const OnNewArc&on_new_arc // on_new_arc(xy) is called for ever arc xy with x in the set
		){
			assert(can_grow());

			auto see_node = [&](int x){
#ifdef SLOW_DEBUG
				assert(!inside_flag(x));
#endif
				inside_flag.set(x, true);
				++this->node_count_inside_;
				return on_new_node(x);
			};

			auto was_node_seen = [&](int x){
				return inside_flag(x);
			};

			search_algo(graph, tmp, extra_node, was_node_seen, see_node, should_follow_arc, on_new_arc);
			extra_node = -1;
		}

		template<class Graph>
		void set_extra_node(const Graph& /*graph*/, int x){
#ifdef SLOW_DEBUG
			assert(!inside_flag(x));
			assert(extra_node == -1);
#endif
			inside_flag.set(x, true);
			++node_count_inside_;
			extra_node = x;
		}

		bool is_inside(int x) const {
			return inside_flag(x);
		}

		int node_count_inside() const {
			return node_count_inside_;
		}

		int max_node_count_inside() const {
			return inside_flag.preimage_count();
		}

	private:
		int node_count_inside_;
		BitIDFunc inside_flag;
		int extra_node;
	};

	class ReachableNodeSet;

	class AssimilatedNodeSet{
		friend class ReachableNodeSet;
	public:
		template<class Graph>
		explicit AssimilatedNodeSet(const Graph&graph):
			node_set(graph){}

		void clear(){
			node_set.clear();
			front.clear();
		}

		template<class Graph>
		void set_extra_node(const Graph&graph, int x){
			node_set.set_extra_node(graph, x);
		}

		bool can_grow()const{
			return node_set.can_grow();
		}

		template<class Graph, class SearchAlgorithm, class OnNewNode, class ShouldFollowArc, class OnNewArc, class HasFlow>
		void grow(
			const Graph&graph,
			TemporaryData&tmp,
			const SearchAlgorithm&search_algo,
			const OnNewNode&on_new_node, // on_new_node(x) is called for every node x. If it returns false then the search is stopped, if it returns true it continues
			const ShouldFollowArc&should_follow_arc, // is called for a subset of arcs and must say whether the arc sould be followed
			const OnNewArc&on_new_arc, // on_new_arc(xy) is called for ever arc xy with x in the set
			const HasFlow&has_flow
		){
			auto my_on_new_arc = [&](int xy){
				if(has_flow(xy))
					front.push_back(xy);
				on_new_arc(xy);
			};

			node_set.grow(graph, tmp, search_algo, on_new_node, should_follow_arc, my_on_new_arc);
		}

		bool is_inside(int x) const {
			return node_set.is_inside(x);
		}

		int node_count_inside() const {
			return node_set.node_count_inside();
		}

		int max_node_count_inside() const {
			return node_set.max_node_count_inside();
		}

		template<class Graph>
		void shrink_cut_front(const Graph&graph){
			front.erase(
				std::remove_if(
					front.begin(), front.end(),
					[&](int xy){ return node_set.is_inside(graph.head(xy)); }
				),
				front.end()
			);
		}

		const std::vector<int>&get_cut_front() const {
			return front;
		}

	private:
		BasicNodeSet node_set;
		std::vector<int>front;
	};

	class ReachableNodeSet{
	public:
		template<class Graph>
		explicit ReachableNodeSet(const Graph&graph):
			node_set(graph), predecessor(graph.node_count()){}

		void reset(const AssimilatedNodeSet&other){
			node_set = other.node_set;
		}

		void clear(){
			node_set.clear();
		}

		template<class Graph>
		void set_extra_node(const Graph&graph, int x){
			node_set.set_extra_node(graph, x);
		}

		bool can_grow()const{
			return node_set.can_grow();
		}

		template<class Graph, class SearchAlgorithm, class OnNewNode, class ShouldFollowArc, class OnNewArc>
		void grow(
			const Graph&graph,
			TemporaryData&tmp,
			const SearchAlgorithm&search_algo,
			const OnNewNode&on_new_node, // on_new_node(x) is called for every node x. If it returns false then the search is stopped, if it returns true it continues
			const ShouldFollowArc&should_follow_arc, // is called for a subset of arcs and must say whether the arc sould be followed
			const OnNewArc&on_new_arc // on_new_arc(xy) is called for ever arc xy with x in the set
		){
			auto my_should_follow_arc = [&](int xy){
				predecessor[graph.head(xy)] = xy;
				return should_follow_arc(xy);
			};

			node_set.grow(graph, tmp, search_algo, on_new_node, my_should_follow_arc, on_new_arc);
		}

		bool is_inside(int x) const {
			return node_set.is_inside(x);
		}

		int node_count_inside() const {
			return node_set.node_count_inside();
		}

		int max_node_count_inside() const {
			return node_set.max_node_count_inside();
		}

		template<class Graph, class IsSource, class OnNewArc>
		void forall_arcs_in_path_to(const Graph&graph, const IsSource&is_source, int target, const OnNewArc&on_new_arc){
			int x = target;
			while(!is_source(x)){
				on_new_arc(predecessor[x]);
				x = graph.tail(predecessor[x]);
			}
		}

	private:
		BasicNodeSet node_set;
		ArrayIDFunc<int>predecessor;
	};

	struct SourceTargetPair{
		int source, target;
	};

	struct CutterStateDump{
		BitIDFunc source_assimilated, target_assimilated, source_reachable, target_reachable, flow;
	};

	class BasicCutter{
	public:
		template<class Graph>
		explicit BasicCutter(const Graph&graph):
			assimilated{AssimilatedNodeSet(graph), AssimilatedNodeSet(graph)},
			reachable{ReachableNodeSet(graph), ReachableNodeSet(graph)},
			flow(graph.arc_count()),
			cut_available(false)
		{}

		template<class Graph, class SearchAlgorithm>
		void init(const Graph&graph, TemporaryData&tmp, const SearchAlgorithm&search_algo, SourceTargetPair p){
			assimilated[source_side].clear();
			reachable[source_side].clear();
			assimilated[target_side].clear();
			reachable[target_side].clear();
			flow.clear();

			assimilated[source_side].set_extra_node(graph, p.source);
			reachable[source_side].set_extra_node(graph, p.source);
			assimilated[target_side].set_extra_node(graph, p.target);
			reachable[target_side].set_extra_node(graph, p.target);

			grow_reachable_sets(graph, tmp, search_algo, source_side);
			grow_assimilated_sets(graph, tmp, search_algo);

			cut_available = true;
			check_invariants(graph);
		}

		CutterStateDump dump_state()const{
			return {
				id_func(
					assimilated[source_side].max_node_count_inside(),
					[&](int x){
						return assimilated[source_side].is_inside(x);
					}
				),
				id_func(
					assimilated[target_side].max_node_count_inside(),
					[&](int x){
						return assimilated[target_side].is_inside(x);
					}
				),
				id_func(
					assimilated[source_side].max_node_count_inside(),
					[&](int x){
						return reachable[source_side].is_inside(x);
					}
				),
				id_func(
					assimilated[target_side].max_node_count_inside(),
					[&](int x){
						return reachable[target_side].is_inside(x);
					}
				),
				id_func(
					flow.preimage_count(),
					[&](int xy){
						return flow(xy) != 0;
					}
				)
			};
		}

		//! Returns true if a new cut was found. Returns false if no cut was found. False implies that no cut
		//! will be found in the future. Repeatly calling this function after it returned false does not do
		//! anything.
		template<class Graph, class SearchAlgorithm, class ScorePierceNode>
		bool advance(const Graph&graph, TemporaryData&tmp, const SearchAlgorithm&search_algo, const ScorePierceNode&score_pierce_node){
			assert(cut_available);

			check_invariants(graph);
			int side = get_current_cut_side();
			if(assimilated[side].node_count_inside() >= graph.node_count()/2){
				cut_available = false;
				return false;
			}

			int pierce_node = select_pierce_node(graph, side, score_pierce_node);

			if(pierce_node == -1){
				cut_available = false;
				return false;
			}

#ifdef SLOW_DEBUG
			assert(!assimilated[1-side].is_inside(pierce_node));
#endif

			assimilated[side].set_extra_node(graph, pierce_node);
			reachable[side].set_extra_node(graph, pierce_node);

			grow_reachable_sets(graph, tmp, search_algo, side);
			grow_assimilated_sets(graph, tmp, search_algo);
			check_invariants(graph);
			cut_available = true;
			return true;
		}

		bool is_cut_available()const{
			return cut_available;
		}

		template<class Graph, class ScorePierceNode>
		bool does_next_advance_increase_cut(const Graph&graph, const ScorePierceNode&score_pierce_node){
			int side = get_current_cut_side();

			if(assimilated[side].node_count_inside() >= graph.node_count()/2){
				return true;
			}


			int pierce_node = select_pierce_node(graph, side, score_pierce_node);

			if(pierce_node == -1)
				return true;
			else if(reachable[1-side].is_inside(pierce_node))
				return true;
			else
				return false;
		}

		bool is_on_smaller_side(int x)const{
			return assimilated[get_current_cut_side()].is_inside(x);
		}

		static const int source_side = 0;
		static const int target_side = 1;

		int get_current_cut_side()const{
			if(
				reachable[source_side].node_count_inside() == assimilated[source_side].node_count_inside() && (
					reachable[target_side].node_count_inside() != assimilated[target_side].node_count_inside() ||
					assimilated[source_side].node_count_inside() <= assimilated[target_side].node_count_inside()
				)
			)
				return source_side;
			else
				return target_side;
		}

		int get_current_smaller_cut_side_size()const{
			return assimilated[get_current_cut_side()].node_count_inside();
		}

		const std::vector<int>&get_current_cut()const{
			return assimilated[get_current_cut_side()].get_cut_front();
		}

		int get_assimilated_node_count()const{
			return assimilated[source_side].node_count_inside() + assimilated[target_side].node_count_inside();
		}

	private:
		template<class Graph, class ScorePierceNode>
		int select_pierce_node(const Graph&graph, int side, const ScorePierceNode&score_pierce_node){

			int pierce_node = -1;
			int max_score = std::numeric_limits<int>::min();
			for(auto xy : assimilated[side].get_cut_front()){
				int y = graph.head(xy);
				if(!assimilated[1-side].is_inside(y)){
					int score = score_pierce_node(y, side, reachable[1-side].is_inside(y));

					if(score > max_score){
						max_score = score;
						pierce_node = y;
					}
				}
			}

			return pierce_node;
		}

		template<class Graph>
		bool is_saturated(const Graph&graph, int direction, int xy){
			if(direction == target_side)
				xy = graph.back_arc(xy);
			return graph.capacity(xy) == flow(xy);
		}


		template<class Graph, class SearchAlgorithm>
		void grow_reachable_sets(const Graph&graph, TemporaryData&tmp, const SearchAlgorithm&search_algo, int pierced_side){

			int my_source_side = pierced_side;
			int my_target_side = 1-pierced_side;

#ifdef SLOW_DEBUG
			assert(reachable[pierced_side].can_grow());
#endif

			auto is_forward_saturated = [&,this](int xy){
				return this->is_saturated(graph, my_source_side, xy);
			};

			auto is_backward_saturated = [&,this](int xy){
				return this->is_saturated(graph, my_target_side, xy);
			};

			auto is_source = [&](int x){
				return assimilated[my_source_side].is_inside(x);
			};

			auto is_target = [&](int x){
				return assimilated[my_target_side].is_inside(x);
			};

			auto increase_flow = [&](int xy){
				if(pierced_side == source_side)
					flow.increase(graph, xy);
				else
					flow.decrease(graph, xy);
			};

			bool was_flow_augmented = false;

			int target_hit;
			do{
				target_hit = -1;
				auto on_new_node = [&](int x){
					if(is_target(x)){
						target_hit = x;
						return false;
					} else
						return true;
				};
				auto should_follow_arc = [&](int xy){ return !is_forward_saturated(xy); };
				auto on_new_arc = [](int /*xy*/){};
				reachable[my_source_side].grow(graph, tmp, search_algo, on_new_node, should_follow_arc, on_new_arc);

				if(target_hit != -1){
					check_flow_conservation(graph);
					reachable[my_source_side].forall_arcs_in_path_to(graph, is_source, target_hit, increase_flow);
					check_flow_conservation(graph);
					reachable[my_source_side].reset(assimilated[my_source_side]);

					was_flow_augmented = true;
					check_flow_conservation(graph);
				}
			}while(target_hit != -1);

			if(was_flow_augmented){
				reachable[my_target_side].reset(assimilated[my_target_side]);
				auto on_new_node = [&](int /*x*/){return true;};
				auto should_follow_arc = [&](int xy){ return !is_backward_saturated(xy); };
				auto on_new_arc = [](int/* xy*/){};
				reachable[my_target_side].grow(graph, tmp, search_algo, on_new_node, should_follow_arc, on_new_arc);
			}

		}

		template<class Graph, class SearchAlgorithm>
		void grow_assimilated_sets(const Graph&graph, TemporaryData&tmp, const SearchAlgorithm&search_algo){
			auto is_forward_saturated = [&,this](int xy){
				return this->is_saturated(graph, source_side, xy);
			};

			auto is_backward_saturated = [&,this](int xy){
				return this->is_saturated(graph, target_side, xy);
			};

			if(reachable[source_side].node_count_inside() <= reachable[target_side].node_count_inside()){
				auto on_new_node = [&](int/* x*/){return true;};
				auto should_follow_arc = [&](int xy){ return !is_forward_saturated(xy); };
				auto on_new_arc = [](int/* xy*/){};
				auto has_flow = [&](int xy){ return flow(xy) != 0; };
				assimilated[source_side].grow(graph, tmp, search_algo, on_new_node, should_follow_arc, on_new_arc, has_flow);
				assimilated[source_side].shrink_cut_front(graph);
			}else{
				auto on_new_node = [&](int/* x*/){return true;};
				auto should_follow_arc = [&](int xy){ return !is_backward_saturated(xy); };
				auto on_new_arc = [](int/* xy*/){};
				auto has_flow = [&](int xy){ return flow(xy) != 0; };
				assimilated[target_side].grow(graph, tmp, search_algo, on_new_node, should_follow_arc, on_new_arc, has_flow);
				assimilated[target_side].shrink_cut_front(graph);
			}
		}

		template<class Graph>
		void check_flow_conservation([[maybe_unused]] const Graph& graph){
			#ifdef SLOW_DEBUG
			for(int x=0; x<graph.node_count(); ++x)
				if(!assimilated[source_side].is_inside(x) && !assimilated[target_side].is_inside(x)){
					int flow_surplus = 0;
					for(auto xy : graph.out_arc(x))
						flow_surplus += flow(xy);
					assert(flow_surplus == 0 && "Flow must be conserved outside of the assimilated sides");
				}
			#endif
		}

		template<class Graph>
		void check_invariants([[maybe_unused]] const Graph&graph){
			#ifdef SLOW_DEBUG
			for(int side = 0; side < 2; ++side)
				assert(assimilated[side].node_count_inside() > 0 && "Each side must contain at least one node");

			for(int x=0; x<graph.node_count(); ++x)
				assert((!assimilated[source_side].is_inside(x) || !assimilated[target_side].is_inside(x)) && "a node can not be assimilated by both sides");

			for(int side = 0; side < 2; ++side)
				for(int x=0; x<graph.node_count(); ++x)
					if(assimilated[side].is_inside(x))
						assert(reachable[side].is_inside(x) && "assimilated must be a subset of reachable");

			check_flow_conservation(graph);

			int smaller_reachable_side;
			if(reachable[source_side].node_count_inside() <= reachable[target_side].node_count_inside())
				smaller_reachable_side = source_side;
			else
				smaller_reachable_side = target_side;
			assert(reachable[smaller_reachable_side].node_count_inside() == assimilated[smaller_reachable_side].node_count_inside());
			for(int x=0; x<graph.node_count(); ++x)
				assert(reachable[smaller_reachable_side].is_inside(x) == assimilated[smaller_reachable_side].is_inside(x));

			assert(!reachable[source_side].can_grow());
			assert(!reachable[target_side].can_grow());
			assert(!assimilated[smaller_reachable_side].can_grow());
			#endif
		}

		AssimilatedNodeSet assimilated[2];
		ReachableNodeSet reachable[2];
		UnitFlow flow;
		bool cut_available;
	};



	enum class DistanceType{
		no_distance,
		hop_distance,
		weighted_distance
	};

	class DistanceAwareCutter{
	private:
		template<class Graph>
		static void compute_hop_distance_from(const Graph&graph, TemporaryData&tmp, int source, ArrayIDFunc<int>&dist){
			dist.fill(std::numeric_limits<int>::max());
			dist[source] = 0;

			auto was_node_seen = [&](int/* x*/){return false;};
			auto see_node = [](int/* x*/){ return true; };
			auto should_follow_arc = [&](int xy){
				if(dist(graph.tail(xy)) < dist(graph.head(xy)) - 1){
					dist[graph.head(xy)] = dist(graph.tail(xy))+1;
					return true;
				}else{
					return false;
				}
			};
			auto on_new_arc = [&](int/* xy*/){};
			BreadthFirstSearch()(graph, tmp, source, was_node_seen, see_node, should_follow_arc, on_new_arc);
		}
	public:
		template<class Graph>
		DistanceAwareCutter(const Graph&graph):
			cutter(graph),
			node_dist{ArrayIDFunc<int>{graph.node_count()}, ArrayIDFunc<int>{graph.node_count()}}{}

		template<class Graph, class SearchAlgorithm>
		void init(const Graph&graph, TemporaryData&tmp, const SearchAlgorithm&search_algo, DistanceType dist_type, SourceTargetPair p, int random_seed){
			cutter.init(graph, tmp, search_algo, p);

			rng.seed(random_seed);

			switch(dist_type){
			case DistanceType::hop_distance:
				compute_hop_distance_from(graph, tmp, p.source, node_dist[source_side]);
				compute_hop_distance_from(graph, tmp, p.target, node_dist[target_side]);
				break;
			case DistanceType::no_distance:
				break;
			default:
				assert(false);
				break;
			}
		}

		CutterStateDump dump_state()const{
			return cutter.dump_state();
		}

		template<class Graph, class SearchAlgorithm, class ScorePierceNode>
		bool advance(const Graph&graph, TemporaryData&tmp, const SearchAlgorithm&search_algo, const ScorePierceNode&score_pierce_node){
			auto my_score_pierce_node = [&](int x, int side, bool causes_augmenting_path){
				return score_pierce_node(x, side, causes_augmenting_path, node_dist[side](x), node_dist[1-side](x));
			};
			return cutter.advance(graph, tmp, search_algo, my_score_pierce_node);
		}

		bool is_cut_available()const{
			return cutter.is_cut_available();
		}

		template<class Graph, class ScorePierceNode>
		bool does_next_advance_increase_cut(const Graph&graph, const ScorePierceNode&score_pierce_node){
			auto my_score_pierce_node = [&](int x, int side, bool causes_augmenting_path){
				return score_pierce_node(x, side, causes_augmenting_path, node_dist[side](x), node_dist[1-side](x));
			};
			return cutter.does_next_advance_increase_cut(graph, my_score_pierce_node);
		}

		static const int source_side = BasicCutter::source_side;
		static const int target_side = BasicCutter::target_side;

		int get_current_cut_side()const{
			return cutter.get_current_cut_side();
		}

		int get_current_smaller_cut_side_size()const{
			return cutter.get_current_smaller_cut_side_size();
		}

		const std::vector<int>&get_current_cut()const{
			return cutter.get_current_cut();
		}

		int get_assimilated_node_count()const{
			return cutter.get_assimilated_node_count();
		}

		bool is_on_smaller_side(int x)const{
			return cutter.is_on_smaller_side(x);
		}

		bool is_empty()const{
			return node_dist[0].preimage_count() == 0;
		}
	private:
		BasicCutter cutter;
		ArrayIDFunc<int>node_dist[2];
		std::mt19937 rng;
	};

	class MultiCutter{
	public:
		MultiCutter(){}

		template<class Graph, class SearchAlgorithm,  class ScorePierceNode>
		void init(
			const Graph&graph, TemporaryData&tmp,
			const SearchAlgorithm&search_algo, const ScorePierceNode&score_pierce_node, DistanceType dist_type,
			const std::vector<SourceTargetPair>&p, int random_seed, bool should_skip_non_maximum_sides = true
		){
			while(cutter_list.size() > p.size())
				cutter_list.pop_back(); // can not use resize because that requires default constructor...
			while(cutter_list.size() < p.size())
				cutter_list.emplace_back(graph);

			for(int i=0; i<(int)p.size(); ++i){
				auto&x = cutter_list[i];
				auto my_score_pierce_node = [&](int x2, int side, bool causes_augmenting_path, int source_dist, int target_dist){
					return score_pierce_node(x2, side, causes_augmenting_path, source_dist, target_dist, i);
				};

				x.init(graph, tmp, search_algo, dist_type, p[i], random_seed+1+i);
				if(should_skip_non_maximum_sides)
					while(!x.does_next_advance_increase_cut(graph, my_score_pierce_node))
						x.advance(graph, tmp, search_algo, my_score_pierce_node);
			}

			int best_cutter_id = -1;
			int best_cut_size = std::numeric_limits<int>::max();
			int best_cutter_weight = 0;

			for(int i=0; i<(int)p.size(); ++i){
				auto&x = cutter_list[i];
				if(
					(int)x.get_current_cut().size() < best_cut_size
					|| (
						(int)x.get_current_cut().size() == best_cut_size &&
						x.get_current_smaller_cut_side_size() > best_cutter_weight
					)
				){
					best_cutter_id = i;
					best_cut_size = x.get_current_cut().size();
					best_cutter_weight = x.get_current_smaller_cut_side_size();
				}
			}

			current_cutter_id = best_cutter_id;
			current_smaller_side_size = cutter_list[current_cutter_id].get_current_smaller_cut_side_size();
		}

		CutterStateDump dump_state()const{
			if(cutter_list.size() != 1)
				throw std::runtime_error("Can only dump the cutter state if a single instance is run");
			return cutter_list[0].dump_state();
		}

		template<class Graph, class SearchAlgorithm, class ScorePierceNode>
		bool advance(const Graph&graph, TemporaryData&tmp, const SearchAlgorithm&search_algo, const ScorePierceNode&score_pierce_node, bool should_skip_non_maximum_sides = true){
			if(graph.node_count() /2 == get_current_smaller_cut_side_size())
				return false;

			int current_cut_size = cutter_list[current_cutter_id].get_current_cut().size();
			for(;;){
				for(int i=0; i<(int)cutter_list.size(); ++i){
					auto x = std::move(cutter_list[i]);
					auto my_score_pierce_node = [&](int x2, int side, bool causes_augmenting_path, int source_dist, int target_dist){
						return score_pierce_node(x2, side, causes_augmenting_path, source_dist, target_dist, i);
					};
					if(x.is_cut_available()){
						if((int)x.get_current_cut().size() == current_cut_size){
#ifdef SLOW_DEBUG
							assert(x.does_next_advance_increase_cut(graph, my_score_pierce_node));
#endif
							if(x.advance(graph, tmp, search_algo, my_score_pierce_node)){
								assert((int)x.get_current_cut().size() > current_cut_size);
								while(!x.does_next_advance_increase_cut(graph, my_score_pierce_node)){
									if(!x.advance(graph, tmp, search_algo, my_score_pierce_node))
										break;
									if(!should_skip_non_maximum_sides)
										break;
								}
							}
						}
					}

					cutter_list[i] = std::move(x);
				}

				int next_cut_size = std::numeric_limits<int>::max();
				for(auto&x:cutter_list)
					if(x.is_cut_available())
						min_to(next_cut_size, (int)x.get_current_cut().size());

				if(next_cut_size == std::numeric_limits<int>::max())
					return false;


				int best_cutter_weight = 0;
				int best_cutter_id = -1;
				for(int i=0; i<(int)cutter_list.size(); ++i){
					if(cutter_list[i].is_cut_available()){
						if(
							(int)cutter_list[i].get_current_cut().size() == next_cut_size &&
							cutter_list[i].get_current_smaller_cut_side_size() > best_cutter_weight
						){
							best_cutter_id = i;
							best_cutter_weight = cutter_list[i].get_current_smaller_cut_side_size();
						}
					}
				}

				assert(best_cutter_id != -1);

				current_cut_size = next_cut_size;

				if(best_cutter_weight <= current_smaller_side_size)
					continue;

				current_cutter_id = best_cutter_id;
				current_smaller_side_size = cutter_list[current_cutter_id].get_current_smaller_cut_side_size();
				return true;
			}
		}

		int get_current_smaller_cut_side_size()const{
			return current_smaller_side_size;
		}

		bool is_on_smaller_side(int x)const{
			return cutter_list[current_cutter_id].is_on_smaller_side(x);
		}

		const std::vector<int>&get_current_cut()const{
			return cutter_list[current_cutter_id].get_current_cut();
		}

		int get_current_cutter_id()const{
			return current_cutter_id;
		}

	private:
		std::vector<DistanceAwareCutter>cutter_list;
		int current_smaller_side_size;
		int current_cutter_id;
	};

	struct PierceNodeScore{
		static constexpr unsigned hash_modulo = ((1u<<31u)-1u);
		unsigned hash_factor, hash_offset;

		PierceNodeScore(const Config& _config): config(_config){
			std::mt19937 gen;
			gen.seed(config.random_seed);
			gen();
			hash_factor = gen() % hash_modulo;
			hash_offset = gen() % hash_modulo;
		}

		Config config;


		int operator()(int x, int side, bool causes_augmenting_path, int source_dist, int target_dist, int /*cutter_id*/)const{

			auto random_number = [&]{
				if(side == BasicCutter::source_side)
					return (hash_factor * (unsigned)(x<<1) + hash_offset) % hash_modulo;
				else
					return (hash_factor * ((unsigned)(x<<1)+1) + hash_offset) % hash_modulo;
			};

			int score;
			switch(config.pierce_rating){
			case Config::PierceRating::max_target_minus_source_hop_dist:
				score = target_dist - source_dist;
				break;
			case Config::PierceRating::max_target_hop_dist:
				score = target_dist;
				break;
			case Config::PierceRating::min_source_hop_dist:
				score = -source_dist;
				break;
			case Config::PierceRating::oldest:
				score = 0;
				break;
			case Config::PierceRating::random:
				score = random_number();
				break;

			default:
				assert(false);
				score = 0;
			}
			switch(config.avoid_augmenting_path){
			case Config::AvoidAugmentingPath::avoid_and_pick_best:
				if(causes_augmenting_path)
					score -= 1000000000;
				break;
			case Config::AvoidAugmentingPath::do_not_avoid:
				break;
			case Config::AvoidAugmentingPath::avoid_and_pick_oldest:
				if(causes_augmenting_path)
					score = -1000000000;
				break;
			case Config::AvoidAugmentingPath::avoid_and_pick_random:
				if(causes_augmenting_path)
					score = random_number() - 1000000000;
				break;
			default:
				assert(false);
				score = 0;
			}
			return score;
		}
	};

	template<class Graph>
	class SimpleCutter{
	public:
		SimpleCutter(const Graph&_graph, Config _config):
			graph(_graph), tmp(_graph.node_count()), config(_config){
		}

		void init(const std::vector<SourceTargetPair>&p, int random_seed){
			DistanceType dist_type;

			if(
				config.pierce_rating == Config::PierceRating::min_source_hop_dist ||
				config.pierce_rating == Config::PierceRating::max_target_hop_dist ||
				config.pierce_rating == Config::PierceRating::max_target_minus_source_hop_dist
			)
				dist_type = DistanceType::hop_distance;
			else
				dist_type = DistanceType::no_distance;

			switch(config.graph_search_algorithm){
			case Config::GraphSearchAlgorithm::pseudo_depth_first_search:
				cutter.init(graph, tmp, PseudoDepthFirstSearch(), PierceNodeScore(config), dist_type, p, random_seed, config.skip_non_maximum_sides == Config::SkipNonMaximumSides::skip);
				break;

			case Config::GraphSearchAlgorithm::breadth_first_search:
				cutter.init(graph, tmp, BreadthFirstSearch(), PierceNodeScore(config), dist_type, p, random_seed, config.skip_non_maximum_sides == Config::SkipNonMaximumSides::skip);
				break;

			case Config::GraphSearchAlgorithm::depth_first_search:
				throw std::runtime_error("depth first search is not yet implemented");
			default:
				assert(false);

			}
		}

		bool advance(){

			switch(config.graph_search_algorithm){
			case Config::GraphSearchAlgorithm::pseudo_depth_first_search:
				return cutter.advance(graph, tmp, PseudoDepthFirstSearch(), PierceNodeScore(config), config.skip_non_maximum_sides == Config::SkipNonMaximumSides::skip);

			case Config::GraphSearchAlgorithm::breadth_first_search:
				return cutter.advance(graph, tmp, BreadthFirstSearch(), PierceNodeScore(config), config.skip_non_maximum_sides == Config::SkipNonMaximumSides::skip);

			case Config::GraphSearchAlgorithm::depth_first_search:
				throw std::runtime_error("depth first search is not yet implemented");
			default:
				assert(false);
				return false;
			}
		}

		CutterStateDump dump_state()const{
			return cutter.dump_state();
		}

		int get_current_smaller_cut_side_size()const{
			return cutter.get_current_smaller_cut_side_size();
		}

		bool is_on_smaller_side(int x)const{
			return cutter.is_on_smaller_side(x);
		}

		const std::vector<int>&get_current_cut()const{
			return cutter.get_current_cut();
		}

		int get_current_cutter_id()const{
			return cutter.get_current_cutter_id();
		}

	private:
		const Graph&graph;
		TemporaryData tmp;
		MultiCutter cutter;
		Config config;
	};

	template<class Graph>
	SimpleCutter<Graph> make_simple_cutter(const Graph&graph, Config config){
		return SimpleCutter<Graph>(graph, config);
	}

	inline std::vector<SourceTargetPair>select_random_source_target_pairs(int node_count, int cutter_count, int seed){
		std::vector<SourceTargetPair>p(cutter_count);
		std::mt19937 rng(seed);
		for(auto&x:p){
			do{
				// We do not use std::uniform_distribution because it produces difference results for different compilers
				x.source = rng()%node_count;
				x.target = rng()%node_count;
			}while(x.source == x.target);
		}
		return p;
	}
}

#endif

