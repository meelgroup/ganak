#ifndef NODE_FLOW_CUTTER_H
#define NODE_FLOW_CUTTER_H

#include "flow_cutter.hpp"
#include "range.hpp"

namespace flow_cutter{

	//!
	//! An expanded graph is a virtually expanded version of a symmetric input graph.
	//! Expanded means that every original node v is subdivided into two nodes v_in and v_out.
	//!
	//! For each arc x->y two new inter arcs x_out -> y_in and y_in -> x_out are created
	//! For each node x two new intra arcs x_out -> x_in and x_in -> x_out are created
	//! Every in to out arc has capacity 1 and every out to in arc has capacity 0.
	//!

	namespace expanded_graph{
		inline int expanded_node_count(int original_node_count){ return 2*original_node_count; }
		inline int expanded_arc_count(int original_node_count, int original_arc_count){ return 2*(original_node_count+original_arc_count); }
		inline int expanded_node_to_original_node(int x){ return x/2; }
		inline int original_node_to_expanded_node(int x, bool is_out){ return 2*x+is_out; }
		inline bool get_expanded_node_out_flag(int x){ return x&1; }
		inline bool is_expanded_intra_arc(int x, int original_arc_count){ return x >= 2*original_arc_count; }
		inline bool is_expanded_inter_arc(int x, int original_arc_count){ return x < 2*original_arc_count; }
		inline bool get_expanded_arc_tail_out_flag(int x){return x&1;}
		inline int expanded_inter_arc_to_original_arc(int x, int original_arc_count){ (void)original_arc_count; return x/2; }
		inline int original_arc_to_expanded_inter_arc(int x, bool tail_out_flag, int/* original_arc_count*/){ return 2*x+tail_out_flag; }
		inline int expanded_intra_arc_to_original_node(int x, int original_arc_count){ (void)original_arc_count; return x/2-original_arc_count; }
		inline int original_node_to_expanded_intra_arc(int x, bool tail_out_flag, int original_arc_count){ return 2*(original_arc_count+x)+tail_out_flag; }

		template<class OriginalTail>
		struct Tail{
			int original_node_count, original_arc_count;
			OriginalTail original_tail;

			int preimage_count()const{return expanded_arc_count(original_node_count, original_arc_count);}
			int image_count()const{return expanded_node_count(original_node_count);}

			int operator()(int a)const{
				if(is_expanded_intra_arc(a, original_arc_count)){
					return original_node_to_expanded_node(expanded_intra_arc_to_original_node(a, original_arc_count), get_expanded_arc_tail_out_flag(a));
				}else{
					return original_node_to_expanded_node(original_tail(expanded_inter_arc_to_original_arc(a, original_arc_count)), get_expanded_arc_tail_out_flag(a));
				}
			}
		};

		template<class OriginalTail>
		Tail<OriginalTail>tail(int original_node_count, int original_arc_count, OriginalTail original_tail){
			return {original_node_count, original_arc_count, std::move(original_tail)};
		}

		template<class OriginalHead>
		struct Head{
			int original_node_count, original_arc_count;
			OriginalHead original_head;

			int preimage_count()const{return expanded_arc_count(original_node_count, original_arc_count);}
			int image_count()const{return expanded_node_count(original_node_count);}

			int operator()(int a)const{
				if(is_expanded_intra_arc(a, original_arc_count)){
					return original_node_to_expanded_node(expanded_intra_arc_to_original_node(a, original_arc_count), !get_expanded_arc_tail_out_flag(a));
				}else{
					return original_node_to_expanded_node(original_head(expanded_inter_arc_to_original_arc(a, original_arc_count)), !get_expanded_arc_tail_out_flag(a));
				}
			}
		};

		template<class OriginalHead>
		Head<OriginalHead>head(int original_node_count, int original_arc_count, OriginalHead original_head){
			return {original_node_count, original_arc_count, std::move(original_head)};
		}

		template<class OriginalBackArc>
		struct BackArc{
			int original_node_count, original_arc_count;
			OriginalBackArc original_back_arc;

			int preimage_count()const{return expanded_arc_count(original_node_count, original_arc_count);}
			int image_count()const{return expanded_arc_count(original_node_count, original_arc_count);}

			int operator()(int a)const{
				if(is_expanded_intra_arc(a, original_arc_count)){
					return original_node_to_expanded_intra_arc(expanded_intra_arc_to_original_node(a, original_arc_count), !get_expanded_arc_tail_out_flag(a), original_arc_count);
				}else{
					return original_arc_to_expanded_inter_arc(original_back_arc(expanded_inter_arc_to_original_arc(a, original_arc_count)), !get_expanded_arc_tail_out_flag(a), original_arc_count);
				}
			}
		};

		template<class OriginalBackArc>
		BackArc<OriginalBackArc>back_arc(int original_node_count, int original_arc_count, OriginalBackArc original_back_arc){
			return {original_node_count, original_arc_count, std::move(original_back_arc)};
		}

		struct Capacity{
			int original_node_count, original_arc_count;

			int preimage_count()const{return expanded_arc_count(original_node_count, original_arc_count);}

			int operator()(int a)const{
				if(is_expanded_intra_arc(a, original_arc_count)){
					return !get_expanded_arc_tail_out_flag(a);
				}else{
					return get_expanded_arc_tail_out_flag(a);
				}
			}
		};

		inline Capacity capacity(int original_node_count, int original_arc_count){
			return {original_node_count, original_arc_count};
		}

		template<class OriginalOutArc>
		struct OutArcIter{

			typedef typename std::decay<decltype(std::begin(std::declval<const OriginalOutArc>()(0)))>::type OriginalOutArcIter;

			typedef int value_type;
			typedef int difference_type;
			typedef const int* pointer;
			typedef const int& reference;
			typedef std::forward_iterator_tag iterator_category;

			OutArcIter(){}

			OutArcIter(int _intra_arc, bool _node_out_flag, OriginalOutArcIter _base_iter):
				intra_arc(_intra_arc), node_out_flag(_node_out_flag), base_iter(_base_iter){}

			OutArcIter&operator++(){
				if(intra_arc != -1)
					intra_arc = -1;
				else
					++base_iter;
				return *this;
			}

			OutArcIter operator++(int) {
				OutArcIter tmp(*this);
				operator++();
				return tmp;
			}

			int operator*()const{
				if(intra_arc != -1)
					return intra_arc;
				else
					return 2*(*base_iter) + node_out_flag;
			}

			int dummy; // To whom ever defined the ackward operator-> semantics: Skrew you!
			const int*operator->() const {
				dummy = *this;
				return &dummy;
			}

			int intra_arc;
			bool node_out_flag;
			OriginalOutArcIter base_iter;

			friend bool operator==(OutArcIter l, OutArcIter r){
				return l.base_iter == r.base_iter && l.intra_arc == r.intra_arc && l.node_out_flag == r.node_out_flag;
			}

			friend bool operator!=(OutArcIter l, OutArcIter r){
				return !(l == r);
			}
		};

		template<class OriginalOutArc>
		struct OutArc{
			int original_node_count, original_arc_count;
			OriginalOutArc original_out_arc;

			int preimage_count()const{return expanded_node_count(original_node_count);}

			Range<OutArcIter<OriginalOutArc>> operator()(int x)const{
				int original_x = expanded_node_to_original_node(x);
				bool out_flag = get_expanded_node_out_flag(x);
				auto r = original_out_arc(original_x);
				return Range<OutArcIter<OriginalOutArc>>{
					OutArcIter<OriginalOutArc>{original_node_to_expanded_intra_arc(original_x, out_flag, original_arc_count), out_flag, std::begin(r)},
					OutArcIter<OriginalOutArc>{-1, out_flag, std::end(r)}
				};
			}
		};

		template<class OriginalOutArc>
		OutArc<OriginalOutArc>out_arc(int original_node_count, int original_arc_count, OriginalOutArc original_out_arc){
			return {original_node_count, original_arc_count, std::move(original_out_arc)};
		}


		template<class Tail, class Head, class BackArc, class OutArc>
		Graph<
			expanded_graph::Tail<Tail>,
			expanded_graph::Head<Head>,
			expanded_graph::BackArc<BackArc>,
			expanded_graph::Capacity,
			expanded_graph::OutArc<OutArc>
		>
		make_graph(Tail tail, Head head, BackArc back_arc, OutArc out_arc){
			int node_count = tail.image_count(), arc_count = tail.preimage_count();
			return{
				expanded_graph::tail(node_count, arc_count, std::move(tail)),
				expanded_graph::head(node_count, arc_count, std::move(head)),
				expanded_graph::back_arc(node_count, arc_count, std::move(back_arc)),
				expanded_graph::capacity(node_count, arc_count),
				expanded_graph::out_arc(node_count, arc_count, std::move(out_arc))
			};
		}

		struct MixedCut{
			std::vector<int>arcs, nodes;
		};

		inline
		MixedCut expanded_cut_to_original_mixed_cut(const std::vector<int>&expanded_cut, int original_arc_count){
			MixedCut original_cut;

			for(auto x:expanded_cut){
				if(is_expanded_inter_arc(x, original_arc_count)){
					original_cut.arcs.push_back(expanded_inter_arc_to_original_arc(x, original_arc_count));
				}else{
					original_cut.nodes.push_back(expanded_intra_arc_to_original_node(x, original_arc_count));
				}
			}

			return original_cut; // NRVO
		}

		struct Separator{
			std::vector<int>sep;
			int small_side_size;
		};

		template<class Tail, class Head, class FlowCutter>
		Separator extract_original_separator(const Tail&tail, const Head&head, const FlowCutter&cutter){
			int original_node_count = tail.image_count();
			int original_arc_count = tail.preimage_count();

			Separator sep;

			for(auto x:cutter.get_current_cut()){
				if(is_expanded_intra_arc(x, original_arc_count)){
					sep.sep.push_back(expanded_intra_arc_to_original_node(x, original_arc_count));
				}
			}

			int left_side_size = (cutter.get_current_smaller_cut_side_size()-sep.sep.size())/2;
			int right_side_size = original_node_count - sep.sep.size() - left_side_size;

			auto is_original_node_left = [&](int x){
				return cutter.is_on_smaller_side(original_node_to_expanded_node(x, true));
			};

			for(auto x:cutter.get_current_cut()){
				if(is_expanded_inter_arc(x, original_arc_count)){
					auto lr = expanded_inter_arc_to_original_arc(x, original_arc_count);
					auto l = tail(lr), r = head(lr);
					if(is_original_node_left(r))
						std::swap(l, r);

					if(left_side_size > right_side_size){
						sep.sep.push_back(l);
						--left_side_size;
					}else{
						sep.sep.push_back(r);
						--right_side_size;
					}
				}
			}

			sep.small_side_size = std::min(left_side_size, right_side_size);

			std::sort(sep.sep.begin(), sep.sep.end());
			sep.sep.erase(std::unique(sep.sep.begin(), sep.sep.end()), sep.sep.end());

			return sep; // NVRO
		}

		inline
		std::vector<SourceTargetPair>expand_source_target_pair_list(std::vector<SourceTargetPair>p){
			for(auto&x:p){
				x.source = original_node_to_expanded_node(x.source, false);
				x.target = original_node_to_expanded_node(x.target, true);
			}
			return p;
		}
	}


}

#endif
