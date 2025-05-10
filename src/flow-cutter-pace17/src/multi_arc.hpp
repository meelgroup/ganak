#ifndef MULTI_ARC_H
#define MULTI_ARC_H

#include "id_func.hpp"
#include "id_sort.hpp"
#include "array_id_func.hpp"
#include "sort_arc.hpp"

template<class Tail, class Head>
BitIDFunc identify_non_multi_arcs(const Tail&tail, const Head&head){
	const int arc_count = tail.preimage_count();
	auto arc_list = sort_arcs_first_by_tail_second_by_head(tail, head);
	BitIDFunc is_non_multi_arc(arc_count);
	if(arc_count > 0){
		is_non_multi_arc.set(arc_list[0], true);
		for(int i=1; i<arc_count; ++i)
			is_non_multi_arc.set(arc_list[i],
				   head(arc_list[i]) != head(arc_list[i-1])
				|| tail(arc_list[i]) != tail(arc_list[i-1])
			);
	}
	return is_non_multi_arc; // NRVO
}


template<class Tail, class Head>
bool is_symmetric(const Tail&tail, const Head&head){
	const int arc_count = tail.preimage_count();
	auto forward_arc_list = sort_arcs_first_by_tail_second_by_head(tail, head);
	auto backward_arc_list = sort_arcs_first_by_tail_second_by_head(head, tail);

	for(int i=0; i<arc_count; ++i){
		if(tail(forward_arc_list(i)) != head(backward_arc_list(i)) || head(forward_arc_list(i)) != tail(backward_arc_list(i)))
			return false;
	}
	return true;
}

template<class Tail, class Head>
class SymmetricHead{
public:
	SymmetricHead(Tail _tail, Head _head):
		tail(std::move(_tail)), head(std::move(_head)){}

	int preimage_count()const{ return 2*tail.preimage_count(); }
	int image_count()const {return tail.image_count(); }

	int operator()(int x)const{
		if(x < tail.preimage_count())
			return head(x);
		else
			return tail(x - tail.preimage_count());
	}
private:
	Tail tail;
	Head head;
};

template<class Tail, class Head>
SymmetricHead<Tail, Head>make_symmetric_head(Tail tail, Head head){
	return {std::move(tail), std::move(head)};
}

template<class Tail, class Head>
class SymmetricTail{
public:
	SymmetricTail(Tail _tail, Head _head):
		tail(std::move(_tail)), head(std::move(_head)){}

	int preimage_count()const{ return 2*tail.preimage_count(); }
	int image_count()const {return tail.image_count(); }

	int operator()(int x)const{
		if(x < tail.preimage_count())
			return tail(x);
		else
			return head(x - tail.preimage_count());
	}
private:
	Tail tail;
	Head head;
};

template<class Tail, class Head>
SymmetricTail<Tail, Head>make_symmetric_tail(Tail tail, Head head){
	return {std::move(tail), std::move(head)};
}

template<class Tail, class Head>
bool has_multi_arcs(const Tail&tail, const Head&head){
	return count_true(identify_non_multi_arcs(tail, head)) != tail.preimage_count();
}

template<class Tail, class Head>
bool is_loop_free(const Tail&tail, const Head&head){
	for(int i=0; i<tail.preimage_count(); ++i)
		if(tail(i) == head(i))
			return false;
	return true;
}

#endif

