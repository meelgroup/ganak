#ifndef BACK_ARC_H
#define BACK_ARC_H

#include "array_id_func.h"
#include "id_sort.h"

// Input graph must be symmetric
template<class Tail, class Head>
ArrayIDIDFunc compute_back_arc_permutation(const Tail&tail, const Head&head){
	
	const int arc_count = head.preimage_count();
	const int node_count = head.image_count();

	struct D{
		int tail, head, arc_id;
	};

	ArrayIDFunc<D>arc_list(arc_count), tmp(arc_count);
	for(int i=0; i<arc_count; ++i){
		arc_list[i].tail = tail(i);
		arc_list[i].head = head(i);
		if(arc_list[i].tail > arc_list[i].head)
			std::swap(arc_list[i].tail, arc_list[i].head);
		arc_list[i].arc_id = i;
	}

	stable_sort_copy_by_id(
		std::begin(arc_list), std::end(arc_list), 
		std::begin(tmp),
		node_count,
		[](D d){return d.head;}
	);
	stable_sort_copy_by_id(
		std::begin(tmp), std::end(tmp), 
		std::begin(arc_list),
		node_count,
		[](D d){return d.tail;}
	);

	ArrayIDIDFunc back_arc(head.preimage_count(), head.preimage_count());

	for(int i=0; i<arc_count; i+=2){
		back_arc[arc_list(i).arc_id] = arc_list(i+1).arc_id;
		back_arc[arc_list(i+1).arc_id] = arc_list(i).arc_id;
	}
	
	return back_arc;
}

#endif

