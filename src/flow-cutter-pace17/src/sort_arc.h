#ifndef SORT_ARC_H
#define SORT_ARC_H

#include "id_sort.h"
#include "array_id_func.h"
#include "permutation.h"
#include "count_range.h"
#include <cassert>

template<class Tail, class Head>
ArrayIDIDFunc sort_arcs_first_by_tail_second_by_head(const Tail&tail, const Head&head){
	assert(tail.preimage_count() == head.preimage_count());
	assert(tail.image_count() == head.image_count());

	const int arc_count = tail.preimage_count();

	ArrayIDIDFunc 
		x(arc_count, arc_count),
		y(arc_count, arc_count);

	stable_sort_copy_by_id(
		CountIterator{0}, CountIterator{arc_count}, 
		y.begin(),
		head.image_count(),
		head
	);
	stable_sort_copy_by_id(
		y.begin(), y.end(), 
		x.begin(),
		tail.image_count(),
		tail
	);

	return x; //NVRO
}

#endif

