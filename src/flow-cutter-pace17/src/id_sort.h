#ifndef ID_SORT_H
#define ID_SORT_H

#include "array_id_func.h"
#include <cassert>

template<class InIter, class OutIter, class GetID>
void stable_sort_copy_by_id(
	InIter in_begin, InIter in_end, 
	OutIter out_iter, 
	int id_count, const GetID&get_id
){
	ArrayIDFunc<int>pos(id_count);
	pos.fill(0);
	for(InIter i=in_begin; i!=in_end; ++i)
		++pos[get_id(*i)];

	int sum = 0;
	for(int i=0; i<id_count; ++i){
		int tmp = pos[i];
		pos[i] = sum;
		sum += tmp;
	}

	for(InIter i=in_begin; i!=in_end; ++i)
		*(out_iter+pos[get_id(*i)]++) = *i;
}

template<class InIter, class OutIter, class GetID>
void stable_sort_copy_by_id(
	InIter in_begin, InIter in_end, 
	OutIter out_iter, 
	const GetID&get_id
){
	stable_sort_copy_by_id(in_begin, in_end, out_iter, get_id.image_count(), get_id);
}

#endif

