#ifndef ID_MULTI_FUNC_H
#define ID_MULTI_FUNC_H

#include "array_id_func.h"
#include "count_range.h"
#include "range.h"
#include "chain.h"
#include <cassert>

struct RangeIDIDMultiFunc{
	int preimage_count()const{ return range_begin.preimage_count()-1; }
	int image_count()const{ return range_begin(preimage_count()); }
	
	CountRange operator()(int id)const{
		assert(0 <= id && id < preimage_count() && "id out of bounds");
		return count_range(range_begin(id), range_begin(id+1));
	}

	ArrayIDFunc<int>range_begin;
};


template<class T>
struct ArrayIDMultiFunc{
	int preimage_count()const{ return preimage_to_intermediate.preimage_count(); }	

	Range<T*>operator()(int id){
		assert(0 <= id && id < preimage_count() && "id out of bounds");
		return {
			intermediate_to_image.begin() + *std::begin(preimage_to_intermediate(id)), 
			intermediate_to_image.begin() + *std::end(preimage_to_intermediate(id))
		};
	}

	Range<const T*>operator()(int id)const{
		assert(0 <= id && id < preimage_count() && "id out of bounds");
		return {
			intermediate_to_image.begin() + *std::begin(preimage_to_intermediate(id)), 
			intermediate_to_image.begin() + *std::end(preimage_to_intermediate(id))
		};
	}

	RangeIDIDMultiFunc preimage_to_intermediate;
	ArrayIDFunc<T>intermediate_to_image;
};


struct ArrayIDIDMultiFunc{
	int image_count()const{ return intermediate_to_image.image_count(); }
	int preimage_count()const{ return preimage_to_intermediate.preimage_count(); }	

	Range<int*>operator()(int id){
		assert(0 <= id && id < preimage_count() && "id out of bounds");
		return {
			intermediate_to_image.begin() + *std::begin(preimage_to_intermediate(id)), 
			intermediate_to_image.begin() + *std::end(preimage_to_intermediate(id))
		};
	}

	Range<const int*>operator()(int id)const{
		assert(0 <= id && id < preimage_count() && "id out of bounds");
		return {
			intermediate_to_image.begin() + *std::begin(preimage_to_intermediate(id)), 
			intermediate_to_image.begin() + *std::end(preimage_to_intermediate(id))
		};
	}

	RangeIDIDMultiFunc preimage_to_intermediate;
	ArrayIDIDFunc intermediate_to_image;
};

//! Inverts an id-id function f. In this context we have two ID types: preimage IDs
//! and image IDs. f maps preimage IDs onto image IDs. This function computes a
//! id-id multi function g that maps image IDs onto preimage ID ranges.
//! g(x) is a range of all y such that f(y) = x ordered increasing by y.
template<class IDIDFunc>
ArrayIDIDMultiFunc invert_id_id_func(const IDIDFunc&f){
	ArrayIDIDMultiFunc g = {
		RangeIDIDMultiFunc{
			ArrayIDFunc<int>{f.image_count()+1}
		},
		ArrayIDIDFunc{f.preimage_count(), f.preimage_count()}
	};

	auto&begin = g.preimage_to_intermediate.range_begin;
	begin.fill(0);

	for(int i=0; i<f.preimage_count(); ++i)
		++begin[f(i)];


	int sum = 0;
	for(int i=0; i<=f.image_count(); ++i){
		int tmp = begin[i];
		begin[i] = sum;
		sum += tmp;
	}

	assert(f.preimage_count() == sum);

	auto&p = g.intermediate_to_image;
	p.fill(-1);
	
	for(int i=0; i<f.preimage_count(); ++i){
		p[begin[f(i)]++] = i;
	}

	for(int i=0; i<f.preimage_count(); ++i)
		--begin[f(i)];

	return g; // NRVO
}

template<class IDIDFunc>
RangeIDIDMultiFunc invert_sorted_id_id_func(const IDIDFunc&f){
	assert(std::is_sorted(f.begin(), f.end()) && "f is not sorted");

	RangeIDIDMultiFunc h = {ArrayIDFunc<int>{f.image_count()+1}};

	auto&begin = h.range_begin;
	begin.fill(0);

	for(int i=0; i<f.preimage_count(); ++i)
		++begin[f(i)];


	int sum = 0;
	for(int i=0; i<=f.image_count(); ++i){
		int tmp = begin[i];
		begin[i] = sum;
		sum += tmp;
	}

	assert(f.preimage_count() == sum);
	return h; // NVRO
}

template<class Tail, class Head>
ArrayIDIDMultiFunc compute_successor_function(const Tail&tail, const Head&head){
	auto x = invert_id_id_func(tail);
	x.intermediate_to_image = chain(x.intermediate_to_image, head);
	return x; // NRVO
}

#endif

