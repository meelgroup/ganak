#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "id_func_traits.h"
#include "array_id_func.h"
#include <cassert>

template<class IDIDFunc>
typename std::enable_if<
	is_id_id_func<IDIDFunc>::value,
	ArrayIDFunc<int>
>::type compute_histogram(const IDIDFunc&f){
	ArrayIDFunc<int>h(f.image_count());
	h.fill(0);
	for(int i=0; i<f.preimage_count(); ++i)
		++h[f(i)];
	return h; // NVRO
}

template<class IDFunc>
typename std::enable_if<
	is_id_func<IDFunc>::value,
	int
>::type max_histogram_id(const IDFunc&h){
	assert(h.preimage_count() != 0);

	typename id_func_image_type<IDFunc>::type max_element = h(0);
	int max_id = 0;
	for(int i=1; i<h.preimage_count(); ++i){
		auto element = h(i);
		if(max_element < element){
			max_element = element;
			max_id = i;
		}
	}
	#ifndef NDEBUG
	for(int i=0; i<h.preimage_count(); ++i)
		assert(h(i) <= h(max_id));
	#endif
	return max_id;
}

template<class IDFunc>
typename std::enable_if<
	is_id_func<IDFunc>::value,
	int
>::type min_histogram_id(const IDFunc&h){
	assert(h.preimage_count() != 0);

	typename id_func_image_type<IDFunc>::type min_element = h(0);
	int min_id = 0;
	for(int i=1; i<h.preimage_count(); ++i){
		auto element = h(i);
		if(min_element > element){
			min_element = element;
			min_id = i;
		}
	}
	#ifndef NDEBUG
	for(int i=0; i<h.preimage_count(); ++i)
		assert(h(i) >= h(min_id));
	#endif
	return min_id;
}

#endif
