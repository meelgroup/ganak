#ifndef MIN_MAX_H
#define MIN_MAX_H

#include <utility>
#include <algorithm>
#include "id_func_traits.h"

template<class T>
void min_to(T&x, T y){
	if(y < x)
		x = std::move(y);
}

template<class T>
void max_to(T&x, T y){
	if(y > x)
		x = std::move(y);
}

template<class T>
void sort_ref_args(T&x, T&y){
	if(y < x)
		std::swap(x,y);
}

template<class F>
typename id_func_image_type<F>::type min_over_id_func(const F&f){
	assert(f.preimage_count() != 0);
	typename id_func_image_type<F>::type result = f(0);
	for(int i=1; i<f.preimage_count(); ++i)
		min_to(result, f(i));
	return result; //NVRO
}

template<class F>
typename id_func_image_type<F>::type max_over_id_func(const F&f){
	assert(f.preimage_count() != 0);
	typename id_func_image_type<F>::type result = f(0);
	for(int i=1; i<f.preimage_count(); ++i)
		max_to(result, f(i));
	return result; //NVRO
}

template<class F>
int min_preimage_over_id_func(const F&f){
	assert(f.preimage_count() != 0);
	int preimage = 0;
	typename id_func_image_type<F>::type m = f(0);
	for(int i=1; i<f.preimage_count(); ++i){
		auto x = f(i);
		if(x < m){
			m = std::move(x);
			preimage = i;
		}
	}
	return preimage; //NVRO
}

template<class F>
int max_preimage_over_id_func(const F&f){
	assert(f.preimage_count() != 0);
	int preimage = 0;
	typename id_func_image_type<F>::type m = f(0);
	for(int i=1; i<f.preimage_count(); ++i){
		auto x = f(i);
		if(m < x){
			m = std::move(x);
			preimage = i;
		}
	}
	return preimage; //NVRO
}

#endif

