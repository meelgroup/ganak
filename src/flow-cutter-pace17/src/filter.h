#ifndef FILTER_H
#define FILTER_H

#include "tiny_id_func.h"

template<class Pred>
int count_true(const Pred&p){
	int sum = 0;
	for(int i=0; i<p.preimage_count(); ++i)
		if(p(i))
			++sum;
	return sum;
}

template<class Pred, class IDFunc>
typename std::enable_if<
	is_only_id_func<IDFunc>::value,
	ArrayIDFunc<typename id_func_image_type<IDFunc>::type>
>::type keep_if(const Pred&p, int new_preimage_count, const IDFunc&f){
	assert(p.preimage_count() == f.preimage_count());
	assert(new_preimage_count == count_true(p));
	
	ArrayIDFunc<typename id_func_image_type<IDFunc>::type>result(new_preimage_count);

	int out = 0;
	for(int in=0; in<f.preimage_count(); ++in)
		if(p(in))
			result[out++] = f(in);
	assert(out == new_preimage_count);
	return result; // NVRO
}

template<class Pred, class IDFunc>
typename std::enable_if<
	is_id_id_func<IDFunc>::value,
	ArrayIDIDFunc
>::type keep_if(const Pred&p, int new_preimage_count, const IDFunc&f){
	assert(p.preimage_count() == f.preimage_count());
	assert(new_preimage_count == count_true(p));
	
	ArrayIDIDFunc result(new_preimage_count, f.image_count());

	int out = 0;
	for(int in=0; in<f.preimage_count(); ++in)
		if(p(in))
			result[out++] = f(in);
	assert(out == new_preimage_count);
	return result; // NRVO
}

template<class Pred>
ArrayIDIDFunc compute_keep_function(const Pred&pred, int new_image_count){
	ArrayIDIDFunc f(pred.preimage_count(), new_image_count);
	int next_id = 0;
	for(int i=0; i<pred.preimage_count(); ++i)
		if(pred(i))
			f[i] = next_id++;
		else
			f[i] = -1;
	return f; // NRVO
}

template<class Pred>
ArrayIDIDFunc compute_inverse_keep_function(const Pred&pred, int new_image_count){
	ArrayIDIDFunc f(new_image_count, pred.preimage_count());
	int next_id = 0;
	for(int i=0; i<pred.preimage_count(); ++i)
		if(pred(i))
			f[next_id++] = i;
	assert(new_image_count == next_id);
	return f; // NRVO
}
#endif
