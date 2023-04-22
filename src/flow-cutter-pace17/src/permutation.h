#ifndef PERMUTATION_H
#define PERMUTATION_H

#include "tiny_id_func.h"

template<class IDIDFunc>
bool is_permutation(const IDIDFunc&f){
	if(f.preimage_count() != f.image_count())
		return false;

	int id_count = f.preimage_count();

	BitIDFunc already_seen(id_count);
	already_seen.fill(false);
	for(int i=0; i<id_count; ++i){
		int x = f(i);
		if(x < 0 || x >= id_count)
			return false;
		if(already_seen(x))
			return false;
		already_seen.set(x, true);
	}
	return true;
}



template<class IDIDFunc>
ArrayIDIDFunc inverse_permutation(const IDIDFunc&f){
	assert(is_permutation(f));

	int id_count = f.preimage_count();

	ArrayIDIDFunc inv_f(id_count, id_count);
	for(int i=0; i<id_count; ++i)
		inv_f[f(i)] = i;
	return inv_f; // NVRO
}

inline
ArrayIDIDFunc identity_permutation(int id_count){
	ArrayIDIDFunc f(id_count, id_count);
	for(int i=0; i<id_count; ++i)
		f[i] = i;
	return f; // NVRO
}

#endif
