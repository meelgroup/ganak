#ifndef CHAIN_H
#define CHAIN_H

#include <type_traits>
#include "id_func.hpp"
#include "array_id_func.hpp"

// chain(IDIDFunc, IDFunc)
template<class L, class R>
typename std::enable_if<
	   is_id_id_func<L>::value
	&& is_id_func<R>::value
	&& !is_id_id_func<R>::value,
	ArrayIDFunc<typename id_func_image_type<R>::type>
>::type
chain(const L&l, const R&r){
	ArrayIDFunc<typename id_func_image_type<R>::type>result(l.preimage_count());
	for(int i=0; i<l.preimage_count(); ++i)
		result[i] = r(l(i));
	return result; // NVRO
}

// chain(IDIDFunc, IDIDFunc)
template<class L, class R>
typename std::enable_if<
	   is_mutable_id_id_func<L>::value
	&& is_id_id_func<R>::value,
	L
>::type
chain(L l, const R&r){
	assert(l.image_count() == r.preimage_count());
	for(int i=0; i<l.preimage_count(); ++i)
		l.set(i, r(l(i)));
	l.set_image_count(r.image_count());
	return l;
}

template<class L, class R>
typename std::enable_if<
	   is_id_id_func<L>::value
	&& !is_mutable_id_id_func<L>::value
	&& is_id_id_func<R>::value,
	ArrayIDIDFunc
>::type
chain(const L&l, const R&r){
	assert(l.image_count() == r.preimage_count());
	ArrayIDIDFunc result(l.preimage_count(), r.image_count());
	for(int i=0; i<l.preimage_count(); ++i)
		result[i] = r(l(i));
	return result; // NVRO
}

#endif

