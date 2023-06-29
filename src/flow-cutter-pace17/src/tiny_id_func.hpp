#ifndef TINY_INT_ID_FUNC_H
#define TINY_INT_ID_FUNC_H

#include <cstdint>
#include <utility>
#include "id_func.hpp"
#include "array_id_func.hpp"

template<int bit_count>
struct TinyIntIDFunc{
private:
	static constexpr int entry_count_per_uint64 = 64/bit_count;
	static constexpr std::uint64_t entry_mask = (std::uint64_t(1)<<bit_count) - std::uint64_t(1);

public:
	static_assert(1<=bit_count && bit_count <= 64, "an integer with more than 64 bits is not tiny");

	TinyIntIDFunc():preimage_(0){}

	explicit TinyIntIDFunc(int preimage)
		:preimage_(preimage), data_((preimage + entry_count_per_uint64 - 1) / entry_count_per_uint64){}

	int preimage_count()const{
		return preimage_;
	}

	template<class IDFunc>
	TinyIntIDFunc(const IDFunc&other, [[maybe_unused]]
			typename std::enable_if<is_id_func<IDFunc>::value, void>::type* dummy=0)
		:preimage_(other.preimage_count()), data_((other.preimage_count() + entry_count_per_uint64 - 1) / entry_count_per_uint64){
		for(int i=0; i<preimage_count(); ++i)
			set(i, other(i));
	}

	std::uint64_t operator()(int id)const{
		assert(0 <= id && id < preimage_ && "id out of bounds");

		int index = id / entry_count_per_uint64;
		int offset = (id % entry_count_per_uint64)*bit_count;
		return (data_[index] >> offset) & entry_mask;
	}

	void set(int id, std::uint64_t value){
		assert(0 <= id && id < preimage_ && "id out of bounds");
		assert(value <= entry_mask && "value out of bounds");

		int index = id / entry_count_per_uint64;
		int offset = (id % entry_count_per_uint64)*bit_count;

		data_[index] ^= ((((data_[index] >> offset) & entry_mask) ^ value) << offset);
	}

	void fill(std::uint64_t value){
		assert(value <= entry_mask && "value out of bounds");

		if(bit_count == 1){
			if(value == false)
				data_.fill(0);
			else
				data_.fill((std::uint64_t)-1);
		}else if(value == 0){
			data_.fill(0);
		}else{
			std::uint64_t pattern = value;
			int shift = bit_count;
			while(shift < 64){
				pattern |= pattern << shift;
				shift <<= 1;
			}
			data_.fill(pattern);
		}
	}

	std::uint64_t move(int id){
		return operator()(id);
	}

	void swap(TinyIntIDFunc&other)noexcept{
		std::swap(preimage_, other.preimage_);
		data_.swap(other.data_);
	}

	template<class IDFunc>
	TinyIntIDFunc operator=(const typename std::enable_if<is_id_func<IDFunc>::value, IDFunc>::type & other){
		TinyIntIDFunc(other).swap(*this);
		return *this;
	}


	int preimage_;
	ArrayIDFunc<std::uint64_t> data_;
};



typedef TinyIntIDFunc<1> BitIDFunc;

inline BitIDFunc operator~(BitIDFunc f){
	for(int i=0; i<f.data_.preimage_count(); ++i)
		f.data_[i] = ~f.data_[i];
	return f;
}

inline BitIDFunc&operator^=(BitIDFunc&l, const BitIDFunc&r){
	assert(l.preimage_count() == r.preimage_count());

	for(int i=0; i<l.data_.preimage_count(); ++i)
		l.data_[i] ^= r.data_[i];
	return l;
}

inline BitIDFunc&operator&=(BitIDFunc&l, const BitIDFunc&r){
	assert(l.preimage_count() == r.preimage_count());

	for(int i=0; i<l.data_.preimage_count(); ++i)
		l.data_[i] &= r.data_[i];
	return l;
}

inline BitIDFunc&operator|=(BitIDFunc&l, const BitIDFunc&r){
	assert(l.preimage_count() == r.preimage_count());

	for(int i=0; i<l.data_.preimage_count(); ++i)
		l.data_[i] |= r.data_[i];
	return l;
}

inline BitIDFunc&inplace_and_not(BitIDFunc&l, const BitIDFunc&r){
	assert(l.preimage_count() == r.preimage_count());

	for(int i=0; i<l.data_.preimage_count(); ++i)
		l.data_[i] &= ~r.data_[i];
	return l;
}

inline BitIDFunc operator^(BitIDFunc l, const BitIDFunc&r){
	l ^= r;
	return l;
}

inline BitIDFunc operator|(BitIDFunc l, const BitIDFunc&r){
	l |= r;
	return l;
}

inline BitIDFunc operator&(BitIDFunc l, const BitIDFunc&r){
	l &= r;
	return l;
}

inline BitIDFunc and_not(BitIDFunc l, const BitIDFunc&r){
	inplace_and_not(l, r);
	return l;
}

#endif
