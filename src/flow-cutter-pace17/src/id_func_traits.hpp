#ifndef ID_FUNC_TRAITS_H
#define ID_FUNC_TRAITS_H

#include <type_traits>

template<class IDFunc>
struct id_func_image_type{
	typedef typename std::decay<decltype(std::declval<IDFunc>()(0))>::type type;
};

#define MAKE_TYPED_HAS_TRAIT(HAS, TYPE, EXPR)\
	template<class FT>std::false_type HAS##_impl2(...);\
	template<class FT>std::true_type HAS##_impl2(typename std::decay<decltype(EXPR)>::type*);\
	template<class FT>auto HAS##_impl1()->decltype(HAS##_impl2<FT>(static_cast<TYPE*>(nullptr)));\
	template<class FT>struct HAS : std::enable_if<true, decltype(HAS##_impl1<FT>())>::type{};

#define MAKE_UNTYPED_HAS_TRAIT(HAS, EXPR)\
	template<class FT>std::false_type HAS##_impl2(...);\
	template<class FT>std::true_type HAS##_impl2(typename std::decay<decltype(EXPR)>::type*);\
	template<class FT>auto HAS##_impl1()->decltype(HAS##_impl2<FT>(nullptr));\
	template<class FT>struct HAS : std::enable_if<true, decltype(HAS##_impl1<FT>())>::type{};

#define F std::declval<FT>()
#define IT typename id_func_image_type<FT>::type
#define I std::declval<IT>()

// F is an instance of the function, FT is the function type, I is an instance of the function's image type and IT is the image type

MAKE_TYPED_HAS_TRAIT  (has_preimage_count,        int,  F.preimage_count()          )
MAKE_TYPED_HAS_TRAIT  (has_image_count,           int,  F.image_count()             )
MAKE_UNTYPED_HAS_TRAIT(has_int_call_operator,           F(0)                        )
MAKE_TYPED_HAS_TRAIT  (has_int_int_call_operator, int,  F(0)                        )
MAKE_TYPED_HAS_TRAIT  (has_set,                   void, F.set(0, I)                 )
MAKE_TYPED_HAS_TRAIT  (has_fill,                  void, F.fill(I)                   )
MAKE_TYPED_HAS_TRAIT  (has_move,                  IT,   F.move(0)                   )
MAKE_TYPED_HAS_TRAIT  (has_set_image_count,       void, F.set_image_count(0)        )

#define MAKE_BOOL_TRAIT(NAME, EXPR)\
	template<class FT> struct NAME : std::integral_constant<bool, EXPR>{};

MAKE_BOOL_TRAIT(is_id_func, 
	   has_preimage_count<FT>::value 
	&& has_int_call_operator<FT>::value
)
MAKE_BOOL_TRAIT(is_id_id_func, 
	   has_preimage_count<FT>::value 
	&& has_image_count<FT>::value && 
	has_int_int_call_operator<FT>::value
)
MAKE_BOOL_TRAIT(is_only_id_func, 
	   !is_id_id_func<FT>::value 
	&& is_id_func<FT>::value
)
MAKE_BOOL_TRAIT(is_mutable_id_func, 
	   is_id_func<FT>::value 
	&& has_set<FT>::value 
	&& has_fill<FT>::value 
	&& has_move<FT>::value
)
MAKE_BOOL_TRAIT(is_mutable_id_id_func, 
	   is_id_id_func<FT>::value 
	&& has_set<FT>::value 
	&& has_fill<FT>::value 
	&& has_move<FT>::value
	&& has_set_image_count<FT>::value
)

#undef F
#undef IT
#undef I
#undef MAKE_TYPED_HAS_TRAIT
#undef MAKE_UNTYPED_HAS_TRAIT
#undef MAKE_BOOL_TRAIT


#endif
