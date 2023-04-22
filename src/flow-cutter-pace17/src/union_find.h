#ifndef UNION_FIND_H
#define UNION_FIND_H

#include "array_id_func.h"
#include <cassert>

//! An id-id-function that maps a node onto its components representative
struct UnionFind{
public:
	UnionFind():node_count_(0){}

	explicit UnionFind(int node_count):parent_(node_count), node_count_(node_count), component_count_(node_count){
		parent_.fill(-1);
	}

	void reset(){
		parent_.fill(-1);
		component_count_ = node_count_;
	}

	int preimage_count()const{return node_count_;}
	int image_count()const{return node_count_;}

	void unite(int l, int r){
		assert(0 <= l && l < node_count_);
		assert(0 <= r && r < node_count_);

		l = operator()(l);
		r = operator()(r);
		if(l != r){
			--component_count_;
			if(-parent_[l] < -parent_[r]){
				parent_[r] += parent_[l];
				parent_[l] = r;
			}else{
				parent_[l] += parent_[r];
				parent_[r] = l;
			}
		}
	}

	int operator()(int x)const{
		assert(0 <= x && x < node_count_);

		if(is_representative(x))
			return x;

		int y = x;
		while(!is_representative(y))
			y = parent_[y];

		int z = x;
		while(!is_representative(z)){
			int tmp = parent_[z];
			parent_[z] = y;
			z = tmp;
		}

		return y;
	}

	bool in_same(int x, int y)const{
		return (*this)(x) == (*this)(y);
	}

	bool is_representative(int v)const{
		assert(0 <= v && v < node_count_);
		return parent_(v) < 0;
	}

	int component_size(int v)const{
		assert(0 <= v && v < node_count_);
		if(is_representative(v))
			return -parent_(v);
		else
			return 0;		
	}

	int component_count()const{
		return component_count_;
	}

private:
	mutable ArrayIDFunc<int>parent_;
	int node_count_;
	int component_count_;
};

#endif
