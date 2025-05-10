#ifndef KWAY_HEAP_H
#define KWAY_HEAP_H

#include <utility>
#include <vector>
#include <cassert>
#include <functional>
#include <algorithm>

//#define INTERNAL_HEAP_CHECKS

template<class keyT, int k, class key_orderT = std::less<keyT>>
class kway_min_id_heap{
public:
	typedef keyT key_type;
	typedef key_orderT key_order_type;

	explicit kway_min_id_heap(int id_count, key_order_type _order = key_order_type()):
		heap_end(0), heap(id_count), id_pos(id_count, -1), order(std::move(_order)),
		contained_flags(id_count, false){
		check_id_invariants();
		check_order_invariants();
	}

	explicit kway_min_id_heap(key_order_type _order = key_order_type()):
		heap_end(0), order(std::move(_order)){
		check_id_invariants();
		check_order_invariants();
	}

	//! Takes a range of std::pair<int, key_type> elements (or any other struct with first and second) 
	//! and builds a heap containing these elements.
	template<class Range>
	void fill(const Range&range){
		clear();
		
		for(auto r:range){
			assert(0 <= r.first && r.first < (int)id_pos.size() && "range must contain valid id");

			heap[heap_end].id = r.first;
			heap[heap_end].key = r.second;

			id_pos[r.first] = heap_end;
			contained_flags[r.first] = true;
			++heap_end;
		}

		rebuild_heap();

		check_id_invariants();
		check_order_invariants();
	}


	void clear(){
		//heap_end = 0;
		//std::fill(contained_flags.begin(), contained_flags.end(), false);

		for(int i=0; i<heap_end; ++i)
			contained_flags[heap[i].id] = false;
		heap_end = 0;

		#ifdef INTERNAL_HEAP_CHECKS
		for(auto x:contained_flags)
			assert(!x);
		#endif

		check_id_invariants();
		check_order_invariants();
	}

	void reorder(key_order_type new_order){
		check_id_invariants();
		check_order_invariants();

		order = std::move(new_order);
		rebuild_heap();

		check_id_invariants();
		check_order_invariants();
	}

	void reset(int new_id_count = 0, [[maybe_unused]] key_order_type new_order = key_order_type()){
		heap.resize(new_id_count);
		id_pos.resize(new_id_count);
		contained_flags.resize(new_id_count);
		std::fill(contained_flags.begin(), contained_flags.end(), false);
		
		clear();

		check_id_invariants();
		check_order_invariants();
	}

	void reset(key_order_type new_order){
		clear();
		order = std::move(new_order);
	}

	bool empty()const{
		check_id_invariants();
		check_order_invariants();

		return heap_end == 0;
	}

	int size()const{
		return heap_end;
	}

	bool contains(int id)const{
		assert(0 <= id && id < (int)id_pos.size() && "id is in range");

		check_id_invariants();
		check_order_invariants();

		return contained_flags[id];
	}

	const key_type&get_key(int id)const{
		assert(0 <= id && id < (int)id_pos.size() && "id is in range");
		assert(contains(id) && "id is contained");

		check_id_invariants();
		check_order_invariants();

		return heap[id_pos[id]].key;
	}

	void push(int id, key_type key){
		assert(0 <= id && id < (int)id_pos.size() && "id is in range");
		assert(!contains(id) && "can not push an already existing id");
		
		contained_flags[id] = true;

		int new_pos = heap_end;
		++heap_end;
		heap[new_pos].id = id;
		heap[new_pos].key = std::move(key);
		id_pos[id] = new_pos;
		move_up(new_pos);

		check_id_invariants();
		check_order_invariants();
	}

	bool push_or_decrease_key(int id, key_type key){
		assert(0 <= id && id < (int)id_pos.size() && "id is in range");

		check_id_invariants();
		check_order_invariants();
		
		if(!contains(id)){
			push(id, key);
			return true;
		}else{
			if(order(key, heap[id_pos[id]].key)){
				heap[id_pos[id]].key = std::move(key);
				move_up(id_pos[id]);

				check_id_invariants();
				check_order_invariants();
		
				return true;
			}
		}
		return false;
	}

	bool push_or_increase_key(int id, key_type key){
		assert(0 <= id && id < (int)id_pos.size() && "id is in range");
		
		check_id_invariants();
		check_order_invariants();
		
		if(!contains(id)){
			push(id, key);
			return true;
		}else{
			if(order(heap[id_pos[id]].key, key)){
				heap[id_pos[id]].key = std::move(key);
				move_down(id_pos[id]);
				
				check_id_invariants();
				check_order_invariants();
		
				return true;
			}
		}
		return false;
	}

	bool push_or_set_key(int id, key_type key){
		assert(0 <= id && id < (int)id_pos.size() && "id is in range");
		
		check_id_invariants();
		check_order_invariants();
		
		if(!contains(id)){
			push(id, key);
			return true;
		}else{
			if(order(heap[id_pos[id]].key, key)){
				heap[id_pos[id]].key = std::move(key);
				move_down(id_pos[id]);
				
				check_id_invariants();
				check_order_invariants();
		
				return true;
			}else if(order(key, heap[id_pos[id]].key)){
				heap[id_pos[id]].key = std::move(key);
				move_up(id_pos[id]);

				check_id_invariants();
				check_order_invariants();
		
				return true;
			}
		}
		return false;
	}

	key_type peek_min_key()const{
		assert(!empty() && "heap is not empty");
		
		check_id_invariants();
		check_order_invariants();
		
		return heap[0].key;
	}

	int peek_min_id()const{
		assert(!empty() && "heap is not empty");

		check_id_invariants();
		check_order_invariants();
		
		return heap[0].id;
	}

	int pop(){
		assert(!empty() && "heap is not empty");

		check_id_invariants();
		check_order_invariants();
		
		if(heap_end == 1){
			heap_end = 0;
			id_pos[heap[0].id] = -1;
			contained_flags[heap[0].id] = false;
			
			check_id_invariants();
			check_order_invariants();
			
			return heap[0].id;
		}else{	
			int ret = heap[0].id;
			--heap_end;
			#ifdef INTERNAL_HEAP_CHECKS
			assert(ret != heap[heap_end].id);
			#endif
			heap[0].id = heap[heap_end].id;
			heap[0].key = std::move(heap[heap_end].key);
			id_pos[heap[0].id] = 0;
			id_pos[ret] = -1;
			move_down(0);
			contained_flags[ret] = false;
			
			check_id_invariants();
			check_order_invariants();

			return ret;
		}
	}

	int get_element_id(int pos)const{
		assert(0 <= pos && pos < heap_end && "element pos is in range");
		return heap[pos].id;
	}

	int get_element_key(int pos)const{
		assert(0 <= pos && pos < heap_end && "element pos is in range");
		return heap[pos].key;
	}

private:

	static int parent(int pos){
		assert(pos != 0);
		return (pos-1) / k;
	}

	static int children_begin(int pos){
		return k*pos+1;
	}

	static int children_end(int pos){
		return k*(pos+1)+1;
	}

	void rebuild_heap(){
		for(int i=heap_end-1; i>=0; --i)
			move_down(i);
	}

	void move_up(int pos){
		if(pos != 0){
			key_type key = std::move(heap[pos].key);
			int id = heap[pos].id;

			int parent_pos = parent(pos);
			while(order(key, heap[parent_pos].key)){
				heap[pos].id = heap[parent_pos].id;
				heap[pos].key = std::move(heap[parent_pos].key);
				id_pos[heap[parent_pos].id] = pos;

				pos = parent_pos;
				if(pos == 0)
					break;
				parent_pos = parent(pos);
			}

			heap[pos].id = id;
			heap[pos].key = std::move(key);
			id_pos[id] = pos;
		}
	}

	void move_down(int pos){
		key_type key = std::move(heap[pos].key);
		int id = heap[pos].id;

		for(;;){
			int begin = std::min(heap_end, children_begin(pos));
			int end = std::min(heap_end, children_end(pos));

			if(begin == end)
				break;

			int min_child_pos = begin;
			for(int i=begin+1; i<end; ++i){
				if(order(heap[i].key, heap[min_child_pos].key))
					min_child_pos = i;
			}

			if(!order(heap[min_child_pos].key, key))
				break;

			heap[pos].id = heap[min_child_pos].id;
			heap[pos].key = std::move(heap[min_child_pos].key);
			id_pos[heap[min_child_pos].id] = pos;

			pos = min_child_pos;
		}
		heap[pos].id = id;
		heap[pos].key = std::move(key);
		id_pos[id] = pos;
	}

	struct id_key_pair{
		int id;
		key_type key;
	};

	int heap_end;
	std::vector<id_key_pair>heap;
	std::vector<int>id_pos;
	key_order_type order;
	std::vector<bool> contained_flags;

	void check_id_invariants()const{
		#ifdef INTERNAL_HEAP_CHECKS
		for(int i=0; i<heap_end; ++i){
			assert(heap[i].id != -1);
			assert(0 <= heap[i].id);
			assert(heap[i].id < (int)id_pos.size());
			assert(contained_flags[heap[i].id]);
			assert(id_pos[heap[i].id] == i);
		}

		for(int i=0; i<(int)id_pos.size(); ++i){
			if(contained_flags[i]){
				assert(0 <= id_pos[i]);
				assert(id_pos[i] < heap_end);
				assert(heap[id_pos[i]].id == i);
			}
		}
		#endif
	}

	void check_order_invariants()const{
		#ifdef INTERNAL_HEAP_CHECKS
		for(int i=1; i<heap_end; ++i)
			assert(!order(heap[i].key, heap[parent(i)].key));
		#endif
	}
};


template<class keyT, int k, class key_orderT = std::less<keyT>>
class kway_max_id_heap{
public:
	typedef keyT key_type;
	typedef key_orderT key_order_type;

	explicit kway_max_id_heap(int id_count, key_order_type order = key_order_type())
		:heap(id_count, inverted_order(std::move(order))){}

	explicit kway_max_id_heap(key_order_type order = key_order_type())
		:heap(inverted_order(std::move(order))){}

	//! Removes all elements from the heap.
	void clear(){
		heap.clear();
	}

	//! Initializes the heap with a range over pair<int, keyT> elements.
	//! This function has linear running time, whereas repeated push operations have a total running time in O(n log n)
	template<class Range>
	void fill(const Range&range){
		heap.fill(range);
	}
	
	//! Rebuilds the queue with a different order
	void reorder(key_order_type new_order){
		heap.reorder(inverted_order(std::move(new_order)));
	}

	void reset(int new_id_count = 0, key_order_type new_order = key_order_type()){
		heap.reset(new_id_count, inverted_order(std::move(new_order)));
	}

	void reset(key_order_type new_order){
		heap.reset(std::move(new_order));
	}

	bool empty()const{
		return heap.empty();
	}

	int size()const{
		return heap.size();
	}

	bool contains(int id)const{
		return heap.contains(id);
	}

	const key_type&get_key(int id)const{
		return heap.get_key(id);
	}

	void push(int id, key_type key){
		heap.push(id, std::move(key));
	}

	//! Returns true if the element moved its position within the heap.
	bool push_or_decrease_key(int id, key_type key){
		return heap.push_or_increase_key(id, std::move(key));
	}

	//! Returns true if the element moved its position within the heap.
	bool push_or_increase_key(int id, key_type key){
		return heap.push_or_decrease_key(id, std::move(key));
	}

	//! Returns true if the element moved its position within the heap.
	bool push_or_set_key(int id, key_type key){
		return heap.push_or_set_key(id, std::move(key));
	}

	key_type peek_max_key()const{
		return heap.peek_min_key();
	}

	int peek_max_id()const{
		return heap.peek_min_id();
	}

	int pop(){
		return heap.pop();
	}

private:
	struct inverted_order{
		inverted_order(){}
		inverted_order(key_order_type _order):order(std::move(_order)){}
		bool operator()(const key_type&l, const key_type&r){
			return order(r, l);
		}
		key_order_type order;
	};
	kway_min_id_heap<key_type, k, inverted_order>heap;
};

const int standard_heap_arity = 4;

template<class keyT, class key_orderT = std::less<keyT>>
class min_id_heap : public kway_min_id_heap<keyT, standard_heap_arity, key_orderT>{
private:
	typedef kway_min_id_heap<keyT, standard_heap_arity, key_orderT> super_type;
public:
	typedef typename super_type::key_order_type key_order_type;

	explicit min_id_heap(int id_count, key_order_type _order = key_order_type())
		:super_type(id_count, std::move(_order)){}

	explicit min_id_heap(key_order_type order = key_order_type())
		:super_type(std::move(order)){}
};

template<class keyT, class key_orderT = std::less<keyT>>
class max_id_heap : public kway_max_id_heap<keyT, standard_heap_arity, key_orderT>{
private:
	typedef kway_max_id_heap<keyT, standard_heap_arity, key_orderT> super_type;
public:
	typedef typename super_type::key_order_type key_order_type;

	explicit max_id_heap(int id_count, key_order_type order = key_order_type())
		:super_type(id_count, std::move(order)){}

	explicit max_id_heap(key_order_type order = key_order_type())
		:super_type(std::move(order)){}
};

#endif

