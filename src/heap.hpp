/******************************************************************************************[Heap.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#pragma once

#include "Vec.hpp"

namespace GanakInt {

// A heap implementation with support for decrease/increase key.
template<class Comp>
class Heap {
    Comp     lt;       // The heap is a minimum-heap with respect to this comparator
    vec<int> heap;     // Heap of integers
    vec<int> indices;  // Each integers position (index) in the Heap

    inline int left  (int i) { return i * 2 + 1; }
    inline int right (int i) { return (i + 1) * 2; }
    inline int parent(int i) { return (i - 1) >> 1; }
    void percolateUp(int i)
    {
        int x  = heap[i];
        int p  = parent(i);

        while (i != 0 && lt(x, heap[p])) {
            heap[i]          = heap[p];
            indices[heap[p]] = i;
            i                = p;
            p                = parent(p);
        }
        heap   [i] = x;
        indices[x] = i;
    }

    void percolateDown(int i)
    {
        int x = heap[i];
        while (left(i) < (int)heap.size()) {
            int child = right(i) < (int)heap.size() && lt(heap[right(i)], heap[left(i)]) ? right(i) : left(i);
            if (!lt(heap[child], x)) {
                break;
            }
            heap[i]          = heap[child];
            indices[heap[i]] = i;
            i                = child;
        }
        heap   [i] = x;
        indices[x] = i;
    }

public:
    Heap(const Comp& c) : lt(c) { }
    uint32_t size() const { return heap.size(); }
    bool empty() const { return heap.size() == 0; }
    bool in_heap(int n) const {
        return n < (int)indices.size() && indices[n] >= 0;
    }
    void decrease(int n) {
        assert(in_heap(n));
        percolateUp  (indices[n]);
    }
    void increase(int n) {
        assert(in_heap(n));
        percolateDown(indices[n]);
    }

    // Safe variant of insert/decrease/increase:
    void update(int n) {
        if (!in_heap(n)) insert(n);
        else {
            percolateUp(indices[n]);
            percolateDown(indices[n]);
        }
    }

    void insert(int n) {
        indices.growTo(n + 1, -1);
        assert(!in_heap(n));

        indices[n] = heap.size();
        heap.push_back(n);
        percolateUp(indices[n]);
    }

    int removeMin() {
        int x            = heap[0];
        heap[0]          = heap.back();
        indices[heap[0]] = 0;
        indices[x]       = -1;
        heap.pop_back();
        if (heap.size() > 1) percolateDown(0);
        return x;
    }

    void clear(bool dealloc = false) {
        for (auto i : heap) indices[i] = -1;
        heap.clear(dealloc);
    }

    size_t mem_used() const {
        size_t mem = 0;
        mem += heap.capacity()*sizeof(uint32_t);
        mem += indices.capacity()*sizeof(uint32_t);
        return mem;
    }

    bool heap_property(uint32_t i) const {
        return i >= heap.size()
            || ( (i == 0 || !lt(heap[i], heap[parent(i)]))
                  && heap_property( left(i)  )
                  && heap_property( right(i) )
            );
    }

    bool heap_property() const {
        return heap_property(0);
    }
};

}
