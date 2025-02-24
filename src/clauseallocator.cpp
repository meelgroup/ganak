/******************************************
Copyright (C) 2009-2020 Mate Soos

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

#include "clauseallocator.hpp"

#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <cmath>
#include <gmpxx.h>
#include "counter.hpp"
#include "counter_config.hpp"
#include "time_mem.hpp"
#include "common.hpp"

using std::cout;
using std::endl;

using namespace GanakInt;

constexpr double ALLOC_GROW_MULT = 1.5;

template<typename T>
void* ClauseAllocator<T>::alloc_enough(uint32_t num_lits) {
  //Try to quickly find a place at the end of a data_start
  uint64_t neededbytes = sizeof(Clause) + sizeof(Lit)*num_lits;
  uint64_t needed = neededbytes/sizeof(uint32_t) + (bool)(neededbytes % sizeof(uint32_t));

  if (size + needed > capacity) {
    //Grow by default, but don't go under or over the limits
    uint64_t newcapacity = capacity * ALLOC_GROW_MULT;
    newcapacity = std::max<size_t>(newcapacity, MIN_LIST_SIZE);
    while (newcapacity < size+needed) {
        newcapacity *= ALLOC_GROW_MULT;
    }
    assert(newcapacity >= size+needed);
    newcapacity = std::min<size_t>(newcapacity, MAXSIZE);

    //Oops, not enough space anyway
    if (newcapacity < size + needed) {
      std::cerr << "ERROR: memory manager can't handle the load."
      << " size: " << size << " needed: " << needed << " newcapacity: " << newcapacity
      << endl;
      std::cout << "ERROR: memory manager can't handle the load."
      << " size: " << size << " needed: " << needed << " newcapacity: " << newcapacity
      << endl;
      exit(-1);
    }

    uint32_t* new_data_start = (uint32_t*)realloc(data_start , newcapacity*sizeof(uint32_t));
    if (new_data_start == nullptr) {
      std::cerr << "ERROR: while reallocating clause space" << endl;
      exit(-1);
    }
    data_start = new_data_start;
    capacity = newcapacity;
  }

  //Add clause to the set
  Clause* pointer = (Clause*)(data_start + size);
  size += needed;
  currently_used_sz += needed;
  return pointer;
}

template<typename T>
ClauseOfs ClauseAllocator<T>::get_offset(const Clause* ptr) const {
  return ((uint32_t*)ptr - data_start);
}

/**
@brief Frees a clause

If clause was binary, it frees it in quite a normal way. If it isn't, then it
needs to set the data in the Clause that it has been freed, and updates the
stack it belongs to such that the stack can now that its effectively used size
is smaller

NOTE: The size of claues can change. Therefore, currently_used_size can in fact
be incorrect, since it was incremented by the ORIGINAL size of the clause, but
when the clause is "freed", it is decremented by the POTENTIALLY SMALLER size
of the clause. Therefore, the "currently_used_size" is an overestimation!!
*/
template<typename T>
void ClauseAllocator<T>::clause_free(Clause* cl)
{
    assert(!cl->freed);
    cl->freed = 1;
    uint64_t est_num_cl = cl->sz;
    uint64_t bytes_freed = sizeof(Clause) + est_num_cl*sizeof(Lit);
    uint64_t elems_freed = bytes_freed/sizeof(uint32_t) + (bool)(bytes_freed % sizeof(uint32_t));
    currently_used_sz -= elems_freed;
}

template<typename T>
void ClauseAllocator<T>::clause_free(ClauseOfs offset)
{
  Clause* cl = ptr(offset);
  clause_free(cl);
}

template<typename T>
ClauseOfs ClauseAllocator<T>::move_cl(
    ClauseOfs* new_data_start
    , ClauseOfs*& new_ptr
    , Clause* old
) {
  uint64_t bytes_needed = sizeof(Clause) + old->sz*sizeof(Lit);
  uint64_t size_needed = bytes_needed/sizeof(uint32_t) + (bool)(bytes_needed % sizeof(uint32_t));
  memcpy(new_ptr, old, size_needed*sizeof(uint32_t));

  ClauseOfs new_offset = new_ptr-new_data_start;
  (*old)[0] = Lit::toLit(new_offset & 0xFFFFFFFF);
  old->reloced = true;
  new_sz_while_moving += size_needed;

  new_ptr += size_needed;
  return new_offset;
}

template<typename T>
void ClauseAllocator<T>::move_one_watchlist(
    vec<ClOffsBlckL>& ws, ClauseOfs* new_data_start, ClauseOfs*& new_ptr)
{
  for(auto& w: ws) {
    Clause* old = ptr(w.ofs);
    assert(!old->freed);
    Lit blocked = w.blckLit;
    if (old->reloced) {
      ClauseOfs new_offset = (*old)[0].raw();
      w = ClOffsBlckL(new_offset, blocked);
    } else {
      ClauseOfs new_offset = move_cl(new_data_start, new_ptr, old);
      w = ClOffsBlckL(new_offset, blocked);
    }
  }
}

/**
@brief If needed, compacts stacks, removing unused clauses

Firstly, the algorithm determines if the number of useless slots is large or
small compared to the problem size. If it is small, it does nothing. If it is
large, then it allocates new stacks, copies the non-freed clauses to these new
stacks, updates all pointers and offsets, and frees the original stacks.
*/
template<typename T>
bool ClauseAllocator<T>::consolidate(Counter<T>* solver , const bool force) {
  //If re-allocation is not really neccessary, don't do it
  //Neccesities:
  //1) There is too much memory allocated. Re-allocation will save space
  //2) There is too much empty, unused space (>30%)
  if (!force
      && (float_div(currently_used_sz, size) > 0.8 || capacity < (100ULL*1000ULL))
  ) {
    verb_print(1, "[mem] Not consolidating memory. Used sz/sz: " <<
        float_div(currently_used_sz, size)
        << " Currently used size: " << currently_used_sz/1000 << " K");
    return false;
  }
  const double my_time = cpu_time();
  new_sz_while_moving = 0;

  //Pointers that will be moved along
  uint32_t * const new_data_start = (uint32_t*)malloc(currently_used_sz*sizeof(uint32_t));
  uint32_t * new_ptr = new_data_start;

  assert(sizeof(uint32_t) % sizeof(Lit) == 0);

  for(auto& ws: solver->watches) move_one_watchlist(ws.watch_list_, new_data_start, new_ptr);
  update_offsets(solver->long_irred_cls, new_data_start, new_ptr);
  update_offsets(solver->long_red_cls, new_data_start, new_ptr);

  //Fix up var_data
  for (auto& vdata: solver->var_data) {
    if (vdata.ante.isAnt() && vdata.ante.isAClause()) {
      Clause* old = ptr(vdata.ante.as_cl());
      assert(!old->freed);
      ClauseOfs new_offset = (*old)[0].raw();
      vdata.ante = Antecedent(new_offset);
    }
  }

  //Update sizes
  const uint64_t old_size = size;
  size = new_ptr-new_data_start;
  capacity = currently_used_sz;
  currently_used_sz = new_sz_while_moving;
  free(data_start);
  data_start = new_data_start;

  const double time_used = cpu_time() - my_time;
  if (conf.verb) {
    size_t log_2_size = 0;
    if (size > 0) log_2_size = std::log2(size);
    verb_print(1, "[mem] consolidate "
      << " old-sz: " << print_value_kilo_mega(old_size*sizeof(uint32_t))
      << " new-sz: " << print_value_kilo_mega(size*sizeof(uint32_t))
      << " new bits offs: " << std::fixed << std::setprecision(2) << log_2_size
      << " T: " << time_used);
  }
  return true;
}

template<typename T>
void ClauseAllocator<T>::update_offsets(
    vector<ClauseOfs>& offsets,
    ClauseOfs* new_data_start,
    ClauseOfs*& new_ptr
) {
  for(ClauseOfs& offs: offsets) {
    Clause* old = ptr(offs);
    if (!old->reloced) offs = move_cl(new_data_start, new_ptr, old);
    else offs = (*old)[0].raw();
  }
}

template<typename T>
size_t ClauseAllocator<T>::mem_used() const
{
  uint64_t mem = 0;
  mem += capacity*sizeof(uint32_t);

  return mem;
}

template class GanakInt::ClauseAllocator<complex<mpq_class>>;
template class GanakInt::ClauseAllocator<mpz_class>;
template class GanakInt::ClauseAllocator<mpq_class>;
