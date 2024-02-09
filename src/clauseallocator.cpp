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

#include <stdlib.h>
#include <algorithm>
#include <string.h>
#include <limits>
#include <cassert>
#include <cmath>
#include "counter.hpp"
#include "counter_config.hpp"
#include "time_mem.hpp"
#include "common.hpp"

using std::cout;
using std::endl;

#define MIN_LIST_SIZE (50000 * (sizeof(Clause) + 4*sizeof(Lit))/sizeof(uint32_t))
#define ALLOC_GROW_MULT 1.5
#define MAXSIZE ((1ULL << 32)-1)

ClauseAllocator::ClauseAllocator(const CounterConfiguration& _conf) :
    dataStart(nullptr)
    , size(0)
    , capacity(0)
    , currentlyUsedSize(0)
    , conf(_conf)
{
    assert(MIN_LIST_SIZE < MAXSIZE);
}

ClauseAllocator::~ClauseAllocator() { free(dataStart); }

void* ClauseAllocator::allocEnough( uint32_t num_lits) {
  //Try to quickly find a place at the end of a dataStart
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

    //Reallocate data
    uint32_t* new_dataStart;
    new_dataStart = (uint32_t*)realloc(
      dataStart
      , newcapacity*sizeof(uint32_t)
    );

    //Realloc failed?
    if (new_dataStart == nullptr) {
      std::cerr << "ERROR: while reallocating clause space" << endl;
      exit(-1);
    }
    dataStart = new_dataStart;

    //Update capacity to reflect the update
    capacity = newcapacity;
  }

  //Add clause to the set
  Clause* pointer = (Clause*)(dataStart + size);
  size += needed;
  currentlyUsedSize += needed;
  return pointer;
}

ClauseOfs ClauseAllocator::get_offset(const Clause* ptr) const {
  return ((uint32_t*)ptr - dataStart);
}

/**
@brief Frees a clause

If clause was binary, it frees it in quite a normal way. If it isn't, then it
needs to set the data in the Clause that it has been freed, and updates the
stack it belongs to such that the stack can now that its effectively used size
is smaller

NOTE: The size of claues can change. Therefore, currentlyUsedSizes can in fact
be incorrect, since it was incremented by the ORIGINAL size of the clause, but
when the clause is "freed", it is decremented by the POTENTIALLY SMALLER size
of the clause. Therefore, the "currentlyUsedSizes" is an overestimation!!
*/
void ClauseAllocator::clauseFree(Clause* cl)
{
    assert(!cl->freed);
    cl->freed = 1;
    uint64_t est_num_cl = cl->sz;
    est_num_cl = std::max(est_num_cl, (uint64_t)3); //we sometimes allow gauss to allocate 3-long clauses
    uint64_t bytes_freed = sizeof(Clause) + est_num_cl*sizeof(Lit);
    uint64_t elems_freed = bytes_freed/sizeof(uint32_t) + (bool)(bytes_freed % sizeof(uint32_t));
    currentlyUsedSize -= elems_freed;
}

void ClauseAllocator::clauseFree(ClauseOfs offset)
{
  Clause* cl = ptr(offset);
  clauseFree(cl);
}

ClauseOfs ClauseAllocator::move_cl(
    ClauseOfs* newDataStart
    , ClauseOfs*& new_ptr
    , Clause* old
) const {
  uint64_t bytesNeeded = sizeof(Clause) + old->sz*sizeof(Lit);
  uint64_t sizeNeeded = bytesNeeded/sizeof(uint32_t) + (bool)(bytesNeeded % sizeof(uint32_t));
  memcpy(new_ptr, old, sizeNeeded*sizeof(uint32_t));

  ClauseOfs new_offset = new_ptr-newDataStart;
  (*old)[0] = Lit::toLit(new_offset & 0xFFFFFFFF);
  old->reloced = true;

  new_ptr += sizeNeeded;
  return new_offset;
}

void ClauseAllocator::move_one_watchlist(
    vector<ClOffsBlckL>& ws, ClauseOfs* newDataStart, ClauseOfs*& new_ptr)
{
  for(auto& w: ws) {
    Clause* old = ptr(w.ofs);
    assert(!old->freed);
    Lit blocked = w.blckLit;
    if (old->reloced) {
      ClauseOfs new_offset = (*old)[0].raw();
      w = ClOffsBlckL(new_offset, blocked);
    } else {
      ClauseOfs new_offset = move_cl(newDataStart, new_ptr, old);
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
bool ClauseAllocator::consolidate(Counter* solver , const bool force) {
  //If re-allocation is not really neccessary, don't do it
  //Neccesities:
  //1) There is too much memory allocated. Re-allocation will save space
  //2) There is too much empty, unused space (>30%)
  if (!force
      && (float_div(currentlyUsedSize, size) > 0.8 || currentlyUsedSize < (100ULL*1000ULL))
  ) {
    verb_print(1, "[mem] Not consolidating memory. Used sz/sz: " <<
        float_div(currentlyUsedSize, size)
        << " Currently used size: " << currentlyUsedSize/1000 << " K");
    return false;
  }
  const double myTime = cpuTime();

  //Pointers that will be moved along
  uint32_t * const newDataStart = (uint32_t*)malloc(currentlyUsedSize*sizeof(uint32_t));
  uint32_t * new_ptr = newDataStart;

  assert(sizeof(uint32_t) % sizeof(Lit) == 0);

  for(auto& ws: solver->watches_) move_one_watchlist(ws.watch_list_, newDataStart, new_ptr);
  update_offsets(solver->longIrredCls, newDataStart, new_ptr);
  update_offsets(solver->longRedCls, newDataStart, new_ptr);

  //Fix up variables_
  for (auto& vdata: solver->variables_) {
    if (vdata.ante.isAnt() && vdata.ante.isAClause()) {
      Clause* old = ptr(vdata.ante.asCl());
      assert(!old->freed);
      ClauseOfs new_offset = (*old)[0].raw();
      vdata.ante = Antecedent(new_offset);
    }
  }

  //Update sizes
  const uint64_t old_size = size;
  size = new_ptr-newDataStart;
  capacity = currentlyUsedSize;
  currentlyUsedSize = size;
  free(dataStart);
  dataStart = newDataStart;

  const double time_used = cpuTime() - myTime;
  if (conf.verb) {
    size_t log_2_size = 0;
    if (size > 0) log_2_size = std::log2(size);
    verb_print(1, "[mem] consolidate "
      << " old-sz: " << print_value_kilo_mega(old_size*sizeof(uint32_t))
      << " new-sz: " << print_value_kilo_mega(size*sizeof(uint32_t))
      << " new bits offs: " << std::fixed << std::setprecision(2) << log_2_size
      << " T: " << time_used);;
  }
  return true;
}

void ClauseAllocator::update_offsets(
    vector<ClauseOfs>& offsets,
    ClauseOfs* newDataStart,
    ClauseOfs*& new_ptr
) {
  for(ClauseOfs& offs: offsets) {
    Clause* old = ptr(offs);
    if (!old->reloced) offs = move_cl(newDataStart, new_ptr, old);
    else offs = (*old)[0].raw();
  }
}

size_t ClauseAllocator::mem_used() const
{
  uint64_t mem = 0;
  mem += capacity*sizeof(uint32_t);

  return mem;
}
