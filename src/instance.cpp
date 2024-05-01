/******************************************
Copyright (C) 2023 Authors of GANAK, see AUTHORS file

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

#include "instance.hpp"
#include "structures.hpp"
#include "clauseallocator.hpp"

#include <algorithm>
#include <fstream>
#include <limits>
#include <sys/stat.h>


Instance::Instance(const CounterConfiguration& _conf) : conf(_conf), stats (this, conf) {
  alloc = new ClauseAllocator(_conf);
  lbd_cutoff = conf.base_lbd_cutoff;
}

Instance::~Instance() {
  delete alloc;
}

void Instance::checkWatchLists() const {
  auto red_cls2 = longRedCls;
  // check for duplicates
  std::sort(red_cls2.begin(), red_cls2.end());
  for(uint32_t i = 1; i < red_cls2.size(); i++) {
    assert(red_cls2[i-1] != red_cls2[i]);
  }

  for(const auto& offs: longRedCls) {
    const auto& cl = *alloc->ptr(offs);
    if (!findOfsInWatch(watches[cl[0]].watch_list_, offs)) {
      cout << "ERROR: Did not find watch cl[0]!!" << endl;
      assert(false);
      exit(-1);
    }
    if (!findOfsInWatch(watches[cl[1]].watch_list_, offs)) {
      cout << "ERROR: Did not find watch cl[1]!!" << endl;
      assert(false);
      exit(-1);
    }
  }
}

bool Instance::findOfsInWatch(const vector<ClOffsBlckL>& ws, ClauseOfs off)  const
{
  for (auto& w: ws) if (w.ofs == off) { return true; }
  return false;
}

struct ClSorter {
  ClSorter(ClauseAllocator* _alloc, uint32_t _lbd_cutoff) :
    alloc(_alloc), lbd_cutoff(_lbd_cutoff) {}

  bool operator()(ClauseOfs& a, ClauseOfs& b) const {
    const auto& ah = *alloc->ptr(a);
    const auto& bh = *alloc->ptr(b);
    assert(ah.red);
    assert(bh.red);
    if (ah.lbd <= lbd_cutoff || bh.lbd <= lbd_cutoff) return ah.lbd < bh.lbd;
    if (ah.used != bh.used) return ah.used > bh.used;
    return ah.total_used > bh.total_used;
  }
  ClauseAllocator* alloc;
  const uint32_t lbd_cutoff;
};

void Instance::reduce_db() {
  stats.reduce_db++;
  if (stats.conflicts > (100ULL*1000ULL) && lbd_cutoff == conf.base_lbd_cutoff
      && num_low_lbd_cls < 100) {
    verb_print(1, " [rdb] bumping rdb cutoff to 3");
    lbd_cutoff++;
  }
  const auto cls_before = longRedCls.size();

  vector<ClauseOfs> tmp_red_cls = longRedCls;
  longRedCls.clear();
  num_low_lbd_cls = 0;
  num_used_cls = 0;
  uint32_t cannot_be_del = 0;
  sort(tmp_red_cls.begin(), tmp_red_cls.end(), ClSorter(alloc, lbd_cutoff));
  uint32_t cutoff = conf.rdb_cls_target;

  for(uint32_t i = 0; i < tmp_red_cls.size(); i++){
    const ClauseOfs& off = tmp_red_cls[i];
    auto& h = *alloc->ptr(off);
    if (h.lbd <= lbd_cutoff) num_low_lbd_cls++;
    else if (h.used) num_used_cls++;

    bool can_be_del = red_cl_can_be_deleted(off);
    cannot_be_del += !can_be_del;
    if (can_be_del && h.lbd > lbd_cutoff && (!conf.rdb_keep_used || !h.used) &&
        i > cutoff + num_low_lbd_cls + (conf.rdb_keep_used ? num_used_cls : 0)) {
      markClauseDeleted(off);
      stats.cls_deleted_since_compaction++;
      stats.cls_removed++;
    } else {
      longRedCls.push_back(off);
      h.used = 0;
    }
  }
  verb_print(1, "[rdb] cls before: " << cls_before << " after: " << longRedCls.size()
      << " low lbd: " << num_low_lbd_cls
      << " lbd cutoff: " << lbd_cutoff
      << " cutoff computed: " << cutoff
      << " cannot be del : " << cannot_be_del
      << " used: " << num_used_cls << " rdb: " << stats.reduce_db);
}

bool Instance::red_cl_can_be_deleted(ClauseOfs off){
  // only first literal may possibly have cl_ofs as antecedent
  Clause& cl = *alloc->ptr(off);
  if (isAntecedentOf(off, cl[0])) return false;
  return true;
}

void Instance::markClauseDeleted(const ClauseOfs off){
  Clause& cl = *alloc->ptr(off);
  watches[cl[0]].del_c(off);
  watches[cl[1]].del_c(off);
  alloc->clause_free(off);
}

void Instance::new_vars(const uint32_t n) {
  if (num_vars_set) {
    cout << "ERROR: you can only call new_vars() once!" << endl;
    exit(-1);
  }

  assert(variables_.empty());
  assert(values.empty());
  assert(watches.empty());
  assert(unit_clauses_.empty());
  assert(longRedCls.empty());

  variables_.resize(n + 1);
  values.resize(n + 1, X_TRI);
  watches.resize(n + 1);
  lbdHelper.resize(n+1, 0);
  num_vars_set = true;
}

Clause* Instance::addClause(const vector<Lit> &lits, bool red) {
  if (lits.size() == 1) {
    assert(!existsUnitClauseOf(lits[0].neg()) && "UNSAT is not dealt with");
    if (!existsUnitClauseOf(lits[0])) unit_clauses_.push_back(lits[0]);
    return nullptr;
  }

  if (lits.size() == 2) {
    add_bin_cl(lits[0], lits[1], red);
    return nullptr;
  }

  Clause* cl = alloc->new_cl(red, lits.size());
  for(uint32_t i = 0; i < lits.size(); i ++) (*cl)[i] = lits[i];
  attach_cl(alloc->get_offset(cl), lits);
  return cl;
}
