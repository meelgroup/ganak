/*
 * instance.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#include "instance.h"
#include "structures.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <sys/stat.h>

void Instance::compactConflictLiteralPool(){
  stats.cls_deleted_since_compaction = 0;

  // We must sort, or we can later overwrite offs1 with offs2 and then offs2 can
  // happen to be a new offset and then we overwrite that, etc... A mess.
  std::sort(red_cls.begin(), red_cls.end());

  auto write_pos = conflict_clauses_begin();
  vector<ClauseOfs> tmp_conflict_clauses = red_cls;
  red_cls.clear();
  for(const auto& offs: tmp_conflict_clauses){
    /* size_t origsize = minimize_cl_with_bins(offs); */
    Lit l1 = *beginOf(offs);
    Lit l2 = *beginOf(offs+1);
    auto read_pos = beginOf(offs) - ClHeader::overheadInLits();
    for(uint32_t i = 0; i < ClHeader::overheadInLits(); i++) *(write_pos++) = *(read_pos++);
    ClauseOfs new_ofs =  write_pos - lit_pool_.begin();
    red_cls.push_back(new_ofs);
    // first substitute antecedent if l1 implied something
    if(isAntecedentOf(offs, l1)) var(l1).ante = Antecedent(new_ofs);

    // now redo the watches
    litWatchList(l1).replaceWatchLinkTo(offs,new_ofs);
    litWatchList(l2).replaceWatchLinkTo(offs,new_ofs);

    // next, copy clause data
    /* size_t i = 0; */
    assert(read_pos == beginOf(offs));
    while(*read_pos != SENTINEL_LIT) {*(write_pos++) = *(read_pos++);}
    /* for(; i < origsize; i++) read_pos++; */
    *(write_pos++) = SENTINEL_LIT;
  }
  lit_pool_.erase(write_pos,lit_pool_.end());
  stats.compactions++;
  SLOW_DEBUG_DO(checkWatchLists());
}

void Instance::checkWatchLists() const {
  auto red_cls2 = red_cls;
  // check for duplicates
  std::sort(red_cls2.begin(), red_cls2.end());
  for(uint32_t i = 1; i < red_cls2.size(); i++) {
    assert(red_cls2[i-1] != red_cls2[i]);
  }

  for(const auto& offs: red_cls) {
    const auto& header = getHeaderOf(offs);
    assert(!header.marked_deleted);

    Lit l1 = *beginOf(offs);
    Lit l2 = *beginOf(offs+1);
    if (!findOfsInWatch(litWatchList(l1).watch_list_, offs)) {
      cout << "ERROR: Did not find watch l1!!" << endl;
      assert(false);
      exit(-1);
    }
    if (!findOfsInWatch(litWatchList(l2).watch_list_, offs)) {
      cout << "ERROR: Did not find watch l2!!" << endl;
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
  ClSorter(const vector<Lit>& lit_pool, uint32_t lbd_cutoff) :
    lit_pool_(lit_pool), lbd_cutoff_(lbd_cutoff) {}

  const ClHeader& getHeaderOf(ClauseOfs cl_ofs) const {
    return *reinterpret_cast<const ClHeader *>(
        &lit_pool_[cl_ofs - ClHeader::overheadInLits()]);
  }
  bool operator()(ClauseOfs& a, ClauseOfs& b) const {
    const auto& ah = getHeaderOf(a);
    const auto& bh = getHeaderOf(b);
    if (ah.lbd <= lbd_cutoff_ || bh.lbd <= lbd_cutoff_) return ah.lbd < bh.lbd;
    if (ah.used != bh.used) return ah.used > bh.used;
    return ah.total_used > bh.total_used;
  }
  const vector<Lit>& lit_pool_;
  const uint32_t lbd_cutoff_;
};

// TODO with propagation and long clauses, etc.
size_t Instance::minimize_cl_with_bins(ClauseOfs off) {
  SLOW_DEBUG_DO(for(const auto& s: tmp_seen) assert(s == 0););
  uint32_t rem = 0;
  tmp_minim_with_bins.clear();
  for (auto l = beginOf(off); *l != SENTINEL_LIT; l++) {
    tmp_seen[l->toPosInt()] = 1; tmp_minim_with_bins.push_back(*l);}
  for(const auto& l: tmp_minim_with_bins) {
    if (!tmp_seen[l.toPosInt()]) continue;
    const auto& w = watches_[l].binary_links_;
    for(const auto& l2: w) {
      assert(l.var() != l2.var());
      if (tmp_seen[(l2.neg()).toPosInt()]) { tmp_seen[(l2.neg()).toPosInt()] = 0; rem++; }
    }
  }
  uint32_t at = 0;
  for(uint32_t i = 0; i < tmp_minim_with_bins.size(); i++) {
    Lit l = tmp_minim_with_bins[i];
    if (tmp_seen[l.toPosInt()] || i <= 1) {
      *beginOf(off+at) = l;
      at++;
      tmp_seen[l.toPosInt()] = 0;
    }
  }
  while (at < tmp_minim_with_bins.size()) {*beginOf(off+at) = SENTINEL_LIT; at++;}
  stats.rem_lits_with_bins+=rem;
  stats.rem_lits_tried++;
  return tmp_minim_with_bins.size(); // original size
}

void Instance::reduceDB() {
  stats.reduceDBs++;
  if (stats.conflicts > (100ULL*1000ULL) && lbd_cutoff == 2
      && num_low_lbd_cls < 100) {
    verb_print(1, " [rdb] bumping rdb cutoff to 3");
    lbd_cutoff = 3;
  }
  const auto cls_before = red_cls.size();

  vector<ClauseOfs> tmp_red_cls = red_cls;
  red_cls.clear();
  num_low_lbd_cls = 0;
  num_used_cls = 0;
  uint32_t cannot_be_del = 0;
  sort(tmp_red_cls.begin(), tmp_red_cls.end(), ClSorter(lit_pool_, lbd_cutoff));
  uint32_t cutoff = config_.rdb_cls_target;

  for(uint32_t i = 0; i < tmp_red_cls.size(); i++){
    const ClauseOfs& off = tmp_red_cls[i];
    auto& h = getHeaderOf(off);
    if (h.lbd <= lbd_cutoff) num_low_lbd_cls++;
    else if (h.used) num_used_cls++;

    bool can_be_del = red_cl_can_be_deleted(off);
    cannot_be_del += !can_be_del;
    if (can_be_del && h.lbd > lbd_cutoff && (!config_.rdb_keep_used || !h.used) &&
        i > cutoff + num_low_lbd_cls + (config_.rdb_keep_used ? num_used_cls : 0)) {
      markClauseDeleted(off);
      stats.cls_deleted_since_compaction++;
      stats.cls_removed++;
    } else {
      red_cls.push_back(off);
      h.used = 0;
    }
  }
  verb_print(1, "[rdb] cls before: " << cls_before << " after: " << red_cls.size()
      << " low lbd: " << num_low_lbd_cls
      << " lbd cutoff: " << lbd_cutoff
      << " cutoff computed: " << cutoff
      << " cannot be del : " << cannot_be_del
      << " used: " << num_used_cls << " rdb: " << stats.reduceDBs);
}

bool Instance::red_cl_can_be_deleted(ClauseOfs cl_ofs){
  // only first literal may possibly have cl_ofs as antecedent
  if (isAntecedentOf(cl_ofs, *beginOf(cl_ofs))) return false;
  return true;
}

void Instance::markClauseDeleted(ClauseOfs cl_ofs){
  litWatchList(*beginOf(cl_ofs)).removeWatchLinkTo(cl_ofs);
  litWatchList(*(beginOf(cl_ofs)+1)).removeWatchLinkTo(cl_ofs);
}

void Instance::new_vars(const uint32_t n) {
  assert(variables_.empty());
  assert(lit_values_.empty());
  assert(occ_lists_.empty());
  assert(watches_.empty());
  assert(lit_pool_.empty());
  assert(unit_clauses_.empty());
  assert(red_cls.empty());

  lit_pool_.push_back(SENTINEL_LIT);
  variables_.resize(n + 1);
  lit_values_.resize(n + 1, X_TRI);
  occ_lists_.resize(n + 1);
  watches_.resize(n + 1);
  target_polar.resize(n + 1);
  lbdHelper.resize(n+1, 0);
}

void Instance::add_irred_cl(const vector<Lit>& lits) {
  if (lits.empty()) {
    cout << "ERROR: UNSAT should have been caught by external SAT solver" << endl;
    exit(-1);
  }
  for(const auto& l: lits) assert(l.var() <= nVars());
  stats.incorporateIrredClauseData(lits);
  ClauseOfs cl_ofs = addClause(lits, false);
  if (lits.size() >= 3)
    for (const auto& l : lits)
      occ_lists_[l].push_back(cl_ofs);
}
