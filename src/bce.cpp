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

#include "bce.hpp"
#include "time_mem.hpp"
#include <cassert>

using namespace GanakInt;

void BCE::init(const LiteralIndexedVector<TriValue>* _vals, uint32_t nvars, uint32_t blk_var_start,
    const std::vector<std::vector<Lit>>& cls, uint32_t _first_bin,
    const std::vector<uint32_t>& _bin_pair, uint32_t num_pairs,
    uint64_t occ_prod_lim, uint64_t cand_lim) {
  const double my_time = cpu_time();
  vals = _vals;
  first_bin = _first_bin;
  bin_pair = _bin_pair;
  pair_blocked.assign(num_pairs, 0);
  const uint32_t ncls = cls.size();
  const uint32_t nlits = 2*(nvars+1);

  cl_off.assign(ncls+1, 0);
  lits.clear();
  for (uint32_t id = 0; id < ncls; id++) {
    cl_off[id] = lits.size();
    for (const auto& l : cls[id]) lits.push_back(l);
  }
  cl_off[ncls] = lits.size();

  occ_off.assign(nlits+1, 0);
  for (uint32_t id = 1; id < ncls; id++) for (const auto& l : cls[id]) occ_off[l.raw()+1]++;
  for (uint32_t i = 0; i < nlits; i++) occ_off[i+1] += occ_off[i];
  occ.assign(occ_off[nlits], 0);
  {
    std::vector<uint32_t> at(occ_off.begin(), occ_off.end()-1);
    for (uint32_t id = 1; id < ncls; id++) for (const auto& l : cls[id]) occ[at[l.raw()]++] = id;
  }

  st.assign(ncls, 0);
  wl.assign(ncls, {});
  tri.clear();
  cands.clear();
  proc_grp.assign(nvars+1, NONE);
  evs.clear();
  grp_start.clear();
  min_grp = NONE;
  trail_pos = 0;
  q.clear();
  to_recheck.clear();
  reproc.clear();
  cur_blocked = 0;

  std::vector<uint32_t> mark(nlits, 0);
  uint32_t stamp = 0;
  for (uint32_t v = blk_var_start; v <= nvars; v++) {
    for (const bool sign : {false, true}) {
      const Lit l(v, sign);
      const uint32_t* p_beg = occ.data() + occ_off[l.raw()];
      const uint32_t* p_end = occ.data() + occ_off[l.raw()+1];
      const uint32_t* n_beg = occ.data() + occ_off[l.neg().raw()];
      const uint32_t* n_end = occ.data() + occ_off[l.neg().raw()+1];
      if (p_beg == p_end) continue;
      if ((uint64_t)(p_end-p_beg)*(uint64_t)(n_end-n_beg) > occ_prod_lim
          || cands.size() > cand_lim) {
        stats.skipped_lits++;
        continue;
      }
      for (const uint32_t* a = p_beg; a != p_end; a++) {
        stamp++;
        for (uint32_t i = cl_off[*a]; i < cl_off[*a+1]; i++) mark[lits[i].raw()] = stamp;
        Tri tr;
        tr.l = l;
        tr.cl = *a;
        tr.cs = cands.size();
        for (const uint32_t* b = n_beg; b != n_end; b++) {
          bool taut = false;
          for (uint32_t i = cl_off[*b]; i < cl_off[*b+1]; i++) {
            const Lit m = lits[i];
            if (m.var() != v && mark[m.neg().raw()] == stamp) { taut = true; break; }
          }
          if (!taut) cands.push_back(*b);
        }
        tr.ce = cands.size();
        const uint32_t t = tri.size();
        if (tr.cs == tr.ce) {
          tr.watcher = NONE;
          tr.pos = 0;
          to_recheck.push_back(t);
        } else {
          tr.watcher = cands[tr.cs];
          tr.pos = tr.cs;
          wl[tr.watcher].push_back(t);
        }
        tri.push_back(tr);
      }
    }
  }
  stats.triples = tri.size();
  stats.cand_entries = cands.size();
  stats.init_time = cpu_time() - my_time;
}

uint32_t BCE::find_active(Tri& tr, uint32_t& pos) const {
  if (tr.cs == tr.ce) return NONE;
  uint32_t i = (tr.pos >= tr.cs && tr.pos < tr.ce) ? tr.pos : tr.cs;
  for (uint32_t k = tr.cs; k < tr.ce; k++) {
    if (!(st[cands[i]] & INACT)) { pos = i; return cands[i]; }
    if (++i == tr.ce) i = tr.cs;
  }
  return NONE;
}

void BCE::block(uint32_t cl, uint32_t t) {
  assert(!(st[cl] & INACT));
  st[cl] |= INACT | BLOCKED;
  evs.push_back({EV_BLOCK, cl, t});
  if (is_bin(cl)) pair_blocked[bin_pair[cl-first_bin]]++;
  q.push_back(cl);
  cur_blocked++;
  stats.blocked++;
}

void BCE::process_wl(uint32_t c) {
  auto& w = wl[c];
  uint32_t j = 0;
  for (uint32_t i = 0; i < w.size(); i++) {
    const uint32_t t = w[i];
    Tri& tr = tri[t];
    if (tr.watcher != c) continue; // stale, it moved away during a recheck
    if (proc_grp[tr.l.var()] != NONE || (st[tr.cl] & INACT)) { w[j++] = t; continue; }
    uint32_t pos;
    const uint32_t b = find_active(tr, pos);
    if (b != NONE) {
      assert(b != c);
      tr.watcher = b;
      tr.pos = pos;
      wl[b].push_back(t);
      stats.moves++;
      continue;
    }
    block(tr.cl, t);
    w[j++] = t;
  }
  w.resize(j);
}

void BCE::recheck(uint32_t t) {
  Tri& tr = tri[t];
  if (tr.watcher != NONE && !(st[tr.watcher] & INACT)) return;
  if (proc_grp[tr.l.var()] != NONE || (st[tr.cl] & INACT)) {
    evs.push_back({EV_RECHECK, t, NONE});
    return;
  }
  uint32_t pos;
  const uint32_t b = find_active(tr, pos);
  if (b != NONE) {
    tr.watcher = b;
    tr.pos = pos;
    wl[b].push_back(t);
    return;
  }
  block(tr.cl, t);
}

void BCE::proc_lit(const Lit l, const uint32_t g) {
  const uint32_t v = l.var();
  if (proc_grp[v] != NONE) return;
  proc_grp[v] = g;
  evs.push_back({EV_PROC, v, NONE});
  for (uint32_t k = occ_off[l.raw()]; k < occ_off[l.raw()+1]; k++) {
    const uint32_t c = occ[k];
    if (st[c] & INACT) continue;
    st[c] |= INACT;
    evs.push_back({EV_SAT, c, NONE});
    q.push_back(c);
  }
}

void BCE::propagate(const std::vector<Lit>& trail) {
  if (trail_pos >= trail.size() && to_recheck.empty() && reproc.empty()) return;
  const uint32_t g = grp_start.size();
  grp_start.push_back(evs.size());

  for (const Lit l : reproc) if ((*vals)[l] == T_TRI) proc_lit(l, g);
  reproc.clear();
  for (uint32_t i = trail_pos; i < trail.size(); i++) proc_lit(trail[i], g);
  trail_pos = trail.size();

  tmp_recheck.swap(to_recheck);
  for (const auto& t : tmp_recheck) recheck(t);
  tmp_recheck.clear();

  for (uint32_t i = 0; i < q.size(); i++) process_wl(q[i]);
  q.clear();

  if (evs.size() == grp_start.back()) grp_start.pop_back();
}

void BCE::backtrack(uint32_t keep_from) {
  if (keep_from < trail_pos) trail_pos = keep_from;
  if (min_grp == NONE) return;
  assert(min_grp < grp_start.size());
  const uint32_t until = grp_start[min_grp];
  stats.pops++;
  while (evs.size() > until) {
    const Ev e = evs.back();
    evs.pop_back();
    switch (e.kind) {
      case EV_PROC:
        proc_grp[e.id] = NONE;
        // chrono BT kept it on the trail, possibly before trail_pos
        if ((*vals)[Lit(e.id, true)] != X_TRI) {
          reproc.push_back(Lit(e.id, (*vals)[Lit(e.id, true)] == T_TRI));
          stats.reproc++;
        }
        break;
      case EV_SAT:
        st[e.id] = 0;
        break;
      case EV_BLOCK:
        st[e.id] = 0;
        if (is_bin(e.id)) pair_blocked[bin_pair[e.id-first_bin]]--;
        cur_blocked--;
        to_recheck.push_back(e.t);
        break;
      case EV_RECHECK:
        to_recheck.push_back(e.id);
        break;
    }
  }
  grp_start.resize(min_grp);
  min_grp = NONE;
}
