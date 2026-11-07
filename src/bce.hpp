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

#pragma once

// Dynamic blocked clause elimination (Lagniez, Marquis, Biere, SAT'24).
//
// Soundness in Ganak: a clause may only be blocked on a var >= opt_indep_support_end
// (non-indep, never branched on). Blocking on an optional-indep var would make it
// undefined by the indep support, and Ganak sums both branches of those vars.
//
// Undo is LIFO in "groups" (one propagate() call each). Ganak's chrono backtracking
// can unset a var processed in group g while vars of later groups stay assigned, so
// the whole suffix from g is popped and the still-assigned vars are re-processed.
// Because pops are by group, a triple kept in the list of the clause being processed
// (for any reason) is restored together with that clause, as in the paper. The only
// exception is re-examining a triple whose watcher is still inactive (recheck), and
// there the result is pinned by an event of the current group.

#include <cstdint>
#include <vector>
#include "lit.hpp"
#include "containers.hpp"
#include "structures.hpp"

namespace GanakInt {

struct BCEStats {
  uint64_t blocked = 0;
  uint64_t triples = 0;
  uint64_t cand_entries = 0;
  uint64_t rel_cls = 0;
  uint64_t skipped_lits = 0;
  uint64_t pops = 0;
  uint64_t reproc = 0;
  uint64_t moves = 0;
  double init_time = 0;
};

class BCE {
public:
  static constexpr uint32_t NONE = 0xffffffffU;
  static constexpr uint8_t INACT = 1;
  static constexpr uint8_t BLOCKED = 2;

  // cls[id] for id in [1, cls.size()), cls[0] unused. Ids >= first_bin are binary
  // clauses, and bin_pair[id-first_bin] is the analyzer's var-pair id of that clause.
  void init(const LiteralIndexedVector<TriValue>* vals, uint32_t nvars, uint32_t blk_var_start,
      const std::vector<std::vector<Lit>>& cls, uint32_t first_bin,
      const std::vector<uint32_t>& bin_pair, uint32_t num_pairs,
      uint64_t occ_prod_lim, uint64_t cand_lim);

  const uint8_t* state() const { return st.data(); }
  const uint32_t* pair_blocked_data() const { return pair_blocked.data(); }
  bool is_blocked(uint32_t id) const { return st[id] & BLOCKED; }

  void note_unset(uint32_t v) {
    if (proc_grp[v] != NONE && proc_grp[v] < min_grp) min_grp = proc_grp[v];
  }
  // After the trail has been truncated/compacted. Kept lits were moved to keep_from.
  void backtrack(uint32_t keep_from);
  void propagate(const std::vector<Lit>& trail);

  const BCEStats& get_stats() const { return stats; }
  bool useless() const { return tri.empty(); }
  uint32_t num_blocked_now() const { return cur_blocked; }

private:
  enum EvKind : uint8_t { EV_PROC, EV_SAT, EV_BLOCK, EV_RECHECK };
  struct Ev {
    EvKind kind;
    uint32_t id;   // var for EV_PROC, clause for SAT/BLOCK, triple for RECHECK
    uint32_t t;    // triple for EV_BLOCK
  };
  struct Tri {
    Lit l;            // blocking literal, in clause cl
    uint32_t cl;
    uint32_t cs, ce;  // candidate partner clauses: cands[cs..ce)
    uint32_t watcher; // NONE if no candidates at all
    uint32_t pos;     // index of watcher in cands
  };

  void proc_lit(Lit l, uint32_t g);
  void process_wl(uint32_t c);
  void recheck(uint32_t t);
  void block(uint32_t cl, uint32_t t);
  uint32_t find_active(Tri& tr, uint32_t& pos) const;
  bool is_bin(uint32_t id) const { return id >= first_bin; }

  const LiteralIndexedVector<TriValue>* vals = nullptr;
  uint32_t first_bin = 0;
  std::vector<uint32_t> cl_off;    // lits of clause id: lits[cl_off[id] .. cl_off[id+1])
  std::vector<Lit> lits;
  std::vector<uint32_t> occ_off;   // by lit raw
  std::vector<uint32_t> occ;
  std::vector<uint8_t> st;
  std::vector<uint32_t> bin_pair;
  std::vector<uint32_t> pair_blocked;

  std::vector<Tri> tri;
  std::vector<uint32_t> cands;
  std::vector<std::vector<uint32_t>> wl;

  std::vector<uint32_t> proc_grp;
  std::vector<uint8_t> var_rel;
  std::vector<Ev> evs;
  std::vector<uint32_t> grp_start;
  uint32_t min_grp = NONE;
  uint32_t trail_pos = 0;

  std::vector<Lit> reproc; // popped but still assigned: must be processed again
  std::vector<uint32_t> q;
  std::vector<uint32_t> to_recheck;
  std::vector<uint32_t> tmp_recheck;
  uint32_t cur_blocked = 0;
  BCEStats stats;
};

}
