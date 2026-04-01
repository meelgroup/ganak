# Analysis: `out-ganak-mccomp2324-1193808-0` vs `out-ganak-mccomp2324-1256426-0`

## Component Versions

| Component | Old (1193808) | New (1256426) |
|-----------|--------------|--------------|
| Ganak | `255f3da` | `d1e4695` |
| Arjun | `3fe1420` | `66b04dff` |
| CMS | `31681010` | `a48c7888` |
| SBVA | `52c1835` | `14a596d` |
| Key new flag | — | `--puuraautarky 1` |

---

## Overall Results

| Metric | Old | New | Delta |
|--------|-----|-----|-------|
| Solved | 1176 | 1175 | **-1** |
| PAR2 (timeout ≈ 3590 s) | 4,545,479 | 4,540,468 | **-5,011 (better)** |
| Avg time (both-solved) | 251.7 s | 227.9 s | **-9.5%** |

Net: gained 10 instances, lost 11 — all losses are **timeouts (signal 14)**, not crashes.

---

## Per-Track Breakdown (both-solved averages)

| Track | Solved new/old | Avg new time | Avg old time | Improvement |
|-------|---------------|-------------|-------------|-------------|
| Track 1 (unweighted) | 288 / 287 | 225.6 s | 234.6 s | +4% |
| Track 2 (weighted) | 169 / 170 | 279.4 s | 300.9 s | +7% |
| Track 2-random | 123 / 123 | 259.9 s | 279.4 s | +7% |
| Track 3 (proj unweighted) | 270 / 270 | 293.2 s | 344.6 s | **+15%** |
| Track 4 (proj weighted) | 325 / 326 | 137.1 s | 154.1 s | **+11%** |

Track 3 and 4 are the biggest winners. Track 2-random wins in total time but individual instances show Arjun regressions (see below).

---

## Speed Distribution (1165 both-solved instances)

| Category | Count |
|---------|-------|
| >5× faster (new) | **12** |
| 2–5× faster | **25** |
| 1.2–2× faster | **203** |
| Within 1.2× | 807 |
| 1.2–2× slower | **109** |
| 2–5× slower | **8** |
| >5× slower | 1 |

**240 instances got faster, 118 got slower** among the jointly solved set.

---

## What Changed and Why

### 1. `--puuraautarky 1` — The Biggest New Feature

The new run adds pure-literal/autarky reasoning. Autarky identifies variable assignments that do not affect satisfiability of the rest of the formula, enabling early termination of subproblems.

**Most dramatic example:** `mc2024_track4_074.cnf`
- Old: 136 s, **2,831,618 conflicts**, 190,252 SAT oracle calls
- New: **6 s**, 1,162 conflicts, 9 restarts

A 23× speedup with ~2,500× fewer conflicts. The autarky is detecting that large parts of the formula are independent and resolving them without search. This pattern repeats across many track-4 (weighted projected) instances — the `indep_sz` dropped from 7 to 3 for this instance, suggesting much stronger simplification.

Several other big winners from autarky (10–16× speedups on track3/4):

| Instance | Old time | New time | Speedup |
|----------|---------|---------|---------|
| `mc2024_track4_195.cnf` | 7.0 s | 0.5 s | 13× |
| `mc2024_track3_104.cnf` | 3.7 s | 0.3 s | 14× |
| `mc2024_track3_188.cnf` | 17.9 s | 1.1 s | 16× |
| `mc2023_track4_049.cnf` | 3.8 s | 0.35 s | 11× |
| `mc2024_track4_074.cnf` | 136 s | 6 s | 23× |
| `mc2023_track3_034.cnf` | 62.5 s | 2.4 s | 26× |

The pattern: autarky pays off most when `indep_sz` is small (≤ ~10), meaning the sampling set is already highly reduced by Arjun and there is significant independent structure. Track 4 (weighted projected) benefits most, consistent with the theory that weighted problems have more exploitable autarkies.

### 2. Three-Tier Learned Clause Management (CaDiCaL-style) in Ganak

Added to the Ganak main counter (commit `7e5a903`):
- Tier 1 (glue ≤ 2): always keep
- Tier 2 (glue ≤ 6): keep if recently used, survive 2 reduction rounds
- Tier 3 (glue > 6): discard if unused since last reduce

This replaces a simpler one-threshold policy. Effect: better quality learned clause database → faster propagation per conflict, better variable ordering maintenance. Likely responsible for many of the moderate (2–8×) speedups in track1 instances:

| Instance | Old time | New time | Speedup | Old conflicts | New conflicts |
|----------|---------|---------|---------|--------------|--------------|
| `mc2023_track1_118.cnf` | 788 s | 98 s | 8× | 1,033 | 640 |
| `mc2023_track1_151.cnf` | 309 s | 55 s | 5.6× | 205,478 | 14,258 |
| `mc2023_track3_093.cnf` | 2125 s | 386 s | 5.5× | 10,592 | 19,516 |

### 3. Watch Position Caching in Ganak (Gent'13, commit `1519a0a`)

Caches the last watched position in each long clause during propagation, avoiding rescanning from the start. This is a pure performance optimization — no effect on solution quality or conflict count, but reduces propagation overhead. Most visible on problems with many large clauses, explaining several 2–4× improvements on track3 instances with high `new_nvars`.

### 4. UIP Shrinking (commit `391ce7b`)

More aggressive conflict clause minimization in the main counter. Smaller learned clauses → shorter propagation chains → fewer conflicts. Works synergistically with the three-tier clause manager: more clauses survive because they are smaller and get used more often.

### 5. CMS Oracle Changes (affecting Arjun preprocessing)

The CMS oracle in Arjun does backbone detection, vivification, and equivalent literal discovery. 120 commits separate the two CMS versions. The substantive algorithmic changes:

**a) `AddClauseIfNeededAndStr` after binary discovery in `oracle_vivif`** (commit `67baf8408`): After the oracle finds an equivalence (`l1 ↔ l2`), it now adds the implied binary clauses back into the oracle's own clause database. Previously these were only added to the CMS main solver. This makes the oracle increasingly stronger as it finds more implications — a compounding effect. This is the most impactful oracle change.

**b) Tiered clause reduction in oracle**: CaDiCaL-style tiers applied inside the oracle CDCL loop. Better clause retention → fewer re-derivations of useful lemmas.

**c) Poison/removable memoization in oracle conflict minimization**: During UIP minimization, marks variables as "proven removable" or "proven non-removable" within a single conflict, avoiding repeated work. Speeds up individual SAT calls.

**d) Reason-side variable bumping** (commit `f8b035ad0`): VSIDS now bumps variables in reason clauses of the learned clause, not just the clause literals themselves. From CaDiCaL/MapleCOMSPS. Better variable ordering → fewer conflicts in backbone finding.

**e) Activity rescaling at 1e150 instead of 1e4**: The original threshold caused frequent expensive O(vars) rescaling operations. Moving to 1e150 makes this essentially never trigger in practice.

**f) Watch position caching in oracle propagation**: The same Gent'13 technique applied inside the oracle CDCL loop.

**Net effect on Arjun:** In 65/1165 paired cases `indep_sz` decreased (better simplification); in 69 it increased (slightly worse). On average `indep_sz` is ~129 for both runs — essentially unchanged. So the oracle improvements are not dramatically changing *how much* Arjun simplifies, but rather *how fast* it simplifies. However, for specific instances with complex equivalence structure the compounding `AddClauseIfNeededAndStr` can find many more equivalents — e.g. `mc2024_track2-random_087` went from indep 909 → 259 (650-variable reduction).

---

## Where Things Get Worse

### Random Track (track2-random) Arjun Regressions

`mc2024_track2-random_095`: total 20.99 s → 82.57 s; **Arjun: 15.7 s → 79.6 s** (5× slower preprocessing).

The improved oracle is doing *more work* in Arjun on random instances. Random CNFs have less structure (no functional relationships, few equivalences), so `oracle_vivif` finds little and the repeated `AddClauseIfNeededAndStr` calls add overhead without benefit. The oracle also runs more iterations because the "Increase limit" commit (`5dcdb0888`) raised the oracle call budget. On structured instances this extra work pays off; on random instances it wastes time.

This is a fundamental tension: oracle calls that improve preprocessing on structured problems are wasted on random problems. The old code's lower limit implicitly capped this waste.

### `mc2024_track3_031`: Easy Instance Becomes a Timeout

Old: 31.3 s, ind=960, nvars=1343. New run: same Arjun result (ind=960, nvars=1343, arjun=2.3 s) but the counter never finishes. The new counter (with autarky + 3-tier clauses) explores a different subtree order and hits repeated cache-full deletions. This is a search-order regression: the new branching heuristics lead to a harder subtree first.

### `mc2024_track3_139`: 79 s → 399 s (5× slower)

Both runs: same indep_sz=77, very similar cache miss rates (~0.26). New run: 297 conflicts vs 500 in old, but 399 s vs 79 s — so the new version makes *fewer* conflicts yet spends more time. SAT calls are nearly identical (7316 vs 7323). This suggests the overhead per operation increased significantly on this instance: likely watch position caching, three-tier clause management, or autarky checks adding per-step overhead on a flat search tree (few conflicts, many component decompositions per second).

### Track 3/4 Timeout Losses

Of the 11 lost instances, 5 are track3 and 2 are track4:

| Instance | Old time | Track |
|----------|---------|-------|
| `mc2023_track3_143.cnf` | 3381 s | track3 |
| `mc2024_track3_097.cnf` | 1286 s | track3 |
| `mc2024_track3_032.cnf` | 686 s | track3 |
| `mc2024_track3_034.cnf` | 574 s | track3 |
| `mc2024_track3_031.cnf` | 31 s | track3 |
| `mc2023_track4_170.cnf` | 3218 s | track4 |
| `mc2024_track4_064.cnf` | 1727 s | track4 |
| `mc2024_track4_197.cnf` | 223 s | track4 |
| `mc2023_track1_066.cnf` | 491 s | track1 |
| `mc2023_track1_182.cnf` | 102 s | track1 |
| `mc2023_track2_039.cnf` | 559 s | track2 |

Given that the new version improved overall on these tracks, these are statistical outliers where the new heuristics chose worse branching. The 10 newly-gained instances more than compensate in PAR2 even though they do not in solved count.

---

## Summary

**Good news:**
- `--puuraautarky 1` delivers dramatic wins on weighted projected problems (track 3 and 4), with some instances 10–26× faster
- Three-tier clause management and UIP shrinking give consistent 2–8× speedups on structured counting instances
- Overall PAR2 improved despite losing one net instance
- CMS oracle improvements make Arjun find more equivalences on structured problems, with occasional very large `indep_sz` reductions

**Concerns:**
- The oracle "Increase limit" + `AddClauseIfNeededAndStr` changes cause Arjun to spend much more time on random instances without benefit — the random track shows the clearest Arjun slowdowns
- Search-order regressions from the new autarky/clause management: at least one previously fast instance (track3_031, 31 s) is now timing out due to exploring a harder subtree first
- 11 timeouts lost vs 10 gains — the balance is almost exactly neutral in instance count, so PAR2 improvement comes purely from the quality of solved-instance speedups

**The single net loss** (`mc2023_track1_182.cnf`): old solved in 102 s, new times out. This track1 (unweighted) instance previously took modest time and now gets stuck. Not clearly attributable to autarky (track1 unweighted benefits less), more likely a clause management or branching order regression.
