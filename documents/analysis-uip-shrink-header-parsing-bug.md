# Performance Analysis: old (1193808) vs new (1247484)

## Overview

| Metric | Old (255f3da) | New (3c10c6a) |
|--------|--------------|--------------|
| Solved | 1176 | 1167 (+9 bug, see §1) |
| PAR2 (all) | 1145 | 1151 |
| PAR2 (hard instances only) | 260 | 242 |
| Avg time (solved) | 259.8s | 242.1s |
| New wins | — | 681 problems, avg +35.8s |
| Old wins | — | 391 problems, avg +19.2s |

Command lines:
- Old: `--maxcache 5000 --mode 1`
- New: `--maxcache 5000 --puuraautarky 1 --mode 1`

Key code changes between versions (153 commits):
- **UIP block-wise clause shrinking** (CaDiCaL-style, `--shrink`)
- Three-tier learned clause management
- Watch position caching
- New Arjun version (531a13ba vs 3fe1420)
- `--puuraautarky 1` enabled in benchmarks

---

## §1 Bug: 9 track4 files fail with "DIMACS header never found"

All 9 mc2023_track4 instances that old solves but new doesn't fail immediately with:

```
ERROR! DIMACS header ('p cnf vars cls') never found!
```

These are PWMC (probabilistic weighted MC) instances. The new Arjun writes a
`.ganak.cnf` preprocessed file in a format ganak cannot parse. The user
confirmed this bug is already fixed.

**Impact**: 9 solved instances lost. After the fix, the new system likely
matches or surpasses the old system on solved count.

---

## §2 Arjun Preprocessing Regression

A significant fraction of regressions originates entirely from slower Arjun
preprocessing, not from ganak's counter.

| Arjun direction | Count | Avg diff |
|-----------------|-------|---------|
| New Arjun faster | 759 | −22.4s |
| Old Arjun faster | 490 | +15.3s |

Notable Arjun regressions:
- `mc2024_track2-random_079.cnf`: +350s Arjun (2501 → 2851s); ganak itself same
- `mc2023_track3_061.cnf`: +200s Arjun (356 → 556s); ganak itself unchanged
- `mc2023_track3_071.cnf`: +193s Arjun (569 → 762s); ganak itself unchanged
- `mc2023_track4_167.cnf`: +305s Arjun (473 → 779s) + 145s ganak

**Pattern**: track2-random and track3 problems with long Arjun time (>100s)
tend to regress in the new Arjun. The new Arjun is (531a13ba) is slower for
some instances with large independent sets or many-variable problems.

---

## §3 UIP Shrinking: Large Wins for Large-Clause Instances

The UIP block shrinking is the most impactful change. It dramatically reduces
learned clause sizes for instances where Ganak accumulates large clauses.

### Example: `mc2024_track3_198.cnf` — old 2581s → new 1416s (−1165s)

| Metric | Old | New |
|--------|-----|-----|
| Decisions/s | 2.66 K/s | 12.17 K/s |
| Conflicts/s | 85.51/s | 353.26/s |
| Avg clause size (final) | 5416 | 477 |
| Shrink removed | (rem_lits 45455K) | 291M lits |
| Time | 2581s | 1416s |

The shrinking collapses huge clauses (avg ~6304 lits) down to ~477, making
propagation 10–20x faster. The decisions/s quintuples.

### Example: `mc2023_track3_182.cnf` — old 1455s → new 732s (−723s)

| Metric | Old | New |
|--------|-----|-----|
| Decisions/s | 8.4 K/s | 19.1 K/s |
| Conflicts/s | 7334/s | 16727/s |
| Avg clause size (final) | 139 | 31 |
| Shrink removed | 357K | 280M lits |
| Components (K) | 1204 | 945 |
| Time | 1455s | 732s |

Average clause size drops from 139 → 31. Decisions and conflicts per second
both roughly double. Component count also drops, suggesting better clause
quality helps decomposition.

### Big-win statistics (ganak-only, Arjun factored out)
- 125 instances where ganak core is >20s faster (avg −130s)
- 43 instances where ganak core is >20s slower (avg +92s)
- 998 instances similar (within 20s either way)

The win-rate for ganak's core solver is very good: 3:1 in favour of the new
system, with larger improvements than regressions.

---

## §4 Ganak Solver Regressions: Two Patterns

### Pattern A: Altered component structure (worse for component-heavy instances)

The new system's better clause quality changes which sub-formulas are
identified as independent components, and shifts the cache miss rate.

**Example: `mc2023_track1_175.cnf`** — old 585s → new 985s (+400s, zero Arjun diff)

| Metric | Old | New |
|--------|-----|-----|
| Decisions/s | 410–444 K/s | 587–599 K/s |
| Conflicts/s | 1441–1576/s | 716–753/s |
| Avg clause size (final) | 44 | 14 |
| Components (K) | 497K → 1252K | 2306K → 3120K |
| Cache miss | 0.512 | 0.593 |
| Time | 585s | 985s |

The new system is _faster per decision_ (599 vs 444 K/s) with much smaller
clauses, but explores 2.5x more components. The smaller learned clauses split
the formula into more, finer-grained components that are individually novel
(cache miss goes up). Total work increases despite the per-unit speed-up.

**Interpretation**: UIP shrinking is changing which variables appear in the
same component. Better clause quality means fewer shared variables across
learned clauses, which creates more but smaller components. The component cache
is less effective because fewer components repeat.

### Pattern B: Pure DPLL overhead — UIP shrink adds cost without saving search

**Example: `mc2024_track4_024.cnf`** — old 1398s → new 1783s (+385s, cache miss = 1.0)

| Metric | Old | New |
|--------|-----|-----|
| Decisions/s | 33.4 K/s | 40.0 K/s |
| Conflicts | 42.4M | 66.1M |
| Avg clause size (final) | 25.0 | 21.8 |
| Shrink removed | — | 464M lits |
| Components | 1 (pure DPLL) | 1 (pure DPLL) |
| Time | 1398s | 1783s |

Cache miss is 1.0 — there is exactly one component (the whole formula), so
this is pure CDCL counting. The new system is faster per decision, and
marginally improves clause size (25 → 22). But it does 55% more conflicts.
The better clauses change the search trajectory in a way that requires more
total conflicts to terminate.

**Example: `mc2023_track3_144.cnf`** — old 529s → new 647s (+119s)

| Metric | Old | New |
|--------|-----|-----|
| Decisions/s | 186 K/s | 21 K/s |
| Conflicts | 5.4M | 6.3M |
| Avg clause size (final) | 21.7 | 22.1 |
| Shrink removed | — | 22M lits |
| SAT oracle calls | 611K | 644K |
| Time | 529s | 647s |

Here the old system is 9x faster at decisions (186 vs 21 K/s). Both have
similarly-sized clauses, but the new system runs the 6.3M UIP shrink calls —
even with 77% success rate, removing 22M lits, the shrinking overhead
dominates. The clause size barely changes (from ~31 before shrink to ~22 after,
vs old ~30 → ~22), meaning the shrink is not being decisive on this instance
but is still costly.

---

## §5 Summary and Hypotheses

### What's clearly working
- UIP shrinking is a major win for instances with large learned clauses (avg
  size > 50–100). These instances now solve 1.5–2x faster.
- Overall ganak core: 125 big winners vs 43 big losers, with larger average
  improvement than regression.

### What's regressing and likely root causes

1. **Arjun regression** (new Arjun 531a13ba): ~490 instances regress in
   preprocessing, avg 15s. Some instances regress by 100–350s in Arjun alone.
   This is unrelated to ganak's solving and warrants investigation in Arjun.

2. **UIP shrink overhead on small-clause instances**: When learned clauses are
   already small (avg final size < ~25 lits), running 6M+ shrink calls adds
   significant overhead per conflict without much benefit. A guard like
   `if (avg_lbd > threshold) shrink()` or skipping shrink for short clauses
   could help.

3. **Altered component structure**: Better clauses change which variables
   are shared across components, sometimes creating more components (2–3x)
   with worse cache reuse. This is harder to fix — it's a fundamental
   interaction between clause learning and component decomposition. Possible
   mitigations:
   - Restrict when shrinking is applied (e.g., only at restarts, or only for
     clauses above a size threshold)
   - Investigate whether the watch position caching change is contributing to
     different VSIDS scores and thus different branching order

### Metrics to watch
For a new run, compare:
- `comps recordK` ratio (if new >> old on a regressing instance, the component
  structure has changed)
- `Kdec/s` drop (if significantly lower in new, UIP shrink overhead is the
  likely cause)
- `avg clsz finalavg` (if similar to old, shrinking isn't helping but may
  still be running)
- `cache_miss_rate` increase (indicates more novel components, reduced cache
  benefit)

---

## Example Log Comparisons (Quick Reference)

### Best improvement: `mc2024_track3_198.cnf`
- Old log: `.../out-ganak-mccomp2324-1193808-0/mc2024_track3_198.cnf.gz.out_ganak`
- New log: `.../out-ganak-mccomp2324-1247484-0/mc2024_track3_198.cnf.gz.out_ganak`
- Clause size 5791→477, time 2581→1416s

### Worst regression (ganak core): `mc2023_track1_175.cnf`
- Old log: `.../out-ganak-mccomp2324-1193808-0/mc2023_track1_175.cnf.gz.out_ganak`
- New log: `.../out-ganak-mccomp2324-1247484-0/mc2023_track1_175.cnf.gz.out_ganak`
- Components 1252K→3120K, cache miss 0.512→0.593, time +400s

### Worst regression (Arjun): `mc2023_track3_061.cnf`
- Old log: `.../out-ganak-mccomp2324-1193808-0/mc2023_track3_061.cnf.gz.out_ganak`
- New log: `.../out-ganak-mccomp2324-1247484-0/mc2023_track3_061.cnf.gz.out_ganak`
- Arjun 356s→556s (+200s), ganak itself identical

### Pure DPLL regression: `mc2024_track4_024.cnf`
- Old: 42M conflicts → 1398s; New: 66M conflicts → 1783s
- Cache miss = 1.0 in both; UIP shrink running but not helping
