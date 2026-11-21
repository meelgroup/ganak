# Speed ideas tried on 2026-11-21 (branch develop-new-ideas)

Local box: 16 cores, `--maxcache 1500`. With 16 jobs running, time on an identical
search moves +-5-25%, so timings alone decide nothing here. `perf record` needs
`perf_event_paranoid` -1, but `perf stat -e instructions:u,cycles:u` works at 2, and
instructions are load-independent. Profiles: `scripts/gdb_sample.py`.

## Where the time goes (run 2898948-0)

* Solved instances: Arjun 101.6Ks of 259.5Ks (39%). In the 60-1000s range it is about half.
  319 instances that take 20-300s spend more than half their time in Arjun. Only 82 that
  take 20-120s spend more than 60% of it counting.
* In those 319: Puura 12Ks (its oracle stage 6Ks, backbone 2.5Ks), minimize/backward 8Ks,
  cover-drop 4.3Ks.
* Counting: the 245 instances with `br_multi_pct` < 5 (components almost never split)
  take 309Ks of counting time, 28% of it all, at an 86% cache miss rate. The 137 instances with
  more than 20K SAT calls take 167Ks.
* Profile of one such instance (mc2023_track3_051, 25 indep vars, 141K SAT calls,
  139K projected models): 46% in SAT mode, 41% in component analysis. Profile of an
  unweighted track1 instance: about 15-20% is GMP/FMpz alloc+free around StackLevel counts.

## Kept

* fc795010: StackLevel reuses its count objects, and backtrack() computes the total once.
  Instructions -0.5% to -4.1%, most on weighted mpq.
* 4ded9177: `-fno-semantic-interposition`. Everything is -fPIC, so GCC did not inline
  exported functions inside their TU. Instructions -0.1% to -3.1%.
* 2e445928: `Counter::weighted()` caches `fg->weighted()`. That is a virtual call, and
  set_lit()/unset_lit() ask it on every assignment. Instructions -0.2% to -0.7%.

Check of fc795010 + 4ded9177 against develop (3489117e): the 91 local instances (32
Arjun-heavy, 32 count-heavy, 27 with > 20K SAT calls). Decisions, conflicts and counts
are identical on all 91. User instructions for the whole run, Arjun included: geomean
0.980, best 0.81 (mc2023_track2_083). The next six best, all count-heavy track1/2, are
at -4% to -5%. ctest passes, and count_fuzzer `--only 40/60/40` found no mismatch.

## Option added, default unchanged: `--arjunoraclemult2` (needs arjun develop-new-ideas)

Puura's second pass in `standalone_elim_to_file` (after cover-drop) re-runs the oracle on
a CNF the first pass already did. Run 2898948-0, solved instances: pass 1 sparsify 7.4Ks
+ vivify 8.0Ks removes 688K clauses and 3.9M lits; pass 2 sparsify 6.0Ks + vivify 4.1Ks
removes 92K clauses and 327K lits. Most of the oracle time is sparsify; on
mc2023_track1_050 pass 2 takes 16.6s to remove 112 clauses.

Local, whole-run user instructions (counts equal everywhere):

| set | `--arjunoraclemult 0.25` (both passes) | `--arjunoraclemult2 0.25` (pass 2 only) |
|---|---|---|
| 16 instances where the oracle stage takes > 20s | x0.647; counting 159 -> 255s, CNFs bigger (2814 -> 4104 cls) | x0.846; CNFs close (2814 -> 2924 cls) |
| 48 random solved, 10-200s | x0.931; counting 1506 -> 1548s, mc2023_track1_127 counting 28.6 -> 47.9s | x0.972; CNFs close (2050 -> 2147, 8636 -> 8818 cls), no counting loss seen |

Timings were taken under 16-32 jobs, so only the instructions and CNF sizes are firm.
Pass 2 at 0.25 looks like a default candidate, but it changes the CNF, so it needs a
cluster run first.

## Tried, dropped

| idea | result |
|---|---|
| SAT mode: pick the first 8 decisions by linear scan, build the VSIDS heap only after that (the heap was 30% of SAT-mode samples) | instructions +1.1/-0.9/+0.2/+2.0%; quiet count time +3.5..-4.9%, mean -0.7% |
| Right branch of the last projected decision counted 1 if flipping it in the left branch's model satisfies the component | 0 of 51762 tries succeed: binary clauses force non-indep vars to change, so this needs propagation, which is most of the branch's cost |
| No counted var left in the comp: count it 1, no analysis or SAT call, if trail + saved phases satisfy all its clauses | right count, but it holds in 701 of 143505 tries on mc2023_track3_051 (0.5%) |
| Pass the level total to the cache store and the parent by pointer (no dup when one side is zero) | instructions -0.03% |
| Brute-force components of <= 10 vars (Cara) | <= 10-var components are < 1% of cache misses on 5 of 6 count-heavy instances, 16% on one |
| d4's untouched-clauses-out-of-the-key | for a fixed var set the untouched clauses are fixed, so this changes no hit, it only shortens the hash input. Not built |
| `--arjunoraclemult 0.3` | 50 instances: total time 2862 -> 2914s, noise |
| Arjun backward with CMS chrono BT (`diffdeclevelchrono` 100 / 20) | mc2023_track1_146 (18496 -> 9 vars): 194.8 / 194.8 / 172.9s at -1 / 100 / 20, under load. Backjumps are short, so it is not assumption re-propagation. `--maxc 500`: 128 indep vars, not faster |
| Stop flowcutter after K iterations with no gain over min-degree | saves 2.8Ks over all 1600 at K=100, loses width on 29 TDs (167 -> 143 on one) |
| No flowcutter when min-degree is already >= `--tdflatpct` | 2.2Ks over all 1600, and flowcutter gets below the gate on 20 of the 354 |
| Cover-drop: give up after R rounds with 0 no-Y cubes | R=1000 saves 2.2Ks, kills 24 pure-drop successes; most give-ups do emit cubes |
| Skip backward when it ends with I == U | 6.3Ks overall, but only 391s on solved instances |
| `-fno-semantic-interposition` in CMS too | instructions -0.0 to -0.5% on 5 Arjun-heavy instances |
| Hoist loop invariants out of find_best_branch / record_comp (member stores, `last_dec_candidates`, `dec_level()`) | instructions -0.15% / -0.8..+0.4% |
| `--analyzecand` 5 / 500 instead of 50 | instructions +3.1% / +0.9% on mc2023_track1_097, +-0.5% elsewhere |
| Skip cover-drop on low-width formulas | in run 2898948-0 at tw <= 30: 45 give-ups cost 1.9Ks, and the 37 that fired gained nothing over cover-off (2634304-0). That is ~0.2% of the total, and the width is not known yet when cover-drop runs |
| Skip Puura's second pass when autarky and cover-drop changed nothing | it still changes the formula in 360 of those 902 solved instances (mc2024_track3_062: 45.7K -> 33.6K cls) |

## Found

* CMS `oracle_sparsify()`: `ccnr.run()` is outside the mems budget. On mc2023_track3_166
  (1.5M binary clauses) each 10% of the pass takes 460s while oracle mems stay at ~60M per
  10%. That instance never gets past it, in any run. It is the only such instance in the
  run.
* Arjun `set_up_for_arjun()` turns off chrono BT, and `diffdeclevelchrono` is not a
  `set_option()` (lib) option.
* Competition binaries on instances where ganak takes > 300s or times out, at 300s:
  gpmc solved 0/16 track1 (most hit the 6 GB limit) and 1/16 track3 (mc2023_track3_057,
  199s vs ganak's 362s). d4 and SharpSAT-TD both time out on mc2023_track3_111 (80 vars,
  tw 51, 128M decisions, ganak 303s). gpmc `-mode=2` reads a file without a show line as
  an empty projection, so its 0.01s/count 1 there means nothing.
* After the kept changes, decide_lit (incl. find_best_branch) went from 19% of the samples
  on mc2023_track1_097 to 8%. What remains is component analysis and the cache (70%). Its hash compares
  are the first touch of the cache entry, so they are memory latency.
