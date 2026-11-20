# Restarts: what the code does, what was measured (2026-11)

develop has `--bjdefer` (default 1) and the cache-keeping restart as the only restart
(no `--rstmode`; the cube restart is only the AppMC timeout fallback). `--rstmode`,
`--rsttype 10`, `--progressevery`, `--bjkeep`, `--satprune` and `--rstpartial` exist only
on branch `ganak-restart`.

Local experiments: 16-core box, `--maxcache 1500`, 300s timeout, several configs run
concurrently so they share the load. Sets, all from run `out-ganak-mccomp2324-2887638-0`:

* setC: 100 instances with >= 30K conflicts that took 5-250s on the cluster
* setD: the 51 of them that took <= 120s
* conflicts/decisions/cache stores are deterministic, so they are compared too (geomean
  ratio over the instances both configs solved, `work.py`)

## The old restart (`--rstmode 0`, now only the AppMC fallback) threw everything away

At a restart, every decision level's counted branches became "cubes" with counts, their
negations were added as clauses, and the counter was rebuilt (`init_and_preproc()` makes a
new `CompManager`, so the component cache and the frequency scores went too). Since #97
`compute_cube()` pins *every* indep var from a SAT model, so each cube counts exactly one
model: a restart kept ~2 models per decision level and lost all other counted work, and the
cache. `extend_cubes`/`cube_resolve`/`cube_flp` cannot win that back. Only the AppMC
timeout fallback still uses this path, where it does not matter.

## Restarts that keep the cache (`--rstmode 1`, the only restart on develop)

Backtrack inside the search, keeping every finished component's cached count. Sound only if
no component that is not counted yet is UNSAT: learnt clauses over lits outside a finished
component can prune its models when the rest is UNSAT under the trail, which is why
abandoning a branch normally drops all that was counted in it. So one CMS call first checks
that the trail's decisions extend to a model; if CMS says UNSAT, the failed assumptions
show the branch of their highest level is UNSAT, and we backtrack out of it like a
conflict would (`--satprune N` also does this check every N conflicts).

On setC, with the old readjust policy (TD on/off alternating at each restart) for both:

| | solved | time on the 50 common |
|---|---|---|
| no restarts | 94/100 | 4335s |
| `--rstmode 0` | (50 common) | 6296s |
| `--rstmode 1` | 68/100 | 5442s |

The readjust alternation is most of the damage: after restart #1 `decide=1` (TD ignored)
for the rest of the run. Without it (`--rstreadjust 0 --tdwrstdecay 1`), setD:
**51 -> 43 solved, time x1.34, decisions x1.87, not one instance faster.**
`--rstpartial` (restart only from the lowest level whose decision the scores no longer
make) changes nothing: the level-1 choice always differs, activity is halved every 256
decisions. With `--rstmargin 2` (must beat it by 3x) half the restarts are skipped and it
is still 5% more cache stores on mc2023_track4_100.

## Why restarts cannot win much on runs that finish

* In SAT a restart loses almost nothing. In a counter the partly counted path is the asset,
  and only what the cache gives back survives a change of the top decisions. That is little:
  top components are unique, mc2024_track2-random_115 has a 99.6% cache miss rate anyway.
* ganak's cost is not heavy-tailed in the branching order. `--shufflevars 1` vs not
  (2425144-0/1): >2x apart on 23 of 772 instances (>20s), >4x on 5. `--tdflatpct 0` vs 50
  (2357294-0/1): >4x on 3 instances. The 30 instances the virtual best of 73 runs solves
  but the default does not are almost all preprocessing changes/regressions.

## Progress estimate and `--rsttype 10`

`search_progress()`: Knuth-style, a decision halves its branch's weight, a finished left
branch is done, a branch's weight is shared by its comps by their number of vars. Printed
with the stat line, `--progressevery N` every N decisions. setC, estimated / actual total
decisions:

| at % of search | done/progress: p10 / median / p90 / max | min(linear, log-linear): p10 / median / p90 / max |
|---|---|---|
| 10% | 0.18 / 1.20 / 90 / 32579 | 0.18 / 0.25 / 4.1 / 32462 (too few samples) |
| 25% | 0.37 / 1.17 / 8.2 / 5226 | 0.32 / 0.45 / 0.77 / 2.8 |
| 50% | 0.72 / 1.00 / 4.9 / 658 | 0.59 / 0.79 / 1.15 / 1.4 |

Progress often grows exponentially with the work (later parts go faster: learnt clauses,
cache), so done/progress is very pessimistic early; extrapolating log(progress) between two
checks mostly fixes it. The first `--rsttype 10` restarted when even this optimistic
estimate was >= 50x the work done (checked from 2M decisions on, at doubling intervals),
flipping whether the TD guides branching. On setD it
never fires and is neutral (2304s vs 2307s). The cluster data show only mc2024_track4_170
where TD-on beats the default's TD-off by >4x.

On 32 instances the default does not solve in 3600s (`setU`, 600s each), the first version
(switch on done/progress alone, at 2M decisions) switched on 30. The last progress estimate
did not move in the median (log10 ratio -0.04), but on some it fell by orders of magnitude:
mc2024_track2-random_173 6.7e10 -> 5.7e7 decisions, random_177 6.8e10 -> 7.3e7,
mc2024_track3_130 6.3e8 -> 4.7e7, mc2024_track4_135 3.3e10 -> 1.7e8, mc2024_track1_103
5e31 -> 1.4e9; on others it rose as much. Nothing was solved in 600s either way. On setC,
1 of the 32 instances past 2M decisions would have switched on done/progress alone
(mc2023_track3_099, which the cluster solves).

The extrapolations were then too far off to decide on (on random_173 the log-linear one
said 17M decisions in total at 6M, then progress stalled for 20 minutes), so the policy is
now a single probe compared at the same amount of work: at 2M decisions, if < 2% of the
tree is done (`--rsthopeless 0.02`), restart into the other strategy (`--rstprobealt 1`,
or the same one with 0); after as many decisions keep it, unless the first strategy's
progress extrapolated linearly to as much work is ahead. At most two restarts.

* **The restart, not the strategy, is what helps** on random_173/177 (unsolved in 3600s on
  the cluster): progress after ~2.4M decisions is 2.4e-5 before the restart and 4.8e-2
  (same strategy) / 4.2e-2 (TD off) after it, while the default and `--decide 1` (TD off
  from the start) are at 7.7e-4 / 5.5e-4 after 7.3M. After 3400s the restarted runs are at
  26% / 29% of the tree (~7e7 decisions in total), the default estimates 1e10. Not solved
  in 3400s here; whether it adds solves needs a cluster run.
* **On solvable runs it costs**: forced on the 25 setC instances past 2.5M decisions,
  decisions x1.35, time x1.25 (17 of the 22 that restarted went back). The 2% gate fires on
  1 of 32 such instances.

### Cluster run 2893602: `--restart 1` (`--rsttype 10`) loses

1600 instances, 3600s. Default 1218 solved, `--restart 1` 1208 (+3 / -13), geomean time
x1.059. The gate is not rare: it fired on 372 of 1600 (not 1 of 32 as on setU).

| | instances | gained/lost | >=1.5x slower/faster | geomean |
|---|---|---|---|---|
| no restart | 1228 | 0 / 0 | 0 / 0 | 0.997 |
| restarted, probe went back | 152 | 0 / 0 | 1 / 0 | 1.024 |
| restarted, kept the other strategy | 220 | 3 / 13 | 39 / 2 | 2.03 |

The cache-keeping restart itself costs ~2%. All the loss is the probe picking the
wrong strategy (both ways: 134 into TD-off, 86 into TD-on). Progress at equal work
misleads: mc2023_track1_118 at 2M decisions shows 2.6e-3 TD-off vs 5.9e-3 TD-on, but
TD-off finishes in 31M decisions and TD-on in 381M. Do not retry this probe.

## Backjumps were the real waste

All cache pollution ganak drops comes from `go_back_to`, none from UNSAT branches in
`backtrack()`, and ~95% of it at the backjump's *target* level: the finished siblings of
the comp the conflict was in. On setC an average of 36% of all cache stores were dropped
as pollution, >30% on 44% of the instances.

Those siblings are valid unless the target level's branch is UNSAT. `--bjkeep N` asks CMS
(when >= N entries would be dropped): setD decisions x0.938 / x0.946 / x0.962 for N = 1 / 5
/ 20, time x1.036 / x0.977 / x0.983 (the SAT calls). `--bjdefer 1` keeps them for free:
the level remembers where the skipped truncation started (`StackLevel::kept_from`), and
removing its pollution later truncates from there. **Default since this work.** setC, both
98/100 solved, same counts: decisions x0.901 (track1 x0.890, track2 x0.851, track2-random
x0.925, track3 x0.881, track4 x0.956), counting time (total minus Arjun) x0.969 geomean,
5382s -> 5145s summed; mc2024_track1_174 5.22M -> 4.03M decisions, 60.4 -> 50.8s.
Checked with `CHECK_COUNT`+`SLOW_DEBUG` (cache on) on ~200 fuzz CNFs, with up to 69
evictions per run, and `tests/cnf-files/restart_bj.cnf`.
What is still dropped with `--bjdefer` (22.7% of stores on setC, from 36.7%) is needed:
levels above the conflict's, and the target of a flipped decision, whose branch the learnt
clause shows UNSAT.

`--satprune` alone: on mc2023_track4_100, 112 UNSAT trails in 956 checks, jumping 0.8
levels on average: ganak's own CDCL finds those conflicts about as soon as CMS does.

## More or less dynamic branching

`--vsadsadjust` (halve all activities every N decisions), setD, all 51 solved by each,
vs 256: 64: decisions x1.013, conflicts x1.081, time x1.002; 1024: x1.042, x0.960,
x1.020; 4096: x1.059, x0.884, x1.034. Slower decay trades conflicts for decisions, and
loses. 256 stays.

`--actscorediv` (activity's share of the score), setD vs 3: 1: decisions x1.053; 10:
x0.986, conflicts x1.033; 30: x1.009, x1.121; 100: x1.061, x1.253. 10 is a shallow
optimum, worth a cluster A/B, not a default change on 51 instances.

## Not done / ideas

* Cluster A/B of `--restart 1` (the probe): it only fires when < 2% is done at 2M
  decisions, so the solvable set should not move; the question is the timeouts.
* Cubes made only of the top "spine" (levels with a single comp above them) would be large
  and exact, and would keep that work across a restart whatever the new order is.
