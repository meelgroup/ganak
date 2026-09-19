# Branching heuristic: what the logs and experiments say (2026-09)

Context: `out-ganak-mccomp2324-1783906-0` (old TD) vs `out-ganak-mccomp2324-2345011-2`
(new TD chooser) differ mostly in the tree decomposition, and that alone moved
solved 1181 -> 1198 and the time on commonly solved instances 282.8Ks -> 242.0Ks.
All local experiments below: 2-core laptop-class box, `--maxcache 1500`, decisions are
the primary metric (CPU time is noisy there), counts always cross-checked.

## How a decision var is picked today

`score = activity/3 + td_weight * tdscore[v] + comp_frequency/25`

* `tdscore` is only the depth of the var's first bag counted from the centroid bag,
  normalized to [0,1]. All vars of one bag tie.
* `td_weight` is in practice binary: 7 (times up to 2 from the split) when
  nodes/width < ~8, 60 (120) when > ~10.

New stats (`br ...` lines, `br_*`/`td_levels`/`td_centroid_bag`/`td_soft_width` DB columns)
show two regimes:

* TD regime (low conflicts): TD obeyed ~100% of decisions, TD is 96-99.7% of the
  winning score, activity ~0%. Among the 3-25 vars tying on the TD score the
  component frequency alone decides.
* Flat regime (one bag, e.g. small projection sets): activity 97%, frequency 3%.

How coarse the TD signal is (1763 TDs in the new run): the centroid bag holds >25%
of ALL vars in 29% of the TDs (>50% in 13%), and 42% of TDs have <= 5 levels.

## Natural experiment in the two cluster runs (same instance, different TD)

182 instances with one TD each and the same node count:

* correlation of log(decision ratio) with log(width ratio): **0.02**. Width change
  does not predict the change in work.
* the biggest wins (x0.2-0.35 decisions: mc2023_track2_161, track3_098, track3_111,
  track1_115, ...; all dense, n=80-100, tw~50) came with TDs that were 6-12% WIDER
  but had **0.2-0.5x the levels**.

## What was tried, and what happened

Geometric mean of decisions vs. baseline on 27 quick instances, unless noted:

| change | result |
|---|---|
| `--tdsepwpct 100`: inside a TD level prefer vars in small adhesions in front of big subtrees | x1.27 on 11 inst., one x6.3. At 25%: x1.03. **Loss**: overriding the frequency score inside a bag hurts |
| `--tdmaxlevels 5` / `10`: merge TD levels | x1.24 / x1.02, up to x13 on sparse graphs. **Loss** when width/nodes < ~0.45, no effect above (already few levels) |
| `--freqscorediv 5`: 5x frequency weight | x1.08 |
| `--freqshortbonus 2` / `6`: more frequency score for clauses down to 2 unknown lits | x1.04 / x1.05 |
| `--cutvars 1 --cutw 50` / `500`: bonus for articulation vars of the component | x0.99 / x1.03, one instance 1.1M -> >13M decisions at 500 |
| **`--tdflatpct 50`**: TD does not guide branching when width >= 50% of nodes | **-17% summed time on the 13 instances in range, see below. Now default** |

Articulation vars are a rich signal that does not pay: on mc2023_track1_014 some
candidate would cut >= 10% off the component at 57% of all decisions (avg best cut 19.7%),
and we pick the best cut var only 34% of the time. But with component caching, picking
it late does not redo the cut-off part's work (it is a cache hit later), so decisions do
not drop, and forcing it overrides the TD order.

Conclusion: inside the TD regime the var choice is at a local optimum for every knob and
signal tried. The leverage is in WHEN the TD should be listened to, and in which TD.

## `--tdflatpct`: dense graphs

When the width is a large part of the graph, the centroid bag is most of the graph, the
TD cannot say where the graph comes apart, yet with `td_minweight` 7 it still makes up
~99% of the score. Base vs. TD ignored (time, decisions):

| instance | width/nodes | base | TD ignored |
|---|---|---|---|
| mc2023_track1_106 | 0.59 | 124.9s 45.2M | 69.3s 19.9M |
| mc2023_track1_097 | 0.53 | 103.2s 34.8M | 76.4s 24.7M |
| mc2023_track2_103 | 0.59 | 154.2s 37.0M | 119.8s 21.9M |
| mc2023_track1_109 | 0.59 | 90.8s 28.1M | 71.6s 23.2M |
| mc2024_track3_012 | 0.79 | 44.6s 3.5M | 28.5s 1.9M |
| mc2023_track3_080 | 0.57 | 180.2s 549K | 150.1s 472K |
| mc2023_track4_174 | 0.77 | 56.4s 7.7M | 54.6s 8.1M |
| mc2023_track3_091 | 0.55 | 148.1s | 147.6s |
| mc2024_track4_148 | 0.85 | 25.0s | 24.2s (same decisions) |
| mc2023_track4_138 | 0.54 | 16.5s 1.41M | 16.6s 1.60M |
| mc2023_track3_018 | 0.68 | 10.8s 937K | 10.5s 1002K |
| mc2023_track4_143 | 0.68 | 37.6s | 40.6s |
| mc2023_track4_083 | 0.63 | 15.9s 2.21M | 17.9s 2.49M |
| below the gate: mc2023_track2_117 | 0.46 | 122.6s | 139.4s |
| below the gate: mc2023_track4_085 | 0.47 | 23.4s | 28.0s |

In the new cluster run 224 of the 1117 instances with a TD have width/nodes in
[0.5, 0.9): 75 CPU-hours, 53 unsolved (above 0.9 the TD is one bag and flat anyway).
**This needs an A/B run on the cluster** (`--tdflatpct 0` is the old behaviour).

## Which TD property predicts the work? (8 quick instances x 4-5 distinct TDs)

Concordant : discordant pairs vs. decisions: width 23:8, soft width
`log2 sum 2^|bag|` 24:7, centroid bag size 22:8, levels 16:9, bags 16:14, centroid split
(largest comp after removing the centroid bag) **7:23** when forced on sparse graphs via
`--tddensepct 0 --tdbandpct 30`. So keeping the split band to dense graphs is right; a
better TD cost proxy than the width was not found at this sample size. Same-width TDs
still differ 2.4x in decisions (mc2024_track1_169: tw 22 flowcutter 873K vs tw 23
min-degree 371K).

## Found on the way

* The restart decay of `td_weight` clamped from below to `td_minweight`, RAISING a
  weight that was set to 0.1 on purpose. Fixed.
* `do_check_td_vs_ind` (`td_weight = 0.1` when the TD is wider than the indep support) was
  dead code: the clamp to `td_minweight` on the next lines undid it. Fixed (moved after the
  clamp), but now **off by default**, which is what all runs so far effectively had. The
  indep support here is Arjun's minimized one, so the rule fires widely: 106 instances of the
  cluster set that `--tdflatpct` does not already cover, including mc2024_track3_105/181
  where the TD is the big win. Off vs on, 9 of those: decisions geomean x1.10, time
  267s -> 276s, worst mc2024_track4_197 22.8s -> 32.4s (397K -> 924K decisions), best
  mc2023_track4_115 114K -> 79K.
  It compared the TD width over the opt-indep graph against the minimized indep set: 202
  of 1117 TDs trip it, 106 of them with width < 50% of the graph. Replaced by
  `--tdindtoppct N` (off by default): TD weight 0.1 when >= N% of the indep vars tie on
  the top TD score, i.e. the TD cannot order the vars the count is over. On the misfires
  above it stays quiet (track3_181: 32 indep vars over 4 levels, 9.4% on top).
  Removed after batch 2357294: on top of `--tdflatpct 50` it fired on only 57 TDs and
  solved 1202 vs 1206 without it; the flat gate already covers where it helped. Every run
  logs `[td] indep vars: ... pct: ... TD weight: ...`, DB columns `td_ind_*`/`td_weight`.
* `--tdlook` >= 0 (TD lookahead) was broken four ways, one behind the other: asserted
  `dl != -1` (probing from inside `decide_lit()` before the level's var is set), gave
  **wrong weighted counts** (`unset_lit()` multiplied the probed lits' weights into the
  level's count), a `setBag` assert (full-size graph given to the TD after contraction),
  and an empty `tdscore` when the toplevel TD was skipped. All fixed, 220 fuzz cases pass.
* The periodic stat line needed both 20M cache lookups AND 150K conflicts, so
  low-conflict timeouts left no stats at all. Either is enough now.
* `../count_fuzzer/fuzz.py` dies in arjun: `--bveplanner 5` is not accepted any more.
