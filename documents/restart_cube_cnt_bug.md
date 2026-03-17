# Bug: Incorrect `c.cnt` in `compute_cube` for Non-Weighted Counting

## Summary

`compute_cube` computed a wrong `c.cnt` for non-weighted counting modes (`--mode 0`).
The root cause was that `c.cnf` did not cover all required indep-support vars, leaving
some vars free in the cube. Because the blocking clause only pins the vars in `c.cnf`,
it banned more models than `c.cnt` claimed, causing an undercount (and an overcount
after an attempted fix via ancestor-multiplication).

Two additional bugs were found in the process:
- `include_solution` for `!is_indep` levels did not propagate UNSAT (zero count)
  correctly: it always set `branch_mc = fg->one()` even when `solutions = nullptr`.
- `init_decision_stack()` passed the (stale) `is_indep` member variable to the root
  level 0, which could leave that level with `is_indep = false`. The root level must
  always be `is_indep = true` so it properly accumulates the final projected count
  (especially important for the EXIT-state cube count with `--satsolver 0`).

---

## Background: How Restart Cube Collection Works

At restart time (`restart_if_needed`), Ganak walks the decision stack from the deepest
level to level 1. For each level and branch side (0 = left, 1 = right), it calls
`compute_cube(cube, side)` to capture the counted subspace. Each cube is then:

1. Added to `mini_cubes` with its count `cube.cnt`.
2. Converted to a blocking clause `sat_solver->add_clause(ganak_to_cms_cl(c.cnf))`.
   This clause is satisfied when at least one pinned var differs from its cube value,
   so it bans exactly those models where ALL pinned vars match the cube assignment.

The invariant: **`cube.cnt` must equal the number of models banned by the blocking
clause**. Too small → permanent undercount. Too large → overcounting (double-count).

`check_count_norestart_cms(cube)` independently verifies the count by forcing the
cube vars as unit clauses and enumerating satisfying indep-support assignments.

---

## The Bug: Incomplete Cube Coverage

`compute_cube` builds `c.cnf` from:
1. **Decision literals** — all trail decisions for vars `< opt_indep_support_end`.
2. **Unprocessed component vars** — indep vars (`< indep_support_end`) from component
   ranges `[remaining_comps_ofs_, unprocessed_comps_end_)` at each decision level,
   taken from the current SAT model.

**What was missing**: Indep vars (`< indep_support_end`) that are:
- In already-*processed* sub-components (not in any unprocessed range), AND
- Not in the trail as decisions (e.g., implied/propagated, or solved by cache hit).

These vars ended up free in `c.cnf`. The blocking clause then banned all satisfying
assignments for those free vars, but `c.cnt` did not account for them.

### Concrete failing example (`fuzzTest_2560.cnf`, `--mode 0 --td 0 --rstfirst 3 --restart 1`)

```
indep_support_end: 17    (indep vars = 1..16)
Cube CNF: 20 -8 -9 -4 10 -1 -2 -3 5 -11 -12 -13 -14 -15 -16
```

Vars 6 and 7 were absent from the cube — they lived in already-processed sub-components
and were implied rather than decided. The blocking clause banned all 3 satisfying
assignments for `{6, 7}` given the 14 pinned vars, but `c.cnt` was wrong.

**Original code** (`c.cnt = fg->one()`): claimed 1, actual was 3 → undercount by 2.

**First fix attempt** (ancestor `branch_mc` multiplication): `96 × 3 × 4 × 2 = 2304`.
Wrong because `branch_mc` values are *aggregate* counts accumulated over many different
assignments in the DPLL tree — not the count for the one specific assignment pinned
by the cube. The specific assignment from the SAT model was already fully pinned for
all the processed sub-comp vars (they were in the trail as decisions), so each
contributed 1 to the actual count, not their historical aggregate count.

---

## Fixes

### Fix 1 — Complete `c.cnf` coverage (`compute_cube` in `counter.cpp`)

After the dedup pass, scan all vars `v` in `[1, indep_support_end)` and add any not
yet present, taken from the SAT model:

```cpp
{
  set<uint32_t> seen_vars;
  auto it = std::remove_if(c.cnf.begin(), c.cnf.end(),
      [&seen_vars](const Lit& l) { return !seen_vars.insert(l.var()).second; });
  c.cnf.erase(it, c.cnf.end());
  for (uint32_t v = 1; v < indep_support_end; v++) {
    if (seen_vars.count(v)) continue;
    c.cnf.push_back(Lit(v, sat_solver->get_model()[v-1] == CMSat::l_False));
  }
}
```

With all indep vars pinned, the blocking clause bans exactly one indep-support
assignment. For non-weighted, `c.cnt = fg->one()` is now correct.

### Fix 2 — UNSAT propagation in `include_solution` (`stack.hpp`)

The old code for `!is_indep` levels:

```cpp
if (!is_indep) branch_mc[act_branch] = fg->one();   // BUG: ignores UNSAT children
else { /* proper accumulation */ }
```

For non-indep levels, Ganak's projected counting intentionally tracks only satisfiability
(0/1), not the exact count — the non-indep variable assignment does not affect the
projected count over indep vars when no indep vars remain in the component. However,
the old code incorrectly set `branch_mc = fg->one()` even when `solutions = nullptr`
(child returned UNSAT). This meant a branch where the child component was UNSAT would
be treated as SAT, causing over-counting.

Fix: propagate zero correctly for the UNSAT case, while preserving the 0/1 invariant:

```cpp
if (!is_indep) {
  if (cnt_is_zero(solutions)) {
    branch_zero[act_branch] = true;
    branch_mc[act_branch] = nullptr;
  } else {
    branch_mc[act_branch] = fg->one();
  }
} else {
  /* existing indep accumulation logic */
}
```

### Fix 3 — Root level always `is_indep = true` (`init_decision_stack` in `counter.cpp`)

`init_decision_stack()` pushed the root level (level 0) with the stale `is_indep`
member variable:

```cpp
decisions.push_back(StackLevel(1, 2, is_indep, tstamp, fg));  // BUG: is_indep may be false
```

Level 0 is a dummy root that accumulates the final total count from all level-1
sub-problems. It must always have `is_indep = true` so that `include_solution` properly
accumulates child counts (products of component sub-counts) rather than collapsing to
0/1. When `is_indep = false`, `include_solution` discards the actual count and stores 1,
making `total_model_count()` at the EXIT state return 1 regardless of the true count.
This was the reason a `check_count_norestart_cms` fallback existed at the EXIT state —
it patched over the wrong `total_model_count()` with a slow one-by-one enumeration.

Fix: always use `is_indep = true` for the root level:

```cpp
decisions.push_back(StackLevel(1, 2, true, tstamp, fg));  // root always accumulates count
```

With this fix, the EXIT-state cube's `c.cnt = decisions.top().total_model_count()` is
correct even for `--satsolver 0`, and the `check_count_norestart_cms` production
fallback (which was the slow debug-only enumeration) is no longer needed.

---

## Steps to Find the Bug

1. **Initial report**: performance bug — some counted subspaces were not collected as
   cubes, causing re-exploration after each restart.

2. **Code review** of `restart_if_needed` and `compute_cube`: found `c.cnt = fg->one()`
   applied unconditionally for non-weighted, regardless of how many models the cube's
   partial assignment covers. Wrote analysis to `documents/restart_cube_analysis.md`.

3. **First fix attempt**: removed the `c.cnt = fg->one()` override; used ancestor
   `branch_mc` multiplication instead.

4. **Test failure**: `CHECK_COUNT` assertion triggered:
   ```
   recorded c.cnt : 2304   (= 96 × 3 × 4 × 2 from branch_mc products)
   actual check   : 3
   ```
   Debug showed: all "excluded comp" vars were covered (0 uncovered), free indep vars
   were `{6, 7}` with 3 satisfying assignments. `branch_mc` aggregate values were wrong
   because they represent counts over many different DPLL-tree assignments, not the one
   specific SAT model assignment pinned by the cube.

5. **Root cause identified**: the cube was not covering all indep vars. Vars 6 and 7
   (in processed sub-components, not decided) were absent from `c.cnf`. Under the cube's
   partial pinning, the blocking clause banned 3 models but `c.cnt` claimed 1 (or 2304).

6. **Correct fix for compute_cube**: pin ALL indep vars from the SAT model in `c.cnf`
   during cube construction. With full coverage, `c.cnt = fg->one()` is exact (Fix 1).

7. **Second regression** (`fuzzTest_2574.cnf`): `assert(*decisions.top().total_model_count() == *cnt_cmp)` failed in `check_count()`. This led to investigation of `include_solution`
   for `!is_indep` levels and the UNSAT propagation bug (Fix 2), plus the discovery that
   `init_decision_stack()` was initializing the root level with `is_indep = false` (Fix 3).
   The slow `check_count_norestart_cms` fallback in the EXIT state was a symptom of this
   root-level `is_indep` bug and was removed once Fix 3 was in place.
