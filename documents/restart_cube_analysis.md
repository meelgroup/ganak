# Restart Cube Collection Analysis

## Summary

The cube collection loop itself has **no performance bug** — all non-zero-count branches are
attempted, and every successful `compute_cube` call is added to `mini_cubes`. However, there is a
**correctness bug** in the count assigned to each cube for non-weighted counting.

---

## The cube collection loop is correct

In `restart_if_needed` (counter.cpp:1398–1432):

```cpp
while (dec_level() > 0) {
    for(auto i: {0, 1}) {
        const FF& models = decisions.top().get_model_side(i);
        if (models == nullptr || models->is_zero()) continue;
        Cube cube;
        if (compute_cube(cube, i)) {
            mini_cubes.push_back(cube);
            *tot_cnt += *cube.cnt;
        }
        else comp_manager->remove_cache_pollutions_of_if_exists(decisions.top());
    }
    reactivate_comps_and_backtrack_trail(false);
    bool ret = propagate(true);
    assert(ret);
    decisions.pop_back();
}
```

- Both branches (`i=0`, `i=1`) are checked at every decision level.
- Every non-null, non-zero branch calls `compute_cube`.
- Every successful `compute_cube` result is added to `mini_cubes`.
- `compute_cube` returning `false` means the SAT solver found the branch infeasible (cache
  pollution). `remove_cache_pollutions_of_if_exists` cleans the bad cache entries. No real
  model is lost.

There is **no scenario** where a cube is computed, counted, but not pushed to `mini_cubes`.

---

## The actual bug: wrong `c.cnt` for non-weighted counting

**Location**: `counter.cpp:1541–1546` (introduced in commit `fba8fa4`):

```cpp
} else {
    // For non-weighted: the cube pins all indep-support vars via the SAT model,
    // so this specific assignment has count = 1.
    // (c.cnt was set to a partial DPLL sub-component count, which is wrong here.)
    c.cnt = fg->one();
}
```

The comment is **incorrect**: the cube does **not** always pin all indep-support variables.

### When the cube is a partial assignment

`compute_cube` builds `c.cnf` from two sources:

1. **Trail decisions** (lines 1473–1486): variables decided at levels `0..D`.
2. **Remaining component variables** (lines 1497–1515): vars in the range
   `[remaining_comps_ofs_[i], unprocessed_comps_end_[i])` at each level `i`.

If at level `D`, sub-component `C1` was **already fully counted and backtracked**:
- `next_unproc_comp()` was called (line 1927), decrementing `unprocessed_comps_end_[D]`.
- `C1`'s variables were `unset_lit`'d during backtracking → not on the trail.
- `C1`'s component index is now `>= unprocessed_comps_end_[D]` → not in the remaining-comps
  range.
- **`C1`'s variables are not in `c.cnf`.**

So when level `D` has: some sub-components `C1..C(n-1)` fully counted
(contributing `K1 × ... × K(n-1)` to `branch_mc[act_branch]`), and `Cn` still being explored
at a deeper level, the cube covers `{decisions, Cn vars}` but **not `C1..C(n-1)` vars**.

### The consequence: models permanently lost

- The **blocking clause** (`ganak_to_cms_cl(c.cnf)` added to the SAT solver) is violated when
  all of `{decisions, Cn=specific}` hold, regardless of `C1..C(n-1)`. This bans
  **`K1 × ... × K(n-1)` models**.
- But `c.cnt = 1`, so only **1 model is counted**.
- The remaining `K1 × ... × K(n-1) - 1` models are **permanently banned but never counted**.
  Future restarts cannot find them because the blocking clause already excludes them.

This is a **soundness bug** (undercounting), not just a performance issue.

### Trigger condition

1. Multiple independent sub-components at some decision level `D`.
2. At least one (`C1`) was fully processed before restart, with `K1 > 1` models.
3. Another sub-component (`Cn`) is still being explored at a deeper level when restart fires.

### Why `CHECK_COUNT_DO` doesn't catch it in typical tests

The assertion at lines 1411–1419 compares `check_count_norestart_cms(cube)` (which returns
`K1`) against `cube.cnt` (which is `1`). However:
- `check_count_norestart_cms` is very expensive (enumerates models via SAT).
- The check is guarded by `fg->exact()`.
- If test formulas decompose into single-variable components or restarts only trigger when
  `branch_mc = 1` everywhere, the assertion never fires.

---

## The correct count

Before the override, `c.cnt` is computed at lines 1460–1467:

```cpp
c.cnt = decisions.top().get_model_side(side)->dup();   // branch_mc[side] at level D
for(int32_t i = 0; i < dec_level(); i++) {
    const auto& mul = decisions[i].get_branch_sols();  // branch_mc[act] at ancestor i
    if (mul == nullptr || mul->is_zero()) continue;
    *c.cnt *= *mul;
}
```

This is `branch_mc[side] × (product of ancestor branch_mc values)`. Each factor accounts for
the already-fully-counted sub-components at that level whose variables are absent from the cube.

**This value is correct for non-weighted counting.** `branch_mc` stores exact integer model
counts when `is_indep=true`. For non-indep components, `include_solution` already forces
`branch_mc = fg->one()`, so their contribution to `c.cnt` is `1` (satisfiability only),
which is also correct.

The comment "DPLL sub-component count is wrong here" applies **only to weighted counting**,
where decision literal weights are applied by `unset_lit` during backtracking — which hasn't
happened yet at restart time. It does **not** apply to non-weighted counting.

---

## Fix

**Remove the non-weighted `c.cnt = fg->one()` override** (lines 1541–1546). Keep `c.cnt` as
computed at lines 1460–1467.

```cpp
// BEFORE (buggy):
} else {
    // For non-weighted: the cube pins all indep-support vars via the SAT model,
    // so this specific assignment has count = 1.
    // (c.cnt was set to a partial DPLL sub-component count, which is wrong here.)
    c.cnt = fg->one();
}

// AFTER (fixed):
// Delete the else block entirely. c.cnt already holds the correct value.
```

The weighted branch is **unaffected** — replacing c.cnt with the SAT model weight product
remains correct for weighted counting.

### Edge case: `!do_use_sat_solver` mode

For the EXIT state (line 1105–1107), `total_model_count()` is recomputed via
`check_count_norestart_cms` when `!do_use_sat_solver`, because non-indep components only track
satisfiability (0/1) and `total_model_count()` can be stale. In `compute_cube`, `c.cnt` is
built from `branch_mc` where non-indep levels already contribute exactly `1` (via
`include_solution` calling `branch_mc[act] = fg->one()` when `!is_indep`). The accumulated
`c.cnt` is therefore correct in both modes without needing a recompute.
