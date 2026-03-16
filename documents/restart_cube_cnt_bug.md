# Bug: Incorrect `c.cnt` in `compute_cube` for Non-Weighted Counting

## Summary

`compute_cube` computed a wrong `c.cnt` value for non-weighted counting modes
(`--mode 0`). Depending on which version of the code was in place, the count
was either heavily overcounted (ancestor-multiplication path) or undercounted
(the `c.cnt = fg->one()` path). Both caused incorrect final model counts when
restarts were enabled.

---

## Background: How Restart Cube Collection Works

At restart time (`restart_if_needed`), Ganak walks the decision stack from the
deepest level back to level 1. For each level and each branch side (0 = left,
1 = right), it calls `compute_cube(cube, side)` to capture the subspace
already counted on that branch. Each cube is then:

1. Added to `mini_cubes` with its associated count `cube.cnt`.
2. Converted into a blocking clause (negation of the cube's CNF) that prevents
   those models from being counted again in the next restart.

The invariant is: the blocking clause bans exactly `cube.cnt` models, and
`cube.cnt` is added to the running total. If `cube.cnt` is wrong, the final
count is wrong.

`compute_cube` builds `cube.cnf` from:
- **Decision literals** (negated, so that `deal_with_irred_cls` restores the
  original values when unit-propagating).
- **Component variables** from the SAT model — for all unprocessed component
  ranges at every decision level, indep-support vars are pinned to their SAT
  model values.
- A deduplication pass to remove vars appearing in multiple component ranges.

`check_count_norestart_cms(cube)` independently counts the models consistent
with the cube by adding the formula's irreducible clauses plus a unit clause
for each literal in `cube.cnf` (negated to force the original value), then
enumerating satisfying assignments via a SAT solver loop.

---

## The Bug

### First (undercounting) form — `c.cnt = fg->one()`

The original code for non-weighted counting was:

```cpp
} else {
    // non-weighted
    c.cnt = fg->one();
}
```

This assumed that the cube uniquely pins **all** indep-support variables,
leaving exactly 1 satisfying assignment. That assumption is wrong whenever
there are indep vars that are neither decided nor covered by any unprocessed
component. Such vars remain free in the cube, and the actual count can be
greater than 1.

**Example from the failing test** (`fuzzTest_2560.cnf`,
`--mode 0 --td 0 --rstfirst 3 --restart 1 --puura 0`):

```
indep_support_end: 17  (indep vars = 1..16)
Cube CNF: 20 -8 -9 -4 10 -1 -2 -3 5 -11 -12 -13 -14 -15 -16
```

Indep vars pinned: {1,2,3,4,5,8,9,10,11,12,13,14,15,16} — 14 out of 16.
Free indep vars: {6, 7}.
Actual model count under the cube: **3** (not 1).

The blocking clause bans all 3 models, but the code credited only 1. The
remaining 2 were permanently lost.

### Second (overcounting) form — ancestor multiplication

Attempting to fix the undercounting, the code was changed to multiply the
top-level `branch_mc[side]` by ancestor `branch_mc` values:

```cpp
c.cnt = decisions.top().get_model_side(side)->dup();
for (int32_t i = 0; i < dec_level(); i++) {
    const auto& mul = decisions[i].get_branch_sols();
    if (mul == nullptr || mul->is_zero()) continue;
    *c.cnt *= *mul;
}
```

This produced `96 × 3 × 4 × 2 = 2304`, while the actual count was 3. The
`branch_mc` values are **aggregate** counts accumulated across many different
variable assignments in the DPLL tree. They cannot be used directly because:

- The cube pins ONE specific partial assignment (from the current SAT model).
- The `branch_mc` aggregate reflects many different assignments, most of which
  are inconsistent with the specific assignment pinned by the cube.

---

## Root Cause

The analytical formula for `c.cnt` is non-trivial for non-weighted counting.
The cube leaves some indep-support vars unpinned, and the count of satisfying
assignments for those free vars under the formula is not captured by the DPLL
`branch_mc` accumulators.

For **weighted** counting, `c.cnt` can be computed analytically as the product
of all literal weights in the SAT model (one specific assignment, so one
specific weight product). This was already correct.

For **non-weighted** counting, the correct count is the number of complete
indep-support assignments extending the cube's partial assignment that satisfy
the formula. This requires actually enumerating models — exactly what
`check_count_norestart_cms` does.

---

## Fix

After building `cube.cnf` (decisions + component vars + dedup), replace the
`c.cnt` computation for non-weighted with a call to `check_count_norestart_cms`:

```cpp
if (weighted()) {
    // Weighted: product of all SAT-model literal weights.
    c.cnt = fg->one();
    for (uint32_t v = 1; v < opt_indep_support_end; v++) {
        Lit l(v, sat_solver->get_model()[v-1] == CMSat::l_True);
        *c.cnt *= *get_weight(l);
    }
    if (c.cnt->is_zero()) return false;
} else {
    // Non-weighted: enumerate models consistent with the cube constraints.
    c.cnt = check_count_norestart_cms(c);
    if (c.cnt->is_zero()) return false;
}
```

This is correct because `check_count_norestart_cms` adds all formula
constraints and forces the cube's pinned vars, then counts satisfying
indep-support assignments. The result matches the number of models that the
blocking clause (derived from this cube) will exclude.

The overhead is acceptable: restarts are infrequent, and the enumeration is
bounded by the model count, which is typically small at individual cube
granularity.

---

## Steps to Find the Bug

1. **Initial report**: performance bug — some counted subspaces were not being
   collected as cubes, so the same subspaces were re-explored after each
   restart.

2. **Code review** of `restart_if_needed` and `compute_cube`:
   - Found that `c.cnt = fg->one()` was applied unconditionally for
     non-weighted, regardless of how many models the cube's partial assignment
     actually covers.
   - Wrote analysis to `documents/restart_cube_analysis.md`.

3. **First fix attempt**: removed the `c.cnt = fg->one()` override and kept
   the ancestor `branch_mc` multiplication, reasoning that processed
   sub-components (not pinned in the cube) contribute their aggregate counts.

4. **Test failure**: running the failing case triggered the `CHECK_COUNT`
   assertion inside `compute_cube`:
   ```
   ERROR [compute_cube]: cnt mismatch for cube: CNF: 20 -8 -9 -4 10 -1 -2 -3 5 -11 -12 -13 -14 -15 -16
     recorded c.cnt : 2304
     actual check   : 3
   ```
   Debug output showed:
   - All "excluded comp" indep vars were covered by the cube (0 uncovered).
   - The ancestor branch_mc product (96 × 3 × 4 × 2 = 2304) was far too large.
   - The actual count (3) matched the free indep vars {6,7} with 3 satisfying
     assignments under the formula.

5. **Root cause identified**: `branch_mc` values are aggregate DPLL counts
   inconsistent with the specific SAT model assignment pinned by the cube.
   There is no simple analytical formula; the count must be computed by
   enumeration.

6. **Correct fix**: use `check_count_norestart_cms(c)` for non-weighted. This
   already existed as a correctness check (via `CHECK_COUNT`) and as a
   fallback path for `!do_use_sat_solver`; it was simply not applied for the
   normal restart path. After the fix, the `CHECK_COUNT` assertion trivially
   passes (the value was just computed the same way), and the model count
   matches the no-restart baseline.
