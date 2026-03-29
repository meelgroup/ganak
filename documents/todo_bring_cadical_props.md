# Plan: Importing CaDiCaL Propagation & Conflict Analysis Ideas into Ganak

## Summary of Findings

**Already in Ganak (don't port):**
- Blocking literals (`blckLit`) — identical to CaDiCaL's `blit`
- Binary clause separation — Ganak uses `binaries` vec, arguably better than CaDiCaL's `size==2` check
- In-place watch list compaction with two-pointer iteration (`it` / `it2`)
- 1-UIP conflict analysis
- Recursive CCMin + binary clause minimization
- `__builtin_prefetch` hints in propagation

**Missing from Ganak (candidate improvements):**

---

## Change 1: Watch Position Caching (Ian Gent 2013)

**What CaDiCaL does:** Stores `clause->pos` — the index where the last replacement watch was found. The next propagation for that clause starts the scan at `pos`, not at 2, then wraps. This turns repeated quadratic re-scanning into amortized linear.

**What Ganak does:** Always scans from index 2 in `counter.cpp`'s propagation loop:
```cpp
uint32_t i = 2;
for(; i < c.sz; i++) if (!is_false(c[i])) break;
```

**What to change:**
1. Add `uint32_t pos = 2` to the `Clause` struct in `structures.hpp`
2. In the propagation loop in `counter.cpp`, change the watch replacement search to start from `c.pos`, scan to end, then wrap from 2 to `c.pos`, update `c.pos` on success

**Files:** `src/structures.hpp`, `src/counter.cpp` (~15 lines changed)
**Risk:** Low. Well-understood optimization with no algorithmic changes.
**Expected gain:** Measurable propagation speedup, especially on long clauses.

---

## Change 2: EVSIDS — Exponential Variable Activity Decay

**What CaDiCaL does:** Uses a growing `score_inc`. Each conflict, `score_inc *= multiplier` (e.g. 1.0/0.95 ≈ 1.053). Bumping a variable adds `score_inc`, not 1.0. Periodic rescaling when any score exceeds ~1e150. This achieves an implicit exponential decay of all older scores relative to newer ones without touching every variable.

**What Ganak does:** `inc_act()` adds 1.0 to `watches[lit].activity`. Decay happens only via `vsads_readjust()` which halves **all** activities every N decisions (`*= 0.5`). This means rescaling all N watches every N decisions — O(|vars|) work periodically.

**What to change:**
1. Add `double act_inc = 1.0` to `Counter`
2. Change `inc_act()` to `watches[lit].activity += act_inc`
3. After each conflict (in `resolve_conflict()` or where `inc_act` is called), multiply `act_inc *= (1.0 / conf.act_decay)` (e.g., decay = 0.95)
4. When any activity exceeds `1e100`, rescale all activities and `act_inc` by `1e-100`
5. Remove or adapt `vsads_readjust()` (halving every N decisions becomes redundant)
6. Add `double act_decay = 0.95` to `CounterConfiguration`

**Files:** `src/counter.hpp`, `src/counter.cpp`, `src/counter_config.hpp` (~30 lines changed)
**Risk:** Low-medium. Well-studied, but the interaction with the `order_heap` and the SAT-mode VSIDS needs care.
**Expected gain:** Better variable ordering → fewer conflicts on hard instances.

---

## Change 3: Three-Tier Learned Clause Management

**What CaDiCaL does:** Strict 3-tier system based on glue (LBD):
- **Tier 1** (lbd ≤ 2): Never deleted
- **Tier 2** (lbd 3–6): Delete only if `used == 0` since last reduce (reset `used` flag each reduce phase)
- **Tier 3** (lbd > 6): Delete unless `used` since last reduce *and* kept under budget

**What Ganak does:** One adaptive threshold (`lbd_cutoff`, starts at 2, can grow). Below threshold = protected. Above = sorted by `total_used`, bottom half deleted. There is no "tier 2" concept: a clause with lbd=4 is treated identically to one with lbd=10.

**What to change:**
1. In `reduce_db()` in `counter.cpp`:
   - Protect lbd ≤ 2 unconditionally (already done via `lbd_cutoff`)
   - Add: keep lbd 3–6 if `cl.used == 1` (reset `used` to 0 after checking)
   - Delete lbd > 6 if `cl.used == 0`, even if `total_used` is high
2. Add a `lbd_tier2_cutoff = 6` to `CounterConfiguration`
3. Reset `used` flag on all surviving clauses after each `reduce_db()` call

**Files:** `src/counter.cpp` (~30 lines changed), `src/counter_config.hpp` (~2 lines)
**Risk:** Low-medium. Changes clause retention policy; may need tuning of `lbd_tier2_cutoff`.
**Expected gain:** Better quality clause database — retains more "medium-glue" useful clauses while being more aggressive with stale ones.

---

## What I'm NOT Porting (and Why)

| CaDiCaL Feature | Reason to skip |
|---|---|
| Clause `size` in Watch entry | Ganak already separates binary clauses into a distinct `binaries` vec — no benefit |
| On-The-Fly Strengthening (OTFS) | Complex, interacts with Ganak's counting semantics (especially component caching); unclear benefit |
| Chronological backtracking | Ganak already has `chrono_check()` and decision flipping — semantically equivalent |
| LRAT proof tracking | Not relevant to model counting |
| VMTF queue | Ganak uses tree-decomposition-guided branching which supersedes VMTF; mixing them is non-trivial |

---

## Implementation Order

1. **Change 1** (watch position caching) — simplest, zero-risk, can be validated immediately
2. **Change 2** (EVSIDS) — independent of #1, validate with fuzz tests
3. **Change 3** (three-tier clause management) — slightly more disruptive, validate last

Each change will be followed by a fuzz test run (`./fuzz.py --only 20 --exact`) and a full `ctest` to verify correctness before proceeding to the next.
