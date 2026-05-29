# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Ganak is a high-performance exact model counter. It counts the number of
satisfying assignments of a CNF formula, supporting multiple algebraic fields
(integers, rationals, complex, modular arithmetic, polynomials over finite
fields).

## Notes about the code
* Never use `check_count_norestart_cms` or `check_count_norestart` in normal code,
  only in `CHECK_COUNT_DO` (i.e. `CHECK_COUNT` defined) code

## Building
Always disable ccache and build the `ganak` target:

```
CCACHE_DISABLE=1 cmake --build build --target ganak -j8 2>&1 | tail -20
```

The build directory is `build/`. Go to that directory, and run:

```
./build_norm.sh
```

Common build options:
- `-DCMAKE_BUILD_TYPE=Debug` — debug symbols, assertions enabled
- `-DSANITIZE=ON` — enable Clang sanitizers
- `-DBUILD_SHARED_LIBS=OFF` — static binary

## Fuzzing (DO THIS AFTER EVERY CHANGE)

Ganak MUST be fuzzed after every change. A wrong count is a silent, severe bug,
and the fuzzers are the main line of defence. Build the binary first — note the
*executable* target is `ganak-bin` (output `build/ganak`); the `ganak` target
only builds the library:

```
CCACHE_DISABLE=1 cmake --build build --target ganak-bin -j4 2>&1 | tail -20
```

### 1. count_fuzzer (`../count_fuzzer`) — the primary fuzzer

Lives at `../count_fuzzer` (a sibling of this repo). Git remote:
`git@github.com:meelgroup/count_fuzzer.git` (`meelgroup/count_fuzzer`). It is a
crude, Ganak/ApproxMC-specific fork of
[SharpVelvet](https://github.com/meelgroup/SharpVelvet).

How it works: it builds random CNF instances with the generators in
`generators/` (build them once: `cd generators && g++ cnf-fuzz-biere.c -o
biere-fuzz`), runs `../ganak/build/ganak` on each with a wide variety of option
combinations (cache on/off, restarts, td, arjun, threads, polarity, modes,
weights, projection, complex, ...), and checks the result. For exact counting it
cross-checks the count against a trusted oracle / alternate configs; it also
catches crashes, assertion failures, OOM, and timeouts.

```
cd ../count_fuzzer
./fuzz.py --only 20 --exact      # quick sanity check after a build
./fuzz.py --only 200 --exact     # thorough
./fuzz.py --unweighted           # only unweighted
./fuzz.py --weighted             # only weighted
./fuzz.py --proj | --unproj      # projected / unprojected only
./fuzz.py --cpx                  # complex field only
./fuzz.py --threads K            # fuzz with --threads K passed to ganak
./fuzz.py --tout T               # per-instance timeout (default 4s)
./fuzz.py                        # non-stop: EVERYTHING (cpx, proj, weighted, ...)
```

`fuzz_session_ganak.sh` runs a longer pre-canned session.

### 2. d-DNNF compilation fuzzer (`tests/ddnnf_fuzz.py`)

Validates the `--compile` (d-DNNF) feature specifically. It generates random
CNFs, brute-forces the true model count + model set as an oracle, and checks the
emitted d4 `.nnf` circuit (parsed/counted by `tests/ddnnf_verify.py`):

```
python3 tests/ddnnf_fuzz.py 200            # circuit count == true count, model set
                                           #   == true models (validates arc literals),
                                           #   strictly decomposable, faithful as a fn
python3 tests/ddnnf_fuzz.py 200 --seed N   # reproduce a specific run
python3 tests/ddnnf_fuzz.py 200 --maxvars 14  # cap var count (oracle is 2^nv, so
                                           #   this is the main runtime/coverage knob;
                                           #   also --minvars N)
```

`ganak --compile out.nnf in.cnf` writes a faithful d-DNNF. `--weak 0` is that
(the default); `--weak 3` is the synthesis share-and-branch (see section 3 and
`--help`). `--ddnfcheck 1` makes Ganak cross-check each decision
level's circuit sub-count against its own count (debugging the compiler). Temp
files are reserved race-proof (atomic `O_CREAT|O_EXCL`) so multiple fuzzer
processes can run concurrently. Failing cases are copied to
`/tmp/ddnnf_fuzz/fail_*.{cnf,nnf}`.

### 3. Functional-synthesis round-trip (`tests/ddnnf_synth.py`)

For a *projected* CNF (inputs X = sampling vars, outputs Y = the rest), the
compiled circuit can be used for Boolean functional synthesis: for an input
assignment X, `ddnnf_verify.synthesize()` reads a witness ψ(X) off the circuit.
The fuzzer checks that for every satisfiable X, `F(X, ψ(X))` holds, and reports
the circuit size vs `--weak 0`.

```
python3 tests/ddnnf_synth.py 200            # --weak 0 (faithful d-DNNF)
python3 tests/ddnnf_synth.py 200 --weak 3   # share-and-branch; also prints size ratio
```

Synthesis notes (empirically established):
- The witness for Y comes from the **SAT oracle** (it stays on in `--compile`),
  which records one example assignment of Y at each SAT leaf (`set_override` in
  `DDNNFCompiler`). This is the main compactness for synthesis (Y is not
  enumerated). `--weak 0` + SAT is a correct, compact synthesis compiler.
- `--weak 3` (share-and-branch): input vars (`< opt_indep_support_end`) may be
  shared across AND children while outputs stay disjoint; cache is selective
  (shared components are not cached, for soundness). It is **sound** (synthesis
  round-trip 0 failures) but **not more compact** than `--weak 0` (~1.05-1.10x
  larger): duplicating shared inputs and losing their caching outweighs the
  finer decomposition. (A weak cut can be smaller only by dropping constraints /
  over-approximating, which is unsound for synthesis -- so `--weak 3` can't
  reproduce that and stays faithful instead.)

## Running Tests
```
cd build && ctest -V
```

Tests use LLVM-lit and are in `tests/cnf-files/`. They check model counts
against known correct values.

To run a single test file:
```
cd build && lit tests/cnf-files/a.cnf
```

### d-DNNF compile tests

`ctest` also runs deterministic d-DNNF compile/cleanup checks. They assert only
*semantic* invariants of the emitted circuit (model count, model set,
decomposability, canonical cleaned form) — never the byte-exact `.nnf` — so they
do not break when search heuristics reshape the circuit:
- `tests/cnf-files/ddnnf/*.cnf` — lit fixtures: `--compile`, then
  `ddnnf_verify.py` (count/model-set/decomposable), `ddnnf-cleanup`, re-verify
  the cleaned file `--strict`, and `ddnnf2dot`. `ddnnf_verify.py` doubles as a
  CLI checker (`--expect-count`, `--cnf`, `--check-decomposable`,
  `--check-weak-decomposable`, `--strict`); bare invocation still prints the count.
- ctest targets `ddnnf_compile_fuzz`, `ddnnf_synth_weak0`, `ddnnf_synth_weak3` —
  the property fuzzers run seeded (so a regression reproduces) and self-check
  against a brute-force oracle.

## Architecture

The counting pipeline flows through these layers:

1. **`main.cpp`** — CLI parsing, reads CNF via Arjun's `SimplifiedCNF`, runs
   Arjun simplification, then calls `setup_ganak()` + `Ganak::count()`
2. **`ganak.hpp/.cpp`** — Public API (pimpl idiom via `CDat`). Delegates to
   `OuterCounter`.
3. **`outer_counter.hpp/.cpp`** — Orchestrates single vs. multi-threaded
   counting. Creates one `Counter` per thread.
4. **`counter.hpp/.cpp`** — Core DPLL-based counting engine. Implements
   DPLL-based branching, unit propagation, clause learning, restarts, and
   component-based decomposition. The main loop is `outer_count()`.
5. **`comp_management.hpp/.cpp`** — Manages connected component decomposition.
   Contains `CompManager` (exactly one per `Counter`) which uses `CompAnalyzer`
   to detect independent sub-problems and `CompCache` to memoize sub-counts.
6. **`comp_cache.hpp`** — Hash-table cache mapping component fingerprints to
   their counts. Template on cache entry type. Controlled by `--cache 0` to
   disable.
7. **`comp_analyzer.hpp/.cpp`** — Detects connected components in the current
   partial assignment using variable/clause stamps.
8. **`comp_types/`** — Component representation types: `Comp` (live component),
   `CacheableComp` (stored in cache), `HashedComp`, `DifferencePackedComp`.
9. **`counter_config.hpp`** — `CounterConfiguration` struct: all tunable
   parameters (restarts, cache size, LBD cutoffs, vivification, probabilistic
   hashing, tree decomposition, etc.).
10. **`statistics.hpp/.cpp`** — `DataAndStatistics`: runtime stats collected
    throughout solving.
11. **`structures.hpp`** — Core data structures: `LitWatchList`,
    `DecisionLevel`, `TriValue` (true/false/undef), `LiteralIndexedVector`.

### Field system
All arithmetic is polymorphic over `CMSat::Field` / `CMSat::FieldGen` (from
cryptominisat5). `FG = unique_ptr<FieldGen>`, `FF = unique_ptr<Field>`. Field
implementations are in:
- `mcomplex.hpp` / `mcomplex-mpfr.hpp` — complex rationals / complex floats
- `mparity.hpp` — parity (mod 2) counting
- `mpoly.hpp` — polynomials over finite fields

### Counting modes (`--mode N`)
| Mode | Field | Notes |
|------|-------|-------|
| 0 | Integer (default) | Supports `--appmct` approximate fallback |
| 1 | Rational (GMP mpq) | Exact weighted counting |
| 2 | Complex rational | Weights specified as `a b` (real, imag) |
| 3 | Polynomial over finite field | Requires `--npolyvars N` |
| 4 | Parity (mod 2) | |
| 5 | Integer mod prime | Requires `--prime X` |
| 6 | Complex float (MPFR) | |

## Debug flags (in `src/common.hpp`)

These are compile-time `#define`s near the top of `src/common.hpp`, all
commented out by default. Uncomment the ones you need, then **rebuild**
(`CCACHE_DISABLE=1 cmake --build build --target ganak-bin -j4`). Each flag gates
a `*_DO(x)` macro (e.g. `SLOW_DEBUG_DO(...)`) that expands to the wrapped check
only when the flag is on, so they cost nothing in normal builds but get *very*
slow when enabled. Turn on the least expensive flag that catches your bug.

| Flag | Gates | What it does | Cost |
|------|-------|--------------|------|
| `VERBOSE_DEBUG` | `debug_print` / `VERBOSE_DEBUG_DO` | Prints a full trace of the search (decisions, propagations, backtracks, component splits, counts). Enable on a *minimized* instance only — output is huge. | very high |
| `SLOW_DEBUG` | `SLOW_DEBUG_DO` | Internal consistency assertions (trail sanity, component/var invariants, cache invariants). First thing to turn on for a suspected logic bug. | high |
| `CHECK_PROPAGATED` | `CHECK_PROPAGATED_DO` | After key steps, asserts every clause is properly propagated / not silently conflicting (`check_all_propagated_conflicted`, `check_trail`). Catches BCP / watch-list bugs. | high |
| `CHECK_IMPLIED` | `CHECK_IMPLIED_DO` | Asserts learnt / 1-UIP clauses are actually implied by the formula (`check_implied`). Catches conflict-analysis bugs. | high |
| `VERY_SLOW_DEBUG` | `VERY_SLOW_DEBUG_DO` | Even more expensive invariant checks than `SLOW_DEBUG`. | very high |
| `ANALYZE_VERBOSE` | `analyze_verb` | Verbose tracing inside component analysis (`comp_analyzer.cpp`) only. | high |
| `CHECK_TRAIL_ENTAILMENT` | trail-entailment check | The slowest check that does *not* verify counts: re-verifies the whole trail is entailed at each step. Last resort for trail/entailment bugs. | extreme |
| `CHECK_COUNT` | `CHECK_COUNT_DO` | **The flag for wrong-count bugs.** Cross-checks every (sub-)count against an independent recount (`check_count*`). Also disables count reuse so cross-checking is valid — additionally pass `--cache 0` on the command line (see the note next to the define). Only `CHECK_COUNT_DO`/`check_count_norestart*` code may call those checkers (see "Notes about the code"). | extreme |
| `BUDDY_ENABLED` | BuDDy code | Not a debug check: compile-time enable of BuDDy BDD-based counting of small components (pairs with the `--buddy` option). | n/a |

## Debugging a fuzzer-found bug

When `count_fuzzer` reports a wrong count or crash, isolate it before staring at
code. `../count_fuzzer` is a fairly complete find-and-isolate system:

1. **Minimize.** `minim_all.py` is the one-stop pipeline — it auto-detects
   crash vs. wrong-count, backs up the file, then minimizes options, weights,
   and clauses, re-verifying at each step:
   ```
   cd ../count_fuzzer
   ./minim_all.py "../ganak/build/ganak --mode 1 --polar 1 file.cnf"
   ```
   Individual tools if you need them:
   - `minim_cnf.py` — **delta debugger** for clauses (hierarchical delta
     debugging); produces `file_min_clauses.cnf`.
   - `minim_opts.py` / `minim_opts_crash.py` / `minim_opts_count.py` — minimize
     the command-line options (auto / crash-preserving / count-preserving;
     count-preserving allows 0.1% tolerance).
   - `minim_weights.py` — minimize `c p weight` lines.
   - `propagate.py` — unit-propagate a CNF to fixpoint to shrink it before
     minimizing.

2. **Pinpoint.** Once you have a tiny instance, rebuild with the relevant
   `common.hpp` flags and run it directly:
   - wrong count → enable `CHECK_COUNT` **and** run with `--cache 0`; it aborts
     at the first component whose recount disagrees.
   - crash / bad invariant → enable `SLOW_DEBUG` (then `CHECK_PROPAGATED` /
     `CHECK_IMPLIED` / `CHECK_TRAIL_ENTAILMENT` as needed).
   - then add `VERBOSE_DEBUG` to read the exact trace around the failure.

3. **Re-fuzz** after the fix (`./fuzz.py --only 200 ...`) before committing.

## Data Analysis

Previous benchmark runs are stored under `build/data/` as `out-ganak-*/`
directories. The relevant ones have already been parsed with
`./get_data_ganak.py` into `build/data/data.sqlite3`.

To view statistics about the data:
```
cd build/data && ./create_graphs_ganak.py --nograph
```

To query the database directly:
```
cd build/data && sqlite3 data.sqlite3
```

Then use SQL, e.g.:
```sql
SELECT dirname, count(*), avg(ganak_time), avg(cache_miss_rate) FROM data GROUP BY dirname;
```

### SQLite schema (`data` table)

| Column | Description |
|--------|-------------|
| `solver` | Solver name (ganak, d4, gpmc, ...) |
| `dirname` | Run directory (identifies the benchmark batch + config) |
| `fname` | CNF filename |
| `ganak_time` | Wall time in seconds (NULL = timeout/unsolved) |
| `ganak_mem_MB` | Peak memory in MB |
| `ganak_call` | Command-line flags used |
| `ganak_ver` | `ganak-<sha>-<cms-sha>` version string |
| `signal` | Exit signal (11=SEGV, 6=ABRT, 8=FPE, 14=ALRM) |
| `mem_out` | 1 if OOM (bad_alloc) |
| `timeout_t` | Timeout wrapper's reported time |
| `conflicts` | Total conflicts during solving |
| `decisionsK` | Decisions in thousands |
| `compsK` | Cache lookups/stores in thousands |
| `cache_miss_rate` | Component cache miss rate |
| `cache_del_time` | Time spent on cache deletion |
| `cache_avg_hit_vars` / `cache_avg_store_vars` | Avg vars in cached components |
| `indep_sz` | Sampling set size after Arjun |
| `opt_indep_sz` | Optional sampling set size |
| `orig_proj_sz` | Original projection set size |
| `new_nvars` | Vars after Arjun simplification |
| `unkn_sz` | Unknown vars at Arjun start |
| `arjun_time` | Time spent in Arjun preprocessing |
| `backbone_time` | Backbone computation time |
| `td_width` | Tree decomposition width |
| `td_time` | Tree decomposition time |
| `restarts` | Number of restarts |
| `cubes_orig` / `cubes_final` | Cubes before/after filtering per restart |
| `sat_called` | Number of SAT oracle calls |
| `gates_extended` / `gates_extend_t` | Gates added by extension + time |
| `padoa_extended` / `padoa_extend_t` | Vars added by Padoa extension + time |
| `primal_density` / `primal_edge_var_ratio` | Primal graph density metrics |

## Dependencies

GMP, MPFR, FLINT, cryptominisat5, arjun, approxmc, treedecomp, zlib (optional). ALL dependencies
are under ../ -- you can find e.g. ../arjun ../appproxmc ../treedecomp ../sbva ../cadical ../cryptominisat5 etc
