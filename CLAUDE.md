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

The build directory is `build/`. If it doesn't exist or CMake needs to be re-run:

```
mkdir -p build && cd build && cmake ..
```

Common build options:
- `-DCMAKE_BUILD_TYPE=Debug` — debug symbols, assertions enabled
- `-DSANITIZE=ON` — enable Clang sanitizers
- `-DBUILD_SHARED_LIBS=OFF` — static binary

## After Building: Fuzz Testing
After a new build, run the fuzzer to check for regressions:
```
cd ~/development/sat_solvers/count_fuzzer && ./fuzz.py --only 20 --exact
```

`--only K` limits fuzz runs to K. ~20 is enough for a quick sanity check; run
more (e.g. `--only 200`) for thorough validation.

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
Uncomment to enable expensive debug checks:
- `VERBOSE_DEBUG` — verbose output
- `SLOW_DEBUG` — slow assertion checks
- `CHECK_COUNT` — verify counts (disables cache to allow cross-checking)
- `CHECK_TRAIL_ENTAILMENT` — slowest trail check

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
