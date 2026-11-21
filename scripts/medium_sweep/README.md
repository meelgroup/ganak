# Bounded preprocessing for medium-duration model counts

This experimental Python driver runs a structural prepass on plain unweighted, unprojected DIMACS before calling a sweep-enabled Ganak. Unit propagation removes forced variables. Strongly connected components of binary implications identify equivalent literals. The driver restores the contribution of free variables by multiplying by a power of two. Certified exactly-one blocks can use bounded adaptive domain search. A restricted residual form can use a five-second SAT enumeration probe. A declined probe falls through to Ganak within the same overall deadline.

Build the two helper programs with a C++20 compiler, GMP and CaDiCaL development headers and libraries. Set CPPFLAGS and LDFLAGS if CaDiCaL is outside the compiler’s search paths. Run from this directory.

```sh
make
export GANAK=/absolute/path/to/sweep-enabled/ganak
python3 check_policy.py
python3 policy.py --seconds 600 instance.cnf
```

The completed local pilot compared stock develop 9a05d35, its sweep port bd4cb7e, and the combined driver on 13 formulas with measured stock runtimes of 60 through 600 seconds. At 60 seconds they solved 0, 4 and 6 respectively. At 600 seconds all three solved 13. Eleven cases came from previously studied inputs. Two came from a separate stock screen. The two additional one-minute solves were in the previously studied group. Per-run data and cohort-level coverage are under results/. The combined arm denotes the full policy.py driver. These measurements used one thread, a 9 GB address-space ceiling and a shared local host.

## Larger local benchmark

On Linux, benchmark/large.py screens the frozen manifest for 200 distinct formulas with a stock develop runtime of 60–600 seconds. Each qualifying input is then run with the existing sweep and combined policy on the same physical core. The manifest contains archive paths and scheduling priorities. Provide the original CNFs separately. For 2022 and 2025 the runner reads heur/reference/mc22.tar and mc25.tar under GANAK_DATA. For other years it reads mcYYYY-track1-mc_competition/mcYYYY_track1_NNN.cnf.xz. The submitted instances come from heur/reference/collected.tar.

Set GANAK_STOCK and GANAK to absolute paths to the stock and sweep executables, GANAK_DATA to that archive root, GANAK_OUTPUT to a writable result directory and GANAK_CORES to a comma-separated list of physical CPU IDs. Run `python3 benchmark/large.py` to start or resume and `python3 benchmark/report.py` to export the current CSVs. The frozen manifest reproduces the input pool and priorities. Actual selection can vary with host speed and completion order.

## Running a paired competition replay

The pilot CSVs describe the previous reduction policy. The current default keeps the unit-simplified formula when units alone remove at least a quarter of the variables and the subsequent equivalence pass leaves more than half of the unit residual's variables. This experimental rule targets cases like MC2026-009 where extra rewriting slowed search. Use `--previous-reduction` on the same driver to test the change.

For ordinary Track 1, use one build of this branch, the same one-thread and memory settings and the same wall limit in each arm. Run directly from the directory above after building the helpers. Both plain and gzip-compressed DIMACS are supported. The driver is opt-in; launching the native executable alone bypasses the prepass.

```sh
# Reference search and sweep in the same executable
"$GANAK" --threads 1 --prob 0 --maxcache 5000 --frontierstates 0 --groupedsweep 0 --frontierprofile 1 input.cnf
"$GANAK" --threads 1 --prob 0 --maxcache 5000 --frontierstates 8192 --groupedsweep 1 --frontierminvars 4 --frontierprofile 1 input.cnf

# Same sweep with the previous and current preprocessing selections
python3 policy.py --seconds 3600 --maxcache 5000 --previous-reduction input.cnf.gz
python3 policy.py --seconds 3600 --maxcache 5000 input.cnf.gz
```

Use an external process-group wall timeout on all arms and retain stdout, stderr, exit status and GNU time output per formula. Fix solver options and rotate arm order across inputs. Compare counts and solved formula IDs at each cutoff. Weighted and projected tracks require their native Ganak configurations.

On diagnostic replays add `--frontierbucket 0` to the native sweep or `--no-buckets` to the driver for the Boolean bucket ablation. Word, sparse and grouped counting can still run. Use `--frontierprofile 2` or driver option `--kernel-timing` for CPU seconds spent in each Boolean kernel, including declined calls. Level 2 makes additional clock reads, so use the same profile level for comparisons.

The default native sweep summary includes completions by kernel. Profiled logs also give exact cache lookups, hits, decisions, kernel calls and returns. Subtract hits from lookups for absolute component-cache misses. Each `c frontier sample` is cumulative. Pair begin/end records by process and call ID, and use their difference for one completed count. Snapshots are intermediate readings; a timeout can leave only a partial interval. Hook CPU includes eligible snapshots and other work outside the three timed Boolean kernels. Root completions are logged separately. The driver emits `c medium` stage, decline and route records, forwards native mechanism diagnostics, and emits the final count after restoring free variables. The local benchmark runner stores compressed stdout/stderr in `GANAK_OUTPUT/logs`, including the available prefix on timeout. Compression is outside the measured wall interval.
