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
