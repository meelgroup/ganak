# bdnnf: queries on Ganak's d-DNNF circuits

`bdnnf` answers queries on the circuits `ganak --compile` writes, in the BDNNF
format ([spec-bdnnf.md](spec-bdnnf.md)). Its query set follows
[ddnnife](https://github.com/practical-kc/ddnnife) by Sundermann, Raab, Heß,
Thüm and Schaefer [[1]](#references); see the [comparison](#comparison-with-ddnnife)
below.

```shell
./ganak --compile out.bdnnf in.cnf
./bdnnf count out.bdnnf
./bdnnf sample out.bdnnf -n 100 --cnf in.cnf --assume 3,-5
```

## Usage

`bdnnf <command> <in.bdnnf> [options]`. Common options:

- `--cnf in.cnf`: literal weights (`c p weight`) for weighted queries. The
  projection set is stored in the circuit; a CNF with a different one is refused.
  Commands that print models also check each one against the CNF's clauses.
- `--assume L1,L2,...`: condition the query on literals of projection variables.

Models are printed as DIMACS lines over the projection variables (`--full`:
all variables). The tools stream the file, keeping a few values per node.

| Command | What it does |
|---|---|
| `count` | (weighted, projected) model count |
| `sample [-n N] [--seed S] [-o out]` | random models, by weight (uniform when unweighted) |
| `enum [-n N] [--start I] [-o out]` | models in a fixed order; `--start` resumes from model I |
| `lits` | for each literal, the count of the models containing it |
| `core` | literals true in every model (core; negated: dead) |
| `batch FILE [--sat]` | one count (or SAT/UNSAT) per line of literals in FILE |
| `atomic [--cross]` | sets of variables equal (with `--cross`: or opposite) in every model |
| `opt --obj FILE [-k K] [--min]` | the k best models under a sum of per-literal values (`literal value` lines) |
| `twise [-t T] [--seed S] [-o out]` | samples covering every satisfiable t-tuple of literals |
| `stream` | the queries above, one per stdin line (`count --assume 3`, ...), circuit held in memory; each answer ends with `c done` |
| `verify` | check the format and every property the file claims |
| `dump`, `dot`, `to-d4` | text form, Graphviz rendering, d4 `.nnf` |

## Comparison with ddnnife

ddnnife [[1]](#references) is a d-DNNF reasoner for feature-model analysis; the
queries of `bdnnf` were modelled on it. The table maps ddnnife's CLI and stream
operations to `bdnnf`.

| ddnnife | bdnnf | Notes |
|---|---|---|
| `count -a A` | `count --assume A` | |
| `count` with a list of literals | `lits --assume A` | `lits` gives the count of every literal |
| `count-features` | `lits` | |
| `count-queries FILE` | `batch FILE` | ddnnife can use several threads; bdnnf evaluates up to 64 queries per pass |
| `sat FILE` | `batch FILE --sat` | |
| `core` | `core` | |
| `atomic-sets [--cross]` | `atomic [--cross]` | ddnnife can restrict the candidates; bdnnf cannot yet |
| `anomalies` | `core`, `atomic` | |
| `urs -n N -s S` | `sample -n N --seed S` | |
| stream `enum` | `enum -n N --start I` | ddnnife continues across stream calls; bdnnf resumes by index |
| `t-wise -t T` | `twise -t T` | different constructions: ddnnife merges partial samples bottom-up, bdnnf covers tuples greedily (as in YASA); ddnnife can restrict the covered literals or variables |
| `t-wise-check` | — | |
| stream fitness / top-k configurations | `opt --obj FILE -k K` | |
| `stream`, `stream-queries FILE` | `stream` (and `stream < FILE`) | different line syntax; ddnnife recounts incrementally (marking), bdnnf recomputes each query |
| stream `v` (per-variable answers) | partly | `lits` and `core` answer for all variables |
| `mermaid` | `dot` | |
| `statistics` | — | |
| `to-cnf` | — | |
| `--save-ddnnf` (c2d) | `to-d4` (d4) | |
| reads c2d and d4 `.nnf`, and CNF via d4 | reads BDNNF from `ganak --compile` | |
| Rust API, Python and Kotlin bindings | — | |

Differences in scope: Ganak's circuits can be weighted and projected, so
`bdnnf` takes literal weights and a projection set, and completes models over
the variables Ganak's preprocessing defined away (spec, section 1.1). ddnnife
works on smooth circuits and smooths d4 input when loading it; `bdnnf` smooths
while evaluating.

As a cross-check, `bdnnf lits` on ddnnife's `example_input/auto1.cnf`
(compiled by Ganak) gives, for all 2513 features, the counts in ddnnife's
reference `tests/data/auto1_sb_fs.csv`. count_fuzzer's `ddnnife_diff.py`
compares the two tools on random unweighted, unprojected CNFs, through
`bdnnf to-d4`: counts, conditioned counts, query files, SAT queries, core
literals, feature counts, atomic sets and enumeration must agree, and both
tools' t-wise samples must cover every satisfiable pair.

## References

[1] Chico Sundermann, Heiko Raab, Tobias Heß, Thomas Thüm and Ina Schaefer.
*Reusing d-DNNFs for Efficient Feature-Model Counting.* ACM Transactions on
Software Engineering and Methodology (TOSEM), 2024.
Preprint: [arXiv:2303.12383](https://arxiv.org/abs/2303.12383).
Tool: [github.com/practical-kc/ddnnife](https://github.com/practical-kc/ddnnife).
