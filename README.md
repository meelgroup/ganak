[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![build](https://github.com/meelgroup/ganak/workflows/build/badge.svg)

# Ganak, A High-Performance Versatile Model Counter
Ganak is a state-of-the art model counter that won _all_
available awards at the 2024 and 2025 Model Counting Competitions. It can count over _any_
field, including but not limited to integers, rationals, complex numbers,
integers modulo prime, and polynomials over a finite field.

It can run in probabilistic mode (default), non-probabilistic mode (`--prob 0`),
and approximate counting mode (`--appmct <timeout>`), where it switches to
[ApproxMC](https://github.com/meelgroup/approxmc/) after a predetermined timeout.

To read more about Ganak-specific ideas, refer to [our newer
paper](https://www.msoos.org/wordpress/wp-content/uploads/2025/05/ganak2.pdf),
and our [older paper](https://www.ijcai.org/proceedings/2019/0163.pdf). You
can also check out our presentation
[PDF here](https://github.com/meelgroup/ganak-presentation/blob/main/ganak2-jul24-cav2025-zagreb.pdf).
Note that Ganak employs many ideas by many people. See AUTHORS for a list.

## Building
You can use Ganak online, straight [from your
browser](https://www.msoos.org/ganak/)

Use of the [release binaries](https://github.com/meelgroup/ganak/releases) is
_strongly_ encouraged, as Ganak requires a specific set of libraries to be
installed. The second best thing to use is Nix. Simply [install
nix](https://nixos.org/download/) and then:
```shell
nix shell github:meelgroup/ganak#ganak
```

Then you will have `ganak` binary available and ready to use.

If this is somehow not what you want, you can also build it. See the [GitHub
Action](https://github.com/meelgroup/ganak/actions/workflows/build.yml) for the
specific set of steps.

### Building statically

To build a static binary, you first need to build GMP with position-independent code
enabled. GMP's hand-optimised assembly is normally compiled without `-fPIC` (fine for
a native static binary), but `-fPIC` is required whenever the static `.a` is linked
into a shared object — for example, the Python extension (`.so`). Passing `--with-pic`
to GMP's `configure` makes both uses work:

```shell
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-static --enable-cxx --enable-shared --with-pic
make -j$(nproc)
sudo make install
cd ..
```

Then point CMake to the installed GMP static libraries (note: use `/usr/local/lib/`,
not a custom build directory, as those may be compiled for the wrong architecture):

```shell
mkdir build && cd build
cmake -DBUILD_SHARED_LIBS=OFF \
    -DGMPXX_LIBRARY=/usr/local/lib/libgmpxx.a \
    -DGMP_INCLUDE_DIR=/usr/local/include \
    ..
make -j$(nproc)
```

## Usage
Ganak takes a CNF in a special, DIMACS-like format as specified by the model
counting competition
[guidelines](https://mccompetition.org/assets/files/mccomp_format_24.pdf).
Basically, the format is as follows:
```plaintext
c t pwmc
p cnf 3 2
c p weight 1 0.3 0
c p weight -1 0.8 0
1 2 3 0
-1 2 0
c p show 1 2 0
```
The first line says it's a projected weighted model counting instance. The
second line says it has 3 variables and 2 clauses. The third and fourth lines
specify the weights of the variables 1. The weight of the literal 1 is 0.3 and
the weight of the literal -1 is 0.8. The weight of all unspecified variables is
1 for both positive and negative literals.

The last line says that the projection set of the counter is only variables 1
and 2. If no projection set is given, then the counter does an unprojected
count, i.e. all variables are assumed to be in the projection set.

Beware to ALWAYS give the weight of both the literal and its negation or
different counters may give different results.

### Accepted header directives

| Directive | Meaning |
|-----------|---------|
| `c t mc` | Unweighted, unprojected model counting (default) |
| `c t pmc` | Unweighted projected model counting |
| `c t wmc` | Weighted model counting |
| `c t wpmc` / `c t pwmc` | Weighted projected model counting (aliases) |
| `c p show v1 v2 ... 0` | Projection (sampling) set |
| `c p weight L W 0` | Weight `W` for literal `L` (positive or negative integer) |

`c t wmc` and `c t pwmc` require a weighted field generator (i.e. not
`--mode 0`); otherwise the parser will reject the file.

## Weights
It is _highly_ encouraged to give both the positive and the negative literal's
weight, e.g. `1` and `-1`:
```plain
c p weight 1 0.3 0
c p weight -1 0.8 0
```

This is important because different counters have different defaults,
so in order to have consistent results, you should must specify the weights of
both the positive and the negative literals.

There are many different ways you can specify the weights, the following are
all valid:
```plain
c p weight 1 1/4 0
c p weight -1 3e-5 0
c p weight 2 -9 0
c p weight -2 -3.4e5 0
```

Notice that in `--mode 1` Ganak will give exact rational counts, so you can
_always_ be sure that the count is exact. This is _very_ unlike all other
counters, which use floating point numbers to represent the weights, and
therefore you can only guess whether the value they show is correct or not.
While many claim to have "infinite precision", in reality, only their error
can be shown to be potentially infinite.

## Approximate Counting
Ganak integrate ApproxMC and can switch to approximate counting after a timeout.
This is done by using the `--appmct` flag. For example, if you want to run
Ganak for 1000 seconds and then switch to approximate counting, you can do:
```plaintext
./ganak --appmct 1000 --mode 0 <file>
```

Note that this can _only_ be used with `--mode 0`, i.e. in unweighted
counting.

## Probabilistic Counting
By default, Ganak uses a probabilistic caching of component counts, which means
that in extremely rare cases, often less than 1 case per billion (depending on the
problem), it can return incorrect count. The probability of the wrong count is displayed
at the end of solving with:
```plaintext
c s pac guarantees epsilon: 0 delta: 6.45757296e-10
```
which means that the probability of the wrong count is at most
`6.45757296e-10`, i.e. less than 1 in a billion.

If you must have a non-probabilistic count, you can use the `--prob 0` flag.

## Python Package (pyganak)

Ganak is available as a Python package on [PyPI](https://pypi.org/project/pyganak/):

```bash
pip install pyganak
```

Pre-built wheels are available for Linux (x86-64, ARM64) and macOS (Apple Silicon, Intel).

### Unweighted counting

```python
from pyganak import Counter

c = Counter()
c.add_clause([1, 2])      # x1 OR x2
c.add_clause([-1, 2])     # NOT x1 OR x2
print(c.count())          # → 2  (exact Python int, arbitrary precision)
```

### Weighted counting

```python
from pyganak import WeightedCounter

c = WeightedCounter()
c.add_clause([1, 2])           # x1 OR x2

# Set weights for both polarities of each variable.
c.set_lit_weight( 1, 0.3)     # weight of  x1 = 0.3
c.set_lit_weight(-1, 0.7)     # weight of ¬x1 = 0.7
c.set_lit_weight( 2, 0.4)     # weight of  x2 = 0.4
c.set_lit_weight(-2, 0.6)     # weight of ¬x2 = 0.6

# Models: (T,T)=0.12  (T,F)=0.18  (F,T)=0.28  → total=0.58
print(c.count())               # → 0.58  (Python float)
```

Weights are supplied as Python `float` (double) values.  Internally,
`WeightedCounter` uses [MPFR](https://www.mpfr.org/) floating-point
arithmetic at configurable precision (default 128 bits, roughly 38 significant
decimal digits).  However, because floating-point arithmetic is **not
associative**, the result is a high-precision approximation rather than an
exact value — the order in which terms are accumulated can affect the last few
bits of the result.  For exact rational weighted counting use `--mode 1` from
the command line (or the C++ library with `FGenMpq`).

The `prec` constructor argument controls the MPFR precision in bits:

```python
c = WeightedCounter(prec=256)   # 256-bit internal precision
```

## Using as a Library
Ganak can be used as a library. The file `src/example.cpp` gives an example of
how to use Ganak as a library. Let's go through it step by step:
```cpp
  std::unique_ptr<CMSat::FieldGen> fg = std::make_unique<ArjunNS::FGenMpq>();
  ArjunNS::SimplifiedCNF cnf(fg);
  cnf.new_vars(10);
  vector<CMSat::Lit> cl;
  cl.push_back(mklit(1));
  cl.push_back(mklit(2));
  cnf.add_clause(cl);
  cnf.set_weighted(true);
  cnf.set_lit_weight(mklit(4), ArjunNS::FMpq(10));
  cnf.set_lit_weight(mklit(-4), ArjunNS::FMpq(1));
  cnf.set_sampl_vars(vector<uint32_t>{mklit(1).var(), mklit(2).var(), mklit(4).var()});

  run_arjun(cnf);
  conf.verb = 0;
  Ganak counter(conf, fg);
  setup_ganak(cnf, counter);

  auto cnt = cnf.get_multiplier_weight()->dup();
  if (!cnf.get_multiplier_weight()->is_zero()) *cnt *= *counter.count();
  cout << "count is: " << std::fixed << *cnt << endl;
```

Here, we first create the mathematical field we'll be counting over. In this
case, we use the field of rationals. Then we create a CNF with 10 variables and
the clause `(x1 v x2)`, with a projection set of `(x1,x2,x4)`. This CNF has a
count of 33, because over (x1,x2) there are 3 solutions, and for each solution,
x4 can be either true or false, with a combined weight of 11 (=10+1). This
makes it a total count of 3 * 11 = 33. We then push this CNF through Arjun's
simplification, and then through Ganak's counting, finally, we combine these
two systems' counts to get the final count.

## Supported Weights

| Mode | Flag | Field | Weight format example | Notes |
|------|------|-------|-----------------------|-------|
| 0 | `--mode 0` | Integer | _unweighted_ | Default; weights not supported. Supports `--appmct` |
| 1 | `--mode 1` | Rational (exact) | `c p weight 1 1/4 0` | No floating-point error |
| 2 | `--mode 2` | Complex rational | `c p weight 9 1/2+4i 0` | Must give both parts: `a+bi` or `a-bi` |
| 3 | `--mode 3` | Polynomial over rationals | `c p weight 9 1/2*x0+4*x1+2 0` | Requires `--npolyvars N` |
| 4 | `--mode 4` | Parity (mod 2) | `c p weight 1 1 0` | |
| 5 | `--mode 5` | Integer mod prime | `c p weight 1 3 0` | Requires `--prime X` |
| 6 | `--mode 6` | Complex float (MPFR) | `c p weight 9 1/2+4i 0` | Must give both parts: `a+bi` or `a-bi`; see `--mpfrprec` |
| 7 | `--mode 7` | Real float (MPFR) | `c p weight 1 0.3 0` | See `--mpfrprec` |
| 13 | `--mode 13` | Laurent polynomial over rationals | `c p weight 9 1/2*z0^2-3*z1^-1 0` | Requires `--npolyvars N`; exponents may be negative |

### Laurent polynomial counting (`--mode 13`)

Mode 13 is like the polynomial mode (`--mode 3`) but the weights are
*multivariate Laurent polynomials* over the rationals, i.e. polynomials whose
variable exponents may be **negative** as well as non-negative. Pass the number
of polynomial variables with `--npolyvars N`, and refer to them as `z0`, `z1`,
..., `z(N-1)` in the weight lines. Terms are written `coeff*zI^e` (the
exponent `e` may be negative, and parentheses around it are optional, e.g.
`z0^-2` or `z0^(-2)`); separate terms with `+`/`-`. A single weight line such as
`c p weight 9 1/2*z0^2 - 3*z1^-1 + z0^(-2) 0` is a valid Laurent polynomial.

#### Worked example

The following CNF has one clause `z0 ∨ z1` and assigns a Laurent-polynomial
weight to each literal:

```
p cnf 2 1
c t wmc
1 2 0
c p weight 1 z0+1 0
c p weight -1 z0^-1 0
c p weight 2 z1+1 0
c p weight -2 z1^-1 0
```

Running it:

```
./ganak --mode 13 --npolyvars 2 example.cnf
```

produces (summing the weight products over the three satisfying assignments):

```
s SATISFIABLE
c s exact laurent z0*z1 + z0 + z1 + z0*z1^-1 + 1 + z0^-1*z1 + z1^-1 + z0^-1
```

The count is printed on the `c s exact laurent` line as the resulting Laurent
polynomial.

You can also write your own field by implementing the `Field` and `FieldGen`
interfaces. Absolutely _any_ field will work, and it's as easy as implementing
`+,-,*` and `/` operators, and the `0` and `1` constants. It's a fun
exercise to do.

## Compiling to a d-DNNF circuit

Ganak can compile its search trace into a
[d-DNNF](https://en.wikipedia.org/wiki/Decision-DNNF) circuit in the d4 `.nnf`
format, for reuse in repeated queries, model enumeration, or functional
synthesis.

**1. Dump** (still prints the count; streamed to disk, so it is not held in
memory). The output is a faithful d-DNNF:

```shell
./ganak --compile out.nnf in.cnf
```

For a projected CNF (inputs `< indep_support_end`, outputs `>= indep_support_end`),
the SAT oracle records one example assignment of the output vars at each SAT
leaf, so the same faithful circuit is also a sound Boolean functional-synthesis
compiler — read a Skolem witness off the circuit with
`tests/ddnnf_verify.synthesize()`.

**2. Fix up.** The raw file is correct from the root but keeps Ganak's internal
ids and may carry unreachable "dead" nodes. `ddnnf-cleanup` drops them and
renumbers `root = 1`, contiguous, yielding a strict d4 file:

```shell
./ddnnf-cleanup out.nnf clean.nnf      # omit 2nd arg (or "-") to write to stdout
```

**3. Verify.** `tests/ddnnf_verify.py` prints the circuit's structural count,
which must match Ganak's. It accepts both the raw and cleaned files:

```shell
python3 tests/ddnnf_verify.py clean.nnf
```

## Fuzzing
We use the [SharpVelvet](https://github.com/meelgroup/SharpVelvet) model counter
fuzzer, developed by [Anna Latour](https://scholar.google.com/citations?user=nf5lfegAAAAJ&hl=nl)
to check correctness of Ganak.

## How to Cite
```
@inproceedings{SM2025,
title={Engineering an Efficient Probabilistic Exact Model Counter},
author={Soos, Mate and Meel, Kuldeep S.},
booktitle={Proceedings of the International Conference on Computer Aided Verification (CAV)}
year={2025}
}

@inproceedings{SRSM19,
title={GANAK: A Scalable Probabilistic Exact Model Counter},
author={Sharma, Shubham and  Roy, Subhajit and  Soos, Mate and  Meel, Kuldeep S.},
booktitle={Proceedings of International Joint Conference on Artificial Intelligence (IJCAI)},
year={2019}
}
```
