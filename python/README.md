# pyganak — Python bindings for Ganak

`pyganak` provides Python bindings for [Ganak](https://github.com/meelgroup/ganak),
a high-performance **exact model counter** for CNF formulas.  It supports
both plain model counting (how many satisfying assignments does the formula
have?) and *projected* model counting (count only over a given set of
variables, existentially quantifying the rest).

Arjun preprocessing is applied automatically before each count, just as the
command-line tool does.

## Installation

```bash
pip install pyganak
```

Pre-built wheels are available for Linux (x86-64, ARM64) and macOS (Apple
Silicon, Intel).

## Quick start

```python
from pyganak import Counter

# Create a counter
c = Counter()

# Add clauses as lists of nonzero signed integers (1-indexed variables).
# Positive literal = positive variable, negative literal = negated variable.
c.add_clause([1, 2])      # x1 OR x2
c.add_clause([-1, 2])     # NOT x1 OR x2

print(c.count())          # → 2
```

## API

### `Counter(verbose=0, seed=0)`

Create a new counter instance.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `verbose` | `int` | `0` | Verbosity level (0 = silent) |
| `seed`    | `int` | `0` | Random seed for the solver |

---

### `add_clause(clause)`

Add a single clause.

- `clause` — iterable of nonzero signed integers.  Positive values are
  positive literals; negative values are negated literals.  Variables are
  **1-indexed**.

```python
c.add_clause([1, -2, 3])   # x1 OR NOT x2 OR x3
```

---

### `add_clauses(clauses)`

Add multiple clauses at once.

```python
c.add_clauses([
    [1, 2],
    [-1, 3],
])
```

---

### `set_sampling_set(vars)`

Set the *sampling set* (independent support / projection set).

- `vars` — iterable of positive integers (1-indexed variable numbers).

When this is called, only the given variables are counted; all other
variables are existentially quantified out.  If you never call this method,
all variables are included (plain model counting).

```python
c.add_clause([1, 2, 3])
c.set_sampling_set([1, 2])   # count projected onto x1, x2
print(c.count())             # → 4  (all four assignments of x1,x2 work)
```

---

### `count() → int`

Run Arjun preprocessing and Ganak model counting.  Returns the exact model
count as a Python `int` (arbitrary precision).

`count()` may only be called **once** per `Counter` instance.  Calling it
a second time raises `RuntimeError`.  To count a modified formula, create a
new `Counter` and add the desired clauses to it.

---

### `new_vars(n)`

Declare `n` new variables, extending the variable universe by `n`.
Useful when you want to count over a known number of variables without
adding any clauses for them.

```python
c = Counter()
c.new_vars(10)
print(c.count())   # 2^10 = 1024  (no constraints → all assignments valid)
```

### `nof_vars() → int`

Return the number of variables seen so far (highest variable index).

### `nof_clauses() → int`

Return the number of clauses added so far.

## Examples

### Plain model counting

```python
from pyganak import Counter

# (x1 XOR x2): 2 models
c = Counter()
c.add_clause([1, 2])
c.add_clause([-1, -2])
print(c.count())   # 2
```

### Projected model counting

```python
from pyganak import Counter

# (x1 OR x2 OR x3) projected onto {x1, x2}: 4 models
c = Counter()
c.add_clause([1, 2, 3])
c.set_sampling_set([1, 2])
print(c.count())   # 4
```

### Large count (arbitrary precision)

```python
from pyganak import Counter

# 30 unconstrained variables → 2^30 models
c = Counter()
for i in range(1, 31):
    c.add_clause([i, -i])   # tautological clause — declares each variable
print(c.count())            # 1073741824
```

## Building from source

```bash
# Standard pip install (scikit-build-core handles everything)
pip install .
```

Requires: GMP ≥ 5, MPFR ≥ 3, FLINT ≥ 2.  On Ubuntu/Debian:

```bash
sudo apt install libgmp-dev libmpfr-dev libflint-dev
```

On macOS:

```bash
brew install gmp mpfr flint
```

## Development testing (after building with CMake)

If you have already built Ganak with CMake (see the top-level `README.md`),
the extension is available as `build/lib/pyganak*.so`.  You can test it
without a full `pip install` by creating a venv and pointing `PYTHONPATH` at
that directory:

```bash
python3 -m venv venv
venv/bin/pip install pytest
PYTHONPATH=build/lib venv/bin/pytest python/tests/ -v
```

## License

MIT — see `LICENSE.txt` in the repository root.
