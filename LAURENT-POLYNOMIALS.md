# Laurent polynomial counting (`--mode 13`)

Mode 13 is like the polynomial mode (`--mode 3`) but the weights are
*multivariate Laurent polynomials* over the rationals, i.e. polynomials whose
variable exponents may be **negative** as well as non-negative. Pass the number
of polynomial variables with `--npolyvars N`, and refer to them as `z0`, `z1`,
..., `z(N-1)` in the weight lines. Terms are written `coeff*zI^e` (the
exponent `e` may be negative, and parentheses around it are optional, e.g.
`z0^-2` or `z0^(-2)`); separate terms with `+`/`-`. A single weight line such as
`c p weight 9 1/2*z0^2 - 3*z1^-1 + z0^(-2) 0` is a valid Laurent polynomial.

## Worked example

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
