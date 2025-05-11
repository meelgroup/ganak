# Ganak2, A Probabilistic Exact Model Counter
Ganak is a state-of-the art probabilistic exact model counter that won _all_
available awards at the 2024 Model Counting Competition. It can count over ANY
field, including but not limited to integers, rationals, complex numbers,
integers modulo prime, and polynomials over a finite field.

To read more about Ganak-specific ideas, please refer to [our
paper](https://www.comp.nus.edu.sg/~meel/Papers/ijcai19srsm.pdf). Note that
-Ganak employs many ideas by many people. See AUTHORS for a list.

## Building
Use of the [release binaries](https://github.com/meelgroup/ganak/releases) is
_strongly_ encouraged, as Ganak requires a specific set of libraries to be
installed. The second best thing to use is Nix. Simply [install
nix](https://nixos.org/download/) and then:
```plaintext
git clone https://github.com/meelgroup/ganak
cd ganak
nix-shell
```

Then you will have `ganak` binary available and ready to use.

If this is somehow not what you want, you can also build it. See the [GitHub
Action](https://github.com/meelgroup/ganak/actions/workflows/build.yml) for the
specific set of steps.

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

## Different Modes of Counting
Ganak supports many different ways of counting:
- For counting over integers, the default mode `--mode 0` works. Here, you can
  even run approximate counting after some timeout, with e.g `--appmct 1000`
  which will automatically switch over to approximate counting after 1000s of
  preprocessing and Ganak.
- For counting over rationals, you can use `--mode 1` which will give you a
  precise count, without _any_ floating point errors.
- For counting over the complex field, you need to specify the weight of the
  literal as a complex number and pass the `--mode 2` flag. For example, if you
  want to give the weight 1/2 + 4i for literal 9, you can specify the weight as
  `c p weight 9 1/2 4 0`.
- For counting over polynomials over a finite field, you can use `--mode 3`
  which will give you a polynomial count. The weight of the literal is given as
  `c p weight 9 1/2*x0+4x1+2 0` which means the weight of literal 9 is
  `1/2*x_0 + 4*x_1+2`. In this mode you MUST specify the number of polynomial
  variables via `--npolyvars N`
- For parity counting, you can use the `--mode 4` flag. This will
  count the number of solutions modulo 2
- For counting over integers modulo a prime, you can use `--mode 5 --prime X`,
  where X is the prime.

You can also write your own field by implementing the `Field` and `FieldGen`
interfaces. Absolutely _any_ field will work, and it's as easy as implementing
+/-/* and / operators, and the 0 and 1 constants.

## Complex Numbers
We can now count the number of solutions of the above formula using Ganak:
```shell
$ ganak --verb 0 --mode 1 a.cnf
c o Total time [Arjun+GANAK]: 0.00
s SATISFIABLE
c s exact arb float 1.8999e+00
c o exact arb 19/10
```

We need to pass `--mode 1` because it's a weighted model counting instance. The count
is presented both in a floating point format and as a fraction. The fraction is
always precise.

## How to Cite
```
@inproceedings{SRSM19,
title={GANAK: A Scalable Probabilistic Exact Model Counter},
author={Sharma, Shubham and  Roy, Subhajit and  Soos, Mate and  Meel, Kuldeep S.},
booktitle={Proceedings of International Joint Conference on Artificial Intelligence (IJCAI)},
month={8},
year={2019}
}
```
