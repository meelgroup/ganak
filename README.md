# Ganak2, A Probabilistic Exact Model Counter
Ganak takes in a CNF formula and returns `count` such that `count` is the
number of solutions of `F` with confidence at least `1 - delta`. Delta is fixed
to approx 2^-40.

To read more about Ganak-specific ideas, please refer to [our
paper](https://www.comp.nus.edu.sg/~meel/Papers/ijcai19srsm.pdf). Note that
-Ganak employs many ideas by many people. See AUTHORS for a list.

## Building
Use of the [release binaries](https://github.com/meelgroup/ganak/releases) is
strongly encouraged, as Ganak requires a specific set of libraries to be
installed. The second best thing to use is Nix. Simply [install
nix](https://nixos.org/download/) and then:
```plaintext
git clone https://github.com/meelgroup/ganak
nix-shell
```

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
