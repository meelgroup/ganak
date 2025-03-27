# Ganak2 -- A Probabilistic Exact Model Counter
Ganak takes in a CNF formula and returns `count` such that `count` is the
number of solutions of `F` with confidence at least `1 - delta`. Delta is fixed to
approx 2^-40.

To read more about Ganak-specific ideas, please refer to [our
paper](https://www.comp.nus.edu.sg/~meel/Papers/ijcai19srsm.pdf). Note that
Ganak employs many ideas by many people. See AUTHORS for a list.

## Building

Use of the [release binaries](https://github.com/meelgroup/ganak/releases) is
strongly encouraged, as Ganak requires a specific set of libraries to be
installed. The second best thing to use is Nix. Simply [install
nix](https://nixos.org/download/) and then:

```plaintext
git clone https://github.com/meelgroup/ganak
nix-shell
```
If this is somehow not what you want, you can also build it. See the GitHub
Action for the specific set of steps. At the end of this README you can
find a detailed set of instructions.

## Usage
Ganak takes a CNF in a special, DIMACS-like format as specified by the
model counting competition [guidelines](https://mccompetition.org/assets/files/mccomp_format_24.pdf).
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

## Benchmarks
Model Counting Competition benchmarks are available on the [Model Counting
Competition website](http://https://mccompetition.org/).

Our benchmarks are available
[here](https://drive.google.com/file/d/15dUJI55drFH_0-4-qWjoF_YR0amb3xnK/view?usp=sharing)

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

## Building from source

```bash
sudo apt-get install -yq libgmp-dev libmpfr-dev

git clone https://github.com/meelgroup/cadical
cd cadical
git checkout add_dynamic_lib
./configure
make
cd ..

git clone https://github.com/meelgroup/cadiback
cd cadiback
git checkout synthesis
./configure
make
cd ..

git clone https://github.com/msoos/cryptominisat
cd cryptominisat
git checkout synthesis
mkdir build && cd build
cmake ..
make
cd ../..

git clone https://github.com/meelgroup/sbva
cd sbva
mkdir build && cd build
cmake ..
make
cd ../..

git clone https://github.com/meelgroup/breakid
cd breakid
mkdir build && cd build
cmake ..
make
cd ../..

git clone https://github.com/meelgroup/approxmc
cd approxmc
mkdir build && cd build
git checkout synthesis
cmake ..
make
cd ../..

git clone https://github.com/meelgroup/ganak
cd ganak
mkdir build && cd build
cmake ..
make
```

if you also need to have the system run on older or different machines, you
need to have it compile for Sandybridge architecture. To do this, you need to
 recompile GMP and MPFR with the right flags. Here is how to do it:

```bash
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar xvf gmp-6.3.0.tar.xz
cd gmp-6.3.0
CFLAGS="-march=sandybridge" ./configure --enable-cxx
make -j12
sudo make install
cd ..
wget https://ftp.gnu.org/gnu/mpfr/mpfr-4.1.0.tar.xz
tar xvf mpfr-4.1.0.tar.xz
cd mpfr-4.1.0
CFLAGS="-march=sandybridge" ./configure
make -j12
sudo make install
cd ..
```
