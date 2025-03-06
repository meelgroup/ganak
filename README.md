# GANAK- A Probabilistic Exact Model Counter
GANAK  takes in a CNF formula `F` and a confidence `delta` as input and returns `count` such that `count` is the number of solutions of `F` with confidence at least `1 - delta`. GANAK supports projected model counting (see below). 

To read more about technical algorithms in Ganak, please refer to [our paper](https://www.cs.toronto.edu/~meel/Papers/ijcai19srsm.pdf) 

## Installation

### Compiling in Linux

To build ganak, issue:

```bash
sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev
git clone https://github.com/meelgroup/ganak
mkdir build && cd build
cmake ..
make
```

For Mac, pass option `cmake -DDOPCC=OFF ...` to cmake

### Model Counting

To count, run:

```bash
cd build
./ganak myfile.cnf
```

### Projected Model Counting
For some applications, one is not interested in solutions over all the variables and instead interested in counting the number of unique solutions to a subset of variables, called sampling set. GANAK allows you to specify the sampling set using the following modified version of DIMACS format:

```
$ cat myfile.cnf
p cnf 500 1
c ind 1 3 4 6 7 8 10 0
3 4 0
```
Above, using the `c ind` line, we declare that only variables 1, 3, 4, 6, 7, 8 and 10 form part of the sampling set out of the CNF's 500 variables `1,2...500`. This line must end with a 0. The solution that GANAK will be giving is essentially answering the question: how many different combination of settings to this variables are there that satisfy this problem? Naturally, if your sampling set only contains e.g. 7 variables, then the maximum number of solutions can only be at most `2**7 = 128`. This is true even if your CNF has thousands of variables. Here, we only have `2**5*3` solutions however, since we ban `-3,-4` solution from the tuple `(3,4)`.

Note: By default if sampling set is present ganak will do projected model counting, to turn off projected model counting use the -noPMC flag.

## Benchmarks
Few toy benchmarks are given in benchmarks directory. Full list of benchmarks used for our experiments is available [here](https://drive.google.com/file/d/15dUJI55drFH_0-4-qWjoF_YR0amb3xnK/view?usp=sharing)


## Issues, questions, bugs, etc.
Please click on "issues" at the top and [create a new issue](https://github.com/meelgroup/ganak/issues). All issues are responded to promptly.

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

## Contributors
  * Shubham Sharma (sharma.shubham736@gmail.com)
  * Mate Soos (soos.mate@gmail.com)
  * Subhajit Roy (subhajit@iitk.ac.in)
  * Kuldeep Meel (meel@comp.nus.edu.sg)
