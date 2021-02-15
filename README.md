# SYMGANAK- Symmetric Component Caching for Model Counting on Combinatorial Instances
SYMGANAK  takes in a CNF formula `F` and a confidence `delta` as input and returns `count` such that `count` is the number of solutions of `F` with confidence at least `1 - delta`. 

To read more about technical algorithms in Symganak, please refer to [our paper](https://www.comp.nus.edu.sg/~meel/Papers/aaai21-bdsrm.pdf) 

## Installation

### Compiling in Linux

To build symganak, issue:

```
bash
sudo apt-get install libgmp-dev
sudo apt-get install libmpfr-dev
sudo apt-get install libmpc-dev
mkdir build && cd build
cmake ..
make
cp symganak ../bin
```


## Usage

You can use the run_symganak.sh script in scripts directory to run symganak. A simple invocation looks as follows:

```bash
./run_symganak.sh <cnffile>
```
To use different settings of parameters modify the helper run_symganak.sh script. The usage instructions and default values to parameters can be found by running:

```bash
../bin/symganak
```

## Benchmarks
Few toy benchmarks are given in benchmarks directory. The benchmarks used in the paper are present [here](https://github.com/VincentDerk/Paper-SymGANAK-benchmark).


## Issues, questions, bugs, etc.
Please click on "issues" at the top and [create a new issue](https://github.com/meelgroup/ganak/issues). All issues are responded to promptly.

## How to Cite
```
@inproceedings{BDSRM21,
author={van Bremen, Timothy and  Derkinderen, Vincent and  Sharma, Shubham and  Roy, Subhajit and  Meel, Kuldeep S},
title={Symmetric Component Caching for Model Counting on Combinatorial Instances},
booktitle={Proceedings of AAAI Conference on Artificial Intelligence (AAAI)},
month={2},
year={2021},
}
```

## Contributors
  * Vincent Derkinderen (vincent.derkinderen@cs.kuleuven.be )
  * Timothy Van Bremen (timothy.vanbremen@cs.kuleuven.be)
  * Shubham Sharma (sharma.shubham736@gmail.com)
  * Subhajit Roy (subhajit@iitk.ac.in)
  * Kuldeep Meel (meel@comp.nus.edu.sg)
