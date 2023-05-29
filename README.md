# GANAK- A Probabilistic Exact Model Counter
GANAK  takes in a CNF formula `F` and a confidence `delta` as input and returns `count` such that `count` is the number of solutions of `F` with confidence at least `1 - delta`. GANAK supports projected model counting (see below). 

To read more about Ganak-specific ideas, please refer to [our paper](https://www.comp.nus.edu.sg/~meel/Papers/ijcai19srsm.pdf). Note that Ganak employs many ideas by many people. See AUTHORS for a list.

## Compiling in Linux

To build ganak, you first need to build [CryptoMiniSat](github.com/msoos/cryptominisat). Then issue:

```bash
sudo apt-get install libgmp-dev
sudo apt-get install libmpfr-dev
sudo apt-get install libmpc-dev
mkdir build && cd build
cmake ..
make
```

## Usage
Ganak takes a CNF in DIMACS form as an input:

```bash
./ganak myfile.cnf
```

## Benchmarks
Few toy benchmarks are given in benchmarks directory. Full list of benchmarks used for our experiments is available [here](https://drive.google.com/file/d/15dUJI55drFH_0-4-qWjoF_YR0amb3xnK/view?usp=sharing)


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

## Copyright, contributors
GANAK licensed under MIT. Flow Cutter licensed under BSD2-style license. Binary built contains both (they are compatible with each other).

See AUTHORS for full list. Ganak-specific contributions by, in alphabetical order: Kuldeep S. Meel, Subhajit Roy, Shubham Sharma, Mate Soos. However, MANY people's ideas contributed to Ganak. It's easy to list everyone from Stockmayer (e.g. probing) to Donald Knuth (e.g. memops).
