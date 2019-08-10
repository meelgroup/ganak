# GANAK- model counter
GANAK is a probabilistic exact model counter. It takes in a CNF formula `F` and a confidence `delta` as input and returns `count` such that `count` is the number of solutions of `F` with confidence at least `1 - delta`.

## Installation
The static binary of ProbsharpSAT and MIS is given in bin directory (tested on ubuntu 18.04, Fedora 21 and CentOS 6.9).

### Compiling in Linux

build MIS from [here](https://github.com/meelgroup/mis) and add `mis.py`, `muser2` and `togmus` to `build` directory.

To build ganak, issue:

```
bash
sudo apt-get install libgmp-dev
sudo apt-get install libmpfr-dev
sudo apt-get install libmpc-dev
mkdir build && cd build
cmake ..
make
cp ../bin/ganak.py .
```


## Usage
You can run GANAK by using 'ganak.py' Python script. A simple invocation looks as follows:
```bash
python3 ganak.py <cnffile>
```

The usage instructions and default values to arguments can be found by running:
```bash
python3 ganak.py -h
```

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
