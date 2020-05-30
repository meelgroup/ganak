# GANAK- A Probabilistic Exact Model Counter
GANAK  takes in a CNF formula `F` and a confidence `delta` as input and returns `count` such that `count` is the number of solutions of `F` with confidence at least `1 - delta`. GANAK supports projected model counting (see below). 

To read more about technical algorithms in Ganak, please refer to [our paper](https://www.comp.nus.edu.sg/~meel/Papers/ijcai19srsm.pdf) 

#### New: Support Added for Weighted Model Counting: Please checkout to wmc branch.

## Installation

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
```


## Usage
You can run GANAK by using 'ganak.py' Python script. A simple invocation looks as follows:
```bash
python3 ../bin/ganak.py <cnffile>
```

The usage instructions and default values to arguments can be found by running:
```bash
python3 ../bin/ganak.py -h
```
## Weight Format
GANAK supports providing weights in CNF itself. Weight of a literal is in [0,1], specified by line starting with 'w' <space> literal <space> weight. For example to give weight 0.5 to literal l use the following format:
`w l 0.5`

While weights for both positive and negative literals should be specified, if weight of only positive literal is specified, GANAK assumes it to be normalized and assigns weight of negative literal as 1 - weight(l). By default, weight of every literal l is set to 1 if the weight of l or -l is not given in CNF. Some examples are available in benchmarks directory for reference.


## Benchmarks
Few toy benchmarks are given in benchmarks directory.


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
