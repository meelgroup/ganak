# GANAK- model counter
GANAK is a probabilistic exact model counter. It takes in a CNF formula `F` and a confidence parameter `delta` as input and returns `count` such that `count` is the number of solutions of `F` with confidence at least `1 - delta`.

## Static Binary
The static binary of GANAK and MIS is given in static_bin directory (tested on ubuntu 18.04, Fedora 21 and CentOS 6.9).

## Running GANAK
You can run GANAK by using 'ganak.py' Python script present in static_bin directory. A simple invocation looks as follows:
```bash
python3 ganak.py <cnffile>
```

The usage instructions and default values to arguments can be found by running:
```bash
python3 ganak.py -h
```

## Benchmarks
Few toy benchmarks are given in benchmarks directory.

## Issues, questions, bugs, etc.
Please click on "issues" at the top and [create a new issue](https://github.com/meelgroup/ganak/issues). All issues are responded to promptly.
