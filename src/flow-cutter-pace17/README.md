# FlowCutter PACE 2017 Submission

This repository contains the FlowCutter code submitted to the [PACE 2017](https://pacechallenge.wordpress.com/2016/12/01/announcing-pace-2017/) tree decomposition challenge. 
FlowCutter was developed at [KIT](https://www.kit.edu) in the [group of Prof. Dorothea Wagner](https://i11www.iti.kit.edu/).

If you are running a Unix-like system, then getting started is very simple. Just clone the repository and build the programs, as follows:

```bash
git clone https://github.com/kit-algo/flow-cutter-pace17.git
cd flow-cutter-pace17
./build.sh
```

There are no dependencies beyond a GCC with version 4.8 or newer. Clang should also work but has not been tested by us. Building the code under Windows probably requires a few code modifications in `pace.cpp`.

After executing the build script, the root directory of the repository should contain the two binary files `flow_cutter_pace17` and `flow_cutter_parallel_pace17`. These are the programs entered into the heuristic, sequential and heuristic, parallel tracks of the competition. The outputted decompositions are guaranteed to be valid but do not necessarily have a minimum width. Both executable have the same interface. 

There are three ways to correctly invoke the program:

```bash
./flow_cutter_pace17 < my_graph.gr 
./flow_cutter_pace17 my_graph.gr
./flow_cutter_pace17 -s 42 < my_graph.gr
```

The first and the last commands read the input graph from the standard input. The second command reads it from a file whose name is given as parameter. The `-s` parameter sets the random seed. By default a seed of 0 is assumed. We tried to make sure that given the same seed, the behaviour of the sequential binary should be the identical even accross compilers.

The executables run until either a SIGINT or SIGTERM signal is sent. Once this signal is encountered the programm prints a tree decomposition to the standard output with the smallest width that it could found and terminates. Note that no decomposition is outputted if you send the signal before any decomposition is found.

The format specification of the input graph and output decompositions follow those of the [PACE 2017](https://pacechallenge.wordpress.com/2016/12/01/announcing-pace-2017/) challenge. 

**Warning:** The FlowCutter PACE 2017 code only optimizes the tree width. If you want to optimize different criteria such as the fill-in, the [FlowCutter PACE 2016 code](https://github.com/ben-strasser/flow-cutter-pace16) will probably give better results.

## Publications

Please cite the following article if you use our code in a publication:

* Graph Bisection with Pareto-Optimization.
  Michael Hamann and Ben Strasser.
  Proceedings of the 18th Meeting on Algorithm Engineering and Experiments (ALENEX'16).

## Contact

Please send an email to Ben Strasser (strasser (at) kit (dot) edu).

