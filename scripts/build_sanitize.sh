#!/bin/bash

rm -rf CMake* src cmake* ganak* sharp* Make*
CC=clang CXX=clang++ cmake -DSANITIZE=ON ..
make -j6


echo "Now run with:"
echo "UBSAN_OPTIONS=print_stacktrace=1 ./ganak"
