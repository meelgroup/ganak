#!/bin/bash

rm -rf CMake* src cmake* ganak* sharp* Make*
CC=clang CXX=clang++ cmake -DSANITIZE=ON \
    -Dcadical_DIR=../cadical \
    -Dcadiback_DIR=../cadiback \
    -Dcryptominisat5_DIR=../cryptominisat \
    -Dsbva_DIR=../sbva \
    -Dtreedecomp_DIR=../treedecomp \
    -DEvalMaxSAT_DIR=../EvalMaxSAT \
    -Darjun_DIR=../arjun \
    -Dapproxmc_DIR=../approxmc \
    ..
make -j$(nproc)


echo "Now run with:"
echo "UBSAN_OPTIONS=print_stacktrace=1 ./ganak"
