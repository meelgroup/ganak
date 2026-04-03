#!/bin/bash

rm -rf make* CMake* src cmake* ganak* sharp* Make*
cmake -DCMAKE_BUILD_TYPE=Debug \
    -Dcadical_DIR=../cadical \
    -Dcadiback_DIR=../cadiback \
    -Dcryptominisat5_DIR=../cryptominisat \
    -Dsbva_DIR=../sbva \
    -Dtreedecomp_DIR=../treedecomp \
    -DEvalMaxSAT_DIR=../EvalMaxSAT \
    -Darjun_DIR=../arjun \
    -Dapproxmc_DIR=../approxmc \
    ..
make -j$(nproc) VERBOSE=1
