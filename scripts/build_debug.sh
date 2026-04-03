#!/bin/bash

rm -rf make*
rm -rf CMake*
rm -rf src
rm -rf cmake*
rm -rf ganak*
rm -rf sharp*
rm -rf Make*
rm -rf deps
rm -rf _deps
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
