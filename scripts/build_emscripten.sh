#!/bin/bash
set -e

rm -rf lib* Test* tests* include tests CM* cmake* arjun Makefile rjun-src
emcmake cmake -DCMAKE_INSTALL_PREFIX=$EMINSTALL -DENABLE_TESTING=OFF \
    -Dcadical_DIR=../cadical \
    -Dcadiback_DIR=../cadiback \
    -Dcryptominisat5_DIR=../cryptominisat \
    -Dsbva_DIR=../sbva \
    -Dtreedecomp_DIR=../treedecomp \
    -DEvalMaxSAT_DIR=../EvalMaxSAT \
    -Darjun_DIR=../arjun \
    -Dapproxmc_DIR=../approxmc \
    ..
emmake make -j$(nproc)
emmake make install
cp ganak.wasm ../html
cp $EMINSTALL/bin/ganak.js ../html
