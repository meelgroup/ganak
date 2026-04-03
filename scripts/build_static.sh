#!/usr/bin/env bash

set -euo pipefail

rm -rf CMake*
rm -rf src
rm -rf cmake*
rm -rf ganak*
rm -rf sharp*
rm -rf Make*
rm -rf Testing*
rm -rf tests*
rm -rf include
rm -rf .cmake
rm -rf lib*
rm -rf deps
rm -rf _deps
cmake -DSTATICCOMPILE=ON \
    -DGMP_LIBRARY=/usr/local/lib/libgmp.a \
    -DGMPXX_LIBRARY=/usr/local/lib/libgmpxx.a \
    -Dcadical_DIR=../../cadical/build \
    -Dcryptominisat5_DIR=../../cryptominisat/build \
    -Dsbva_DIR=../../sbva/build \
    -DEvalMaxSAT_DIR=../../EvalMaxSAT/build \
    -Dtreedecomp_DIR=../../treedecomp/build \
    -Darjun_DIR=../../arjun/build \
    -Dapproxmc_DIR=../../approxmc/build \
    ..
make -j$(nproc)
strip ganak
