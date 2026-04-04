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
SAT_DIR=$(cd ../.. && pwd)
cmake -DBUILD_SHARED_LIBS=OFF \
    -DGMPXX_LIBRARY=/usr/local/lib/libgmpxx.a \
    -Dcryptominisat5_DIR="${SAT_DIR}/cryptominisat/build" \
    -Dsbva_DIR="${SAT_DIR}/sbva/build" \
    -DEvalMaxSAT_DIR="${SAT_DIR}/EvalMaxSAT/build" \
    -Dtreedecomp_DIR="${SAT_DIR}/treedecomp/build" \
    -Darjun_DIR="${SAT_DIR}/arjun/build" \
    -Dapproxmc_DIR="${SAT_DIR}/approxmc/build" \
    -Dcadical_DIR="${SAT_DIR}/cadical/build" \
    -Dcadiback_DIR="${SAT_DIR}/cadiback/build" \
    ..
make -j$(nproc)
strip ganak
