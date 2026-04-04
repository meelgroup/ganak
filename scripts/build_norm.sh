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
rm -rf deps
rm -rf _deps
SAT_DIR=$(cd ../.. && pwd)
cmake -DENABLE_TESTING=ON \
    -Dcryptominisat5_DIR="${SAT_DIR}/cryptominisat/build" \
    -Dsbva_DIR="${SAT_DIR}/sbva/build" \
    -Dtreedecomp_DIR="${SAT_DIR}/treedecomp/build" \
    -DEvalMaxSAT_DIR="${SAT_DIR}/EvalMaxSAT/build" \
    -Darjun_DIR="${SAT_DIR}/arjun/build" \
    -Dapproxmc_DIR="${SAT_DIR}/approxmc/build" \
    -Dcadical_DIR="${SAT_DIR}/cadical/build" \
    -Dcadiback_DIR="${SAT_DIR}/cadiback/build" \
    ..
make -j$(nproc)
make test
