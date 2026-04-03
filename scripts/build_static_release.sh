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
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
    -DGMP_LIBRARY=/usr/local/lib/libgmp.a \
    -DGMPXX_LIBRARY=/usr/local/lib/libgmpxx.a \
    -Dcadical_DIR="${SAT_DIR}/cadical/build" \
    -Dcryptominisat5_DIR="${SAT_DIR}/cryptominisat/build" \
    -Dsbva_DIR="${SAT_DIR}/sbva/build" \
    -DEvalMaxSAT_DIR="${SAT_DIR}/EvalMaxSAT/build" \
    -Dtreedecomp_DIR="${SAT_DIR}/treedecomp/build" \
    -Darjun_DIR="${SAT_DIR}/arjun/build" \
    -Dapproxmc_DIR="${SAT_DIR}/approxmc/build" \
    ..
make -j$(nproc)
strip ganak

# Get SHA1 hashes and build versioned filename
HASHES=$(echo "" | ./ganak 2>&1)
GANAK_SHA=$(echo "$HASHES" | grep "Ganak SHA1" | grep -oP '[0-9a-f]{40}' | cut -c1-8)
ARJUN_SHA=$(echo "$HASHES" | grep "Arjun SHA1" | grep -oP '[0-9a-f]{40}' | cut -c1-8)
CMS_SHA=$(echo "$HASHES" | grep "CMS SHA1" | grep -oP '[0-9a-f]{40}' | cut -c1-8)
APPROXMC_SHA=$(echo "$HASHES" | grep "ApproxMC SHA1" | grep -oP '[0-9a-f]{40}' | cut -c1-8)

DEST="ganak_${GANAK_SHA}_${ARJUN_SHA}_${APPROXMC_SHA}_${CMS_SHA}"
cp ganak "$DEST"
echo "Copied to $DEST"
