#!/usr/bin/env bash
# Builds ganak.js and bdnnf.js for the web page (html/) in build-wasm/ (the
# native build/ stays), against the dependencies' build-wasm/. Run from anywhere.
set -euo pipefail
source ~/emscripten.sh
REPO="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
SAT_DIR="$(cd "$REPO/.." && pwd)"
mkdir -p "$REPO/build-wasm"
cd "$REPO/build-wasm"
emcmake cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=OFF \
    -DCMAKE_INSTALL_PREFIX="${EMINSTALL}" \
    -Dcadical_DIR="${SAT_DIR}/cadical/build-wasm" \
    -Dcadiback_DIR="${SAT_DIR}/cadiback/build-wasm" \
    -Dcryptominisat5_DIR="${SAT_DIR}/cryptominisat/build-wasm" \
    -Dsbva_DIR="${SAT_DIR}/sbva/build-wasm" \
    -Dtreedecomp_DIR="${SAT_DIR}/treedecomp/build-wasm" \
    -Darjun_DIR="${SAT_DIR}/arjun/build-wasm" \
    -Dapproxmc_DIR="${SAT_DIR}/approxmc/build-wasm" \
    "$REPO"
emmake make -j"$(nproc)" ganak-bin bdnnf-bin
cp ganak.js ganak.wasm bdnnf.js bdnnf.wasm "$REPO/html/"
echo "$(git -C "$REPO" rev-parse --short HEAD)$(git -C "$REPO" diff --quiet HEAD -- || echo -dirty)-$(date +%Y%m%d%H%M%S)" > "$REPO/html/version.txt"
cd "$REPO/html" && tsc -p tsconfig.json && tsc -p tsconfig.worker.json
