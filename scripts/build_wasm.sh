#!/usr/bin/env bash
# Builds ganak.js and bdnnf.js into ../html from a separate build-wasm/, against
# the wasm dependencies rebuild_emcripten_all.sh installed into $EMINSTALL; the
# native build/ is left alone. Usage: scripts/build_wasm.sh (from anywhere)
set -euo pipefail
source ~/emscripten.sh
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$ROOT/build-wasm"
cd "$ROOT/build-wasm"
deps=()
for d in cadical cadiback cryptominisat5 sbva treedecomp arjun approxmc; do
    deps+=("-D${d}_DIR=$EMINSTALL/lib/cmake/$d")
done
emcmake cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=OFF "${deps[@]}" "$ROOT"
emmake make -j"$(nproc)" ganak-bin bdnnf-bin
cp ganak.js ganak.wasm bdnnf.js bdnnf.wasm "$ROOT/html/"
echo "$(git -C "$ROOT" rev-parse --short HEAD)-$(date +%Y%m%d%H%M%S)" > "$ROOT/html/version.txt"
cd "$ROOT/html" && tsc -p tsconfig.json && tsc -p tsconfig.worker.json
