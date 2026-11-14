#!/bin/bash
# Builds ganak and its dependencies for wasm, each in its own build-wasm/ (the
# native build/ dirs stay), and puts ganak.js/bdnnf.js into ganak/html/.
set -euxo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
echo "Rebuilding all in $BASE_DIR for emscripten"
source ~/emscripten.sh

"$BASE_DIR/cadical/build_emscripten.sh"
"$BASE_DIR/cadiback/build_emscripten.sh"
"$BASE_DIR/cryptominisat/scripts/build_scripts/build_emscripten.sh"
"$BASE_DIR/sbva/scripts/build_emscripten.sh"
"$BASE_DIR/treedecomp/scripts/build_emscripten.sh"
"$BASE_DIR/arjun/scripts/build_emscripten.sh"
"$BASE_DIR/approxmc/scripts/build_emscripten.sh"
"$BASE_DIR/ganak/scripts/build_emscripten.sh"
