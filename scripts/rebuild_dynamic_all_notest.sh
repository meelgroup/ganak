#!/bin/bash
set -euxo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
echo "Rebuilding all in $BASE_DIR"

cd "$BASE_DIR/cadical/build" || exit 1
./build_notest.sh

cd "$BASE_DIR/cadiback/build" || exit 1
./build_notest.sh

cd "$BASE_DIR/breakid/build" || exit 1
./build_notest.sh

cd "$BASE_DIR/cryptominisat/build" || exit 1
./build_notest.sh

cd "$BASE_DIR/sbva/build" || exit 1
./build_notest.sh

cd "$BASE_DIR/EvalMaxSAT/build" || exit 1
./build_notest.sh

cd "$BASE_DIR/treedecomp/build" || exit 1
./build_notest.sh

cd "$BASE_DIR/arjun/build" || exit 1
./build_notest.sh

cd "$BASE_DIR/approxmc/build" || exit 1
./build_notest.sh

cd "$BASE_DIR/ganak/build" || exit 1
./build_notest.sh
