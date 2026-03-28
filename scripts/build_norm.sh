#!/usr/bin/env bash
set -euo pipefail

rm -rf CMake* src cmake* ganak* sharp* Make* Testing* tests* include .cmake
cmake -DENABLE_TESTING=ON ..
make -j$(nproc)
make test
