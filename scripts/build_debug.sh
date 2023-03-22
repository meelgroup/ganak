#!/bin/bash

rm -rf CMake* src cmake* ganak* sharp* Make* -DMAKE_BUILD_TYPE=Debug
cmake ..
make -j4
