#!/bin/bash

rm -rf CMake* src cmake* ganak* sharp* Make* -DMAKE_BUILD_TYPE=Debug
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make -j4
