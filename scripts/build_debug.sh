#!/bin/bash

rm -rf make* CMake* src cmake* ganak* sharp* Make*
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j12 VERBOSE=1
