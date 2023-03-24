#!/bin/bash

rm -rf CMake* src cmake* ganak* sharp* Make*
cmake -DDOPCC=OFF ..
make -j4
