#!/bin/bash
set -x
set -e
cd ../../

# export CC=gcc-14
# export CXX=g++-14

cd cadical/
make clean
CXXFLAGS=-fPIC ./configure --competition
make -j12
cd ..

cd cadiback/
make clean
./configure
make -j12
cd ..

cd breakid/build/
./build_static.sh
cd ../../

cd cryptominisat/build/
./build_static.sh
cd ../../

cd sbva/build/
./build_static.sh
cd ../../

cd EvalMaxSAT/build/
./build_static.sh
cd ../../

cd arjun/build/
./build_static.sh
cd ../../

cd approxmc/build/
rm -f build*.sh
ln -s ../scripts/build_scripts/build_*.sh .
./build_static.sh
cd ../../

cd ganak/build/
./build_static.sh
