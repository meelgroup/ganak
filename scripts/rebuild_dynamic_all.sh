#!/bin/bash
set -x
set -e
cd ../../

cd cadical/
# git checkout add_dynamic_lib
make clean
CXXFLAGS=-fPIC ./configure --competition
make -j12
cd ..

cd cadiback/
# git checkout synthesis
make clean
./configure
make -j12
cd ..

cd breakid/build/
# git checkout master
./build_norm.sh
cd ../../

cd cryptominisat/build/
# git checkout synthesis
./build_norm.sh
cd ../../

cd sbva/build/
# git checkout master
./build_norm.sh
cd ../../

cd arjun/build/
# git checkout synthesis
./build_norm.sh
cd ../../

cd approxmc/build/
# git checkout synthesis
./build_norm.sh
cd ../../

cd ganak/build/
# git checkout nodonkey
./build_norm.sh
