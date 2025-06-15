#!/bin/bash
set -x
set -e
cd ../../

# export CC=gcc-14
# export CXX=g++-14

cd cadical/
git checkout add_dynamic_lib
make clean
CXXFLAGS=-fPIC ./configure --competition
make -j12
cd ..

cd cadiback/
git checkout synthesis
make clean
./configure
make -j12
cd ..

cd breakid/build/
git checkout master
./build_static.sh
cd ../../

cd cryptominisat/build/
git checkout master
./build_static.sh
cd ../../

cd sbva/build/
git checkout master
./build_static.sh
cd ../../

cd arjun/build/
git checkout master
./build_static.sh
cd ../../

cd approxmc/build/
git checkout master
rm -f build*.sh
ln -s ../scripts/build_scripts/build_*.sh .
./build_static.sh
cd ../../

cd ganak/build/
./build_static.sh
