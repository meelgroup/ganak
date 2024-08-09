#!/bin/bash
set -x
set -e
cd ../../

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
git checkout synthesis
./build_static_release.sh
cd ../../

cd sbva/build/
git checkout master
./build_static.sh
cd ../../

cd arjun/build/
git checkout synthesis
./build_static_release.sh
cd ../../

cd approxmc/build/
git checkout synthesis
./build_static_release.sh
cd ../../

cd ganak/build/
# git checkout nodonkey
./build_static_release.sh
