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
./configure --release
make -j12
cd ..

cd breakid/build/
git checkout master
./build_static.sh
cd ../../

cd cryptominisat/build/
git checkout master
./build_static_release.sh
cd ../../

cd sbva/build/
git checkout master
./build_static.sh
cd ../../

cd arjun/build/
git checkout master
./build_static_release.sh
cd ../../

cd approxmc/build/
git checkout master
rm -f build_*.sh
ln -s ../scripts/build_scripts/build_*.sh .
./build_static_release.sh
cd ../../

cd ganak/build/
./build_static_release.sh
