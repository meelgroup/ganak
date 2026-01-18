#!/bin/bash
set -x
set -e
cd ../../

cd cadical/
make clean
./build_release.sh
make -j12
cd ..

cd cadiback/
make clean
./build_release.sh
make -j12
cd ..

cd breakid/build/
./build_static.sh
cd ../../

cd cryptominisat/build/
./build_static_release.sh
cd ../../

cd sbva/build/
./build_static.sh
cd ../../

cd arjun/build/
./build_static_release.sh
cd ../../

cd approxmc/build/
rm -f build_*.sh
ln -s ../scripts/build_scripts/build_*.sh .
./build_static_release.sh
cd ../../

cd ganak/build/
./build_static_release.sh
