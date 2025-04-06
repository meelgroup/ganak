#!/bin/bash
set -x
cd ../../

cd cadical/
git checkout cav2025
make -i clean
CXXFLAGS=-fPIC ./configure --competition
make -j12
cd ..

cd cadiback/
git checkout cav2025
make clean
./configure
make -j12
cd ..

cd breakid/build/
git checkout cav2025
./build_norm.sh
cd ../../

cd cryptominisat/build/
git checkout cav2025
./build_norm.sh
cd ../../

cd sbva/build/
git checkout cav2025
./build_norm.sh
cd ../../

cd arjun/build/
git checkout cav2025
./build_norm.sh
cd ../../

cd approxmc/build/
git checkout cav2025
rm -f build_*.sh
ln -s ../scripts/build_scripts/build_*.sh .
./build_norm.sh
cd ../../

cd ganak/build/
./build_norm.sh
