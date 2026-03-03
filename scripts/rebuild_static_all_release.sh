#!/bin/bash
set -x
set -e
cd ../../ || exit 1

cd cadical/ || exit 1
make clean
./build_release.sh
make -j12
cd .. || exit 1

cd cadiback/ || exit 1
make clean
./build_release.sh
make -j12
cd .. || exit 1

cd breakid/build/ || exit 1
./build_static.sh
cd ../../ || exit 1

cd cryptominisat/build/ || exit 1
./build_static_release.sh
cd ../../ || exit 1

cd sbva/build/ || exit 1
./build_static.sh
cd ../../ || exit 1

cd EvalMaxSAT/build/ || exit 1
./build_static_release.sh
cd ../../ || exit 1

cd treedecomp/build/ || exit 1
./build_static_release.sh
cd ../../ || exit 1

cd arjun/build/ || exit 1
./build_static_release.sh
cd ../../ || exit 1

cd approxmc/build/ || exit 1
rm -f build_*.sh
ln -s ../scripts/build_scripts/build_*.sh .
./build_static_release.sh
cd ../../ || exit 1

cd ganak/build/ || exit 1
./build_static_release.sh
