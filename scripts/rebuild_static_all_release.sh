#!/bin/bash
cd ../../breakid/build/
./build_static.sh
cd ../../cryptominisat/build/
./build_static.sh
cd ../../sbva/build/
./build_static.sh
cd ../../arjun/build/
./build_static_release.sh
cd ../../approxmc/build/
./build_static_release.sh
cd ../../ganak/build/
./build_static_release.sh
