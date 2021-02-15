set -x
set -e

rm -rf CMake* cmake_install.cmake Makefile src sym_ganak
cmake -DSTATICCOMPILE=ON ..
make -j4 VERBOSE=1
