rm -rf CMake* cmake sharpSAT Makefile cmake_install.cmake
cmake -DSTATICCOMPILE=ON ..
make -j4
strip ganak
