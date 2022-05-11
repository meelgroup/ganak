rm -rf CMake* cmake ganak Makefile cmake_install.cmake
cmake -DSTATICCOMPILE=ON ..
make -j4
strip ganak
