rm -rf CMake* cmake sharpSAT
cmake -DSTATICCOMPILE=ON ..
make -j4
