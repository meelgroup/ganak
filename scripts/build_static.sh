rm -rf CMake* src cmake* ganak* sharp* Make*
cmake -DSTATICCOMPILE=ON ..
make -j$(nproc)
strip ganak
