rm -rf CMake* src cmake* ganak* sharp* Make*
cmake -DCMAKE_BUILD_TYPE=Release -DSTATICCOMPILE=ON ..
make -j14
strip ganak
