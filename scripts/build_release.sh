rm -rf CMake* src cmake* ganak* sharp* Make*
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j12 VERBOSE=1
