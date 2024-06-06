rm -rf CMake* src cmake* ganak* sharp* Make*
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make -j12 VERBOSE=1
