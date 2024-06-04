rm -rf CMake* src cmake* ganak* sharp* Make*
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make -j12
