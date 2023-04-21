rm -rf CMake* src cmake* ganak* sharp* Make*
CXX=clang++ cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make -j4
