rm -rf CMake* src cmake* ganak* sharp* Make*
CXX=clang++ cmake -DUSE_MOLD=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make -j4
