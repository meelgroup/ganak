rm -rf CMake* src cmake* ganak* sharp* Make*
cmake -DUSE_MOLD=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DSTATICCOMPILE=ON ..
make -j4
strip ganak
