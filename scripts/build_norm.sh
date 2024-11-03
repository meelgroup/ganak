rm -rf CMake* src cmake* ganak* sharp* Make* Testing* tests*
cmake -DUSE_MOLD=ON -DENABLE_TESTING=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make -j12
make test
