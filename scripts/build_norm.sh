rm -rf CMake* src cmake* ganak* sharp* Make* Testing* tests* include
cmake -DENABLE_TESTING=ON ..
make -j12
make test
