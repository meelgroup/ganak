rm -rf CMake* src cmake* ganak* sharp* Make*
cmake -DANALYZE=ON -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang ..
make -j16
