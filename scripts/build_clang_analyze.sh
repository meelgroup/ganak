rm -rf CMake*
rm -rf src
rm -rf cmake*
rm -rf ganak*
rm -rf sharp*
rm -rf Make*
rm -rf deps
rm -rf _deps
SAT_DIR=$(cd ../.. && pwd)
cmake -DANALYZE=ON -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang \
    -Dcryptominisat5_DIR="${SAT_DIR}/cryptominisat/build" \
    -Dsbva_DIR="${SAT_DIR}/sbva/build" \
    -Dtreedecomp_DIR="${SAT_DIR}/treedecomp/build" \
    -DEvalMaxSAT_DIR="${SAT_DIR}/EvalMaxSAT/build" \
    -Darjun_DIR="${SAT_DIR}/arjun/build" \
    -Dapproxmc_DIR="${SAT_DIR}/approxmc/build" \
    -Dcadical_DIR="${SAT_DIR}/cadical/build" \
    -Dcadiback_DIR="${SAT_DIR}/cadiback/build" \
    ..
make -j$(nproc)
