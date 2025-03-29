rm -rf CMake* src cmake* ganak* sharp* Make* Testing* tests* include
emcmake cmake -DENABLE_TESTING=OFF -DPOLYS=OFF -DCMAKE_INSTALL_PREFIX=$EMINSTALL ..
emmake make -j12
emmake make install
