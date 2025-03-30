#!/bin/bash
set -e

rm -rf lib* Test* tests* include tests CM* cmake* arjun Makefile rjun-src
emcmake cmake -DCMAKE_INSTALL_PREFIX=$EMINSTALL -DENABLE_TESTING=OFF -DPOLYS=OFF ..
emmake make -j26
emmake make install
cp ganak.wasm ../html
cp $EMINSTALL/bin/ganak.js ../html
