#!/bin/bash
set -e

rm -rf lib* Test* tests* include tests CM* cmake* arjun Makefile rjun-src
emcmake cmake -DCMAKE_INSTALL_PREFIX=$EMINSTALL -DPOLYS=OFF ..
emmake make -j26
emmake make install
