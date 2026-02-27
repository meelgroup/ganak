# How to build for Emscripten

Install emscripten.sh, then have an emscripten.sh in ~/emscripten.sh with the following content:

```bash
export EMSDK_QUIET=1
source "/home/soos/development/emsdk/emsdk_env.sh"
export EMINSTALL="${HOME}/development/emsdk/upstream/emscripten/cache/sysroot"
```

build MPFI for Emscripten

```bash
emconfigure ./configure --with-mpfr=$EMINSTALL --prefix=$EMINSTALL
emmake make -j 4
emmake make install
```

Useful commands for building stuff;
```
EMCONFIGURE_JS=1 emconfigure ./configure ABI=64 --disable-assembly CFLAGS=-m64 CXXFLAGS=-m64 LDFLAGS=-m64 --enable-cxx
EMCONFIGURE_JS=1 emconfigure ./configure ABI=64 --disable-assembly --enable-cxx
emconfigure ./configure ABI=64 --disable-assembly --enable-cxx
emconfigure ./configure --disable-assembly --host wasm32 --enable-cxx --prefix=$EMINSTALL
emconfigure ./configure --host wasm32 --with-gmp=$EMINSTALL --prefix=$EMINSTALL
emconfigure ./configure --with-mpir=$EMINSTALL --with-mpfr=$EMINSTALL --prefix=$EMINSTALL
```
