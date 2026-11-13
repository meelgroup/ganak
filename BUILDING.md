# Building Ganak

Prefer the [release binaries](https://github.com/meelgroup/ganak/releases) or
Nix (`nix shell github:meelgroup/ganak#ganak`). The
[GitHub Action](https://github.com/meelgroup/ganak/actions/workflows/build.yml)
shows the exact steps for Linux, macOS and Windows.

## Dependencies

System libraries: GMP (with C++ support), MPFR, FLINT 3, and optionally zlib
(gzipped CNFs) and zstd (compressed BDNNF; fetched automatically for static
builds). On Debian/Ubuntu:

```shell
sudo apt-get install libgmp-dev libmpfr-dev zlib1g-dev libzstd-dev
```

and FLINT from source if the distribution's is older than 3 (as the CI does):

```shell
wget https://github.com/flintlib/flint/releases/download/v3.2.0-rc1/flint-3.2.0-rc1.tar.gz
tar xzf flint-3.2.0-rc1.tar.gz && cd flint-3.2.0-rc1
./configure --enable-static --enable-shared
make -j$(nproc) && sudo make install
```

The SAT-side dependencies are CryptoMiniSat, Arjun, ApproxMC and treedecomp,
which bring their own (SBVA, CaDiCaL, CaDiBack, EvalMaxSAT). CMake finds each
one in this order:

1. a pre-built tree given as `-D<name>_DIR=...` (e.g. `-Darjun_DIR=../arjun/build`);
2. a sibling checkout next to `ganak/` (e.g. `../arjun`), built as part of Ganak;
3. otherwise it is fetched from GitHub.

So a plain build needs only the system libraries:

```shell
git clone https://github.com/meelgroup/ganak
cd ganak
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Dynamic build against pre-built dependencies

For development, keep all repositories side by side and build each one in its
`build/` directory with its `build_norm.sh`. Then `ganak/build/build_norm.sh`
points CMake at them (`-D<name>_DIR`), builds, and runs the tests.
`scripts/rebuild_dynamic_all.sh` rebuilds the whole chain in order: cadical,
cadiback, breakid, cryptominisat, sbva, EvalMaxSAT, treedecomp, arjun,
approxmc, ganak.

## Static build

First build GMP with position-independent code. GMP's hand-optimised assembly
is normally compiled without `-fPIC` (fine for a native static binary), but
`-fPIC` is required whenever the static `.a` is linked into a shared object,
e.g. the Python extension. `--with-pic` makes both work:

```shell
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --enable-static --enable-cxx --enable-shared --with-pic
make -j$(nproc)
sudo make install
cd ..
```

Then point CMake at the installed static GMP (use `/usr/local/lib/`, not a
custom build directory, which may be compiled for the wrong architecture):

```shell
mkdir build && cd build
cmake -DBUILD_SHARED_LIBS=OFF \
    -DGMPXX_LIBRARY=/usr/local/lib/libgmpxx.a \
    -DGMP_INCLUDE_DIR=/usr/local/include \
    ..
make -j$(nproc)
```

`scripts/build_static.sh` does this against side-by-side pre-built
dependencies (which must themselves be static), and
`scripts/rebuild_static_all.sh` rebuilds the whole chain statically.

## Python package

`pip install pyganak` installs the released wheels. To build the extension from
a checkout (scikit-build-core, `-DBUILD_PYTHON_EXTENSION=ON` is set by
`pyproject.toml`):

```shell
pip install .
```

## WebAssembly

The in-browser version is built with Emscripten. GMP, MPFR and FLINT must
first be built for wasm and installed into `$EMINSTALL` (Emscripten's sysroot);
`.github/workflows/build-emscripten.yml` has the exact steps. Then
`scripts/rebuild_emcripten_all.sh` builds every dependency and Ganak for wasm,
and `scripts/build_emscripten.sh` copies `ganak.js`/`ganak.wasm` and
`bdnnf.js`/`bdnnf.wasm` into `html/` and compiles the page's TypeScript
(`html/src/`, needs `tsc`). Serve `html/` over HTTP (e.g.
`python3 -m http.server -d html`): the page runs each tool in a web worker.

## Useful CMake options

| Option | Effect |
|---|---|
| `-DCMAKE_BUILD_TYPE=Debug` / `Release` / `RelWithDebInfo` | build type |
| `-DBUILD_SHARED_LIBS=OFF` | static binary |
| `-DENABLE_TESTING=ON` | build the tests (`ctest`) |
| `-DSANITIZE=ON` | Clang sanitizers |
