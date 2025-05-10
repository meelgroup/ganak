# Building from source

```bash
sudo apt-get install -yq libgmp-dev libmpfr-dev

git clone https://github.com/meelgroup/cadical
cd cadical
git checkout add_dynamic_lib
./configure
make
cd ..

git clone https://github.com/meelgroup/cadiback
cd cadiback
git checkout synthesis
./configure
make
cd ..

git clone https://github.com/msoos/cryptominisat
cd cryptominisat
git checkout synthesis
mkdir build && cd build
cmake ..
make
cd ../..

git clone https://github.com/meelgroup/sbva
cd sbva
mkdir build && cd build
cmake ..
make
cd ../..

git clone https://github.com/meelgroup/breakid
cd breakid
mkdir build && cd build
cmake ..
make
cd ../..

git clone https://github.com/meelgroup/approxmc
cd approxmc
mkdir build && cd build
git checkout synthesis
cmake ..
make
cd ../..

git clone https://github.com/meelgroup/ganak
cd ganak
mkdir build && cd build
cmake ..
make
```

if you also need to have the system run on older or different machines, you
need to have it compile for Sandybridge architecture. To do this, you need to
 recompile GMP and MPFR with the right flags. Here is how to do it:

```bash
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar xvf gmp-6.3.0.tar.xz
cd gmp-6.3.0
CFLAGS="-march=sandybridge" ./configure --enable-cxx
make -j12
sudo make install
cd ..
wget https://ftp.gnu.org/gnu/mpfr/mpfr-4.1.0.tar.xz
tar xvf mpfr-4.1.0.tar.xz
cd mpfr-4.1.0
CFLAGS="-march=sandybridge" ./configure
make -j12
sudo make install
cd ..
```
