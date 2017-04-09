# fortCMB

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)]()
[![Status](https://img.shields.io/badge/status-dev-ff69b4.svg)]()

Library for statistical testing of the Cosmic Microwave Background radiation.

## Requirements
* [cmake](https://cmake.org)
* [ifort](https://software.intel.com/en-us/intel-compilers)
* [healpix](http://healpix.sourceforge.net)
* [fftw](https://github.com/FFTW/fftw3)
* [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
* [python](https://www.python.org)

## CMake
Debug
```
cmake .. -DCMAKE_BUILD_TYPE=Debug
/ -warn all -check pointer -debug all -g -traceback /
```
and Release options
```
cmake .. -DCMAKE_BUILD_TYPE=Release /
/ -warn all -O3 -m64 -ip -xCORE-AVX2 /
```
