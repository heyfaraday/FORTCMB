cmake_minimum_required(VERSION 3.6)

project(lib)

set(CMAKE_Fortran_COMPILER_ID Intel)

include_directories(/usr/local/Healpix/include)

link_directories(/usr/local/Healpix/lib)
link_directories(/usr/local/Cellar/fftw)
link_directories(/usr/local/Cellar/cfitsio)

set(CMAKE_Fortran_FLAGS
"${CMAKE_Fortran_FLAGS} -warn all -std=f2008 -clang-name=clang")

set(CMAKE_Fortran_FLAGS_DEBUG
"${CMAKE_Fortran_FLAGS} -check pointer -debug-parameters all -g -traceback")

set(CMAKE_Fortran_FLAGS_RELEASE
"${CMAKE_Fortran_FLAGS} -O3 -m64 -ip -xCORE-AVX2")
# or -xHOST

set(CMAKE_EXE_LINKER_FLAGS "-lhealpix -lgif -lcfitsio -lfftw3")
set(STATIC_LIBRARY_FLAGS "-lhealpix -lgif -lcfitsio -lfftw3")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

set(LIB_FILES distance.f90 map_coef.f90 minkowski.f90 write_in_file.f90)

add_library(lib STATIC ${LIB_FILES})
