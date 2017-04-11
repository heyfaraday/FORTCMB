cmake_minimum_required(VERSION 3.6)

project(bin)

set(CMAKE_Fortran_COMPILER_ID Intel)

include_directories(/usr/local/Healpix/include)

link_directories(/usr/local/Healpix/lib)
link_directories(/usr/local/Cellar/fftw)
link_directories(/usr/local/Cellar/cfitsio)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -std=f2008")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -check pointer -debug all -g -traceback")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -m64 -ip -xCORE-AVX2")
# or -xHOST

set(CMAKE_EXE_LINKER_FLAGS "-lhealpix -lgif -lcfitsio -lfftw3")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(BIN_FILES map_test.f90)
add_executable(map_test.a ${BIN_FILES})

target_link_libraries(map_test.a lib nr)
