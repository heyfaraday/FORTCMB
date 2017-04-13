cmake_minimum_required(VERSION 3.6)

project(nr)

set(CMAKE_Fortran_COMPILER_ID Intel)

set(CMAKE_Fortran_FLAGS
"${CMAKE_Fortran_FLAGS} -warn all -std=f2008 -clang-name=clang")

set(CMAKE_Fortran_FLAGS_DEBUG
"${CMAKE_Fortran_FLAGS} -check pointer -debug-parameters all -g -traceback")

set(CMAKE_Fortran_FLAGS_RELEASE
"${CMAKE_Fortran_FLAGS} -O3 -m64 -ip -xCORE-AVX2")
# or -xHOST

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(NR_FILES four1.f90 fourrow.f90 gasdev.f90 nr.f90 nrtype.f90 nrutil.f90
    ran_state.f90 ran1.f90)

add_library(nr STATIC ${NR_FILES})
