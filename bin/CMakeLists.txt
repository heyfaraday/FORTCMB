cmake_minimum_required(VERSION 3.6)

project(bin)

set(CMAKE_Fortran_COMPILER_ID Intel)

include_directories(/usr/local/Healpix/include)

link_directories(/usr/local/Healpix/lib)
link_directories(/usr/local/Cellar/fftw)
link_directories(/usr/local/Cellar/cfitsio)

set(CMAKE_Fortran_FLAGS "-warn all -check all -O2 -m64 -debug all -check all -g -traceback")
set(CMAKE_EXE_LINKER_FLAGS "-lhealpix -lgif -lcfitsio -lfftw3")

set(BIN_FILES hp_test.f90)
add_executable(hp_test.a ${BIN_FILES})

target_link_libraries(hp_test.a lib nr)
