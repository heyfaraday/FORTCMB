cmake_minimum_required(VERSION 3.6)

project(nr)

set(NR_FILES four1.f90 fourrow.f90 gasdev.f90 nr.f90 nrtype.f90 nrutil.f90
    ran_state.f90 ran1.f90)

add_library(nr STATIC ${NR_FILES})
