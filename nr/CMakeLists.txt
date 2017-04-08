cmake_minimum_required(VERSION 3.6)

project(nr)

set(NR_FILES four1.f gasdev.f nr_module.f90 ran1.f)

add_library(nr STATIC ${NR_FILES})
