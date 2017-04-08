cmake_minimum_required(VERSION 3.6)

project(lib)

set(LIB_FILES write_in_file.f90)

add_library(lib STATIC ${LIB_FILES})
