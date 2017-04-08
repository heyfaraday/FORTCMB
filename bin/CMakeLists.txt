cmake_minimum_required(VERSION 3.6)

project(bin)

set(BIN_FILES allocate_test.f90)
add_executable(allocate_test.a ${BIN_FILES})

target_link_libraries(allocate_test.a lib nr)
