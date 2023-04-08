cmake_minimum_required(VERSION 3.24)
project(SLAE)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS  "-g")
add_subdirectory(src)

add_subdirectory(TESTS)
