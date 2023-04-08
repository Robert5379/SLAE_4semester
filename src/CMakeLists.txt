cmake_minimum_required(VERSION 3.24)

file(GLOB_RECURSE COPY_DEMO_SRC *.cpp *.h)
add_library(slae_demo ${COPY_DEMO_SRC})
