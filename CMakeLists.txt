cmake_minimum_required(VERSION 3.1)
project(tfhe++ CXX)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_FLAGS "-march=native -O3 -g -std=c++17 -lomp -lprofiler -ltcmalloc")

set(USE_80BIT_SECURITY OFF BOOL "Use old 80bit security parameters. It is faster, but not recommended.")

add_subdirectory(thirdparties/randen)
add_subdirectory(thirdparties/spqlios)
add_subdirectory(src)
add_subdirectory(test)
