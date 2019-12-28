cmake_minimum_required(VERSION 3.1)
project(tfhe++ CXX)
set(CMAKE_CXX_COMPILER clang++)

add_subdirectory(src)
add_subdirectory(test)