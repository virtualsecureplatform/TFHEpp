cmake_minimum_required(VERSION 3.9)
project(tfhe++ CXX)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_FLAGS "-march=native -O3 -std=c++17")

# For OpenMP
find_package(OpenMP)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(thirdparties/randen)
add_subdirectory(thirdparties/spqlios)
add_subdirectory(src)
add_subdirectory(test)
