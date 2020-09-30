cmake_minimum_required(VERSION 3.9)
project(tfhe++ CXX)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_FLAGS "-march=native -O3 -std=c++17")

OPTION(USE_80BIT_SECURITY "Use 80bit security parameter(faster)" ON)
if(USE_80BIT_SECURITY)
    add_compile_definitions(USE_80BIT_SECURITY)
endif()

# For OpenMP
find_package(OpenMP)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(thirdparties/randen)
add_subdirectory(thirdparties/spqlios)
add_subdirectory(src)
add_subdirectory(test)

install(TARGETS tfhe++ LIBRARY DESTINATION lib)