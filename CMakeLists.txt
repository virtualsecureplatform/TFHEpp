cmake_minimum_required(VERSION 3.16)
project(tfhe++ CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-march=native -O3 -g -funroll-loops")

OPTION(USE_80BIT_SECURITY "Use 80bit security parameter(faster)" OFF)
OPTION(USE_CGGI19 "Use the parameter set proposed in CGGI19" OFF)
OPTION(ENABLE_TEST "Build tests" OFF)
OPTION(ENABLE_BENCHMARK "Build benchmarks" OFF)
OPTION(ENABLE_TUTORIAL "Build tutorial" OFF)

if(USE_80BIT_SECURITY)
    add_compile_definitions(USE_80BIT_SECURITY)
elseif(USE_CGGI19)
    add_compile_definitions(USE_CGGI19)
endif()

# For OpenMP
find_package(OpenMP)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(thirdparties/randen)
add_subdirectory(thirdparties/spqlios)
add_subdirectory(src)

if(ENABLE_TEST)
    add_subdirectory(test)
endif()
if(ENABLE_BENCHMARK)
    add_subdirectory(benchmark)
endif()
if(ENABLE_TUTORIAL)
    add_subdirectory(tutorial)
endif()

install(TARGETS tfhe++ LIBRARY DESTINATION lib)
