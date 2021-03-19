cmake_minimum_required(VERSION 3.9)
project(tfhe++ CXX)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS "-march=native -O3 -std=c++17 -g -funroll-loops")

OPTION(USE_80BIT_SECURITY "Use 80bit security parameter(faster)" OFF)
OPTION(ENABLE_TEST "Build tests" OFF)
OPTION(ENABLE_TUTORIAL "Build tutorial" OFF)

if(USE_80BIT_SECURITY)
    add_compile_definitions(USE_80BIT_SECURITY)
elseif(USE_CGGI16)
    add_compile_definitions(USE_CGGI16)
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
if(ENABLE_TUTORIAL)
    add_subdirectory(tutorial)
endif()

install(TARGETS tfhe++ LIBRARY DESTINATION lib)