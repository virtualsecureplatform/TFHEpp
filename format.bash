#!/bin/bash

clang-format -style=file -i ./src/*.cpp
clang-format -style=file -i ./src/axell/*.cpp
clang-format -style=file -i ./thirdparties/fftw/*.h
clang-format -style=file -i ./thirdparties/fftw/*.cpp
clang-format -style=file -i ./thirdparties/spqliox_aarch64/*.h
clang-format -style=file -i ./thirdparties/spqliox_aarch64/*.cpp
clang-format -style=file -i ./include/*.hpp
clang-format -style=file -i ./include/params/*.hpp
clang-format -style=file -i ./include/axell/*.hpp
clang-format -style=file -i ./test/*.cpp
clang-format -style=file -i ./benchmark/*.cpp
clang-format -style=file -i ./tutorial/*.cpp

cmake-format -i CMakeLists.txt
cmake-format -i src/CMakeLists.txt
cmake-format -i benchmark/CMakeLists.txt
cmake-format -i thirdparties/spqlios/CMakeLists.txt
cmake-format -i thirdparties/fftw/CMakeLists.txt
cmake-format -i thirdparties/spqliox_aarch64/CMakeLists.txt