#!/bin/bash

clang-format -style=file -i ./src/*.cpp
clang-format -style=file -i ./include/*.hpp
clang-format -style=file -i ./include/params/*.hpp
clang-format -style=file -i ./test/*.cpp
clang-format -style=file -i ./benchmark/*.cpp
clang-format -style=file -i ./tutorial/*.cpp