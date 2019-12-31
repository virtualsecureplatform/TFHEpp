#!/bin/bash

clang-format -style=file -i ./src/*.cpp
clang-format -style=file -i ./include/*.hpp
clang-format -style=file -i ./test/*.cpp