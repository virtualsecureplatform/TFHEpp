OPTION(BENCHMARK_ENABLE_GTEST_TESTS "" OFF)
OPTION(CMAKE_BUILD_TYPE "" Release)

add_subdirectory(google-benchmark)
add_executable(bench_nand bench_nand.cpp)
target_include_directories(bench_nand
  PUBLIC ${PROJECT_SOURCE_DIR}/thirdparties/cereal/include
  )
target_link_libraries(bench_nand benchmark::benchmark tfhe++)