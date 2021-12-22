option(BENCHMARK_ENABLE_GTEST_TESTS "" OFF)
option(CMAKE_BUILD_TYPE "" Release)

add_subdirectory(google-benchmark)
add_executable(bench bench.cpp)
target_include_directories(
  bench PUBLIC ${PROJECT_SOURCE_DIR}/thirdparties/cereal/include)
target_link_libraries(bench benchmark::benchmark tfhe++)
