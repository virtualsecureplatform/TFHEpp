option(BENCHMARK_ENABLE_GTEST_TESTS "" OFF)
option(CMAKE_BUILD_TYPE "" Release)

add_subdirectory(google-benchmark)
add_executable(bench bench.cpp)
target_include_directories(
  bench PUBLIC ${PROJECT_SOURCE_DIR}/thirdparties/cereal/include)
target_link_libraries(bench benchmark::benchmark tfhe++)

add_executable(bench_fft bench_fft.cpp)
target_include_directories(
  bench_fft PUBLIC ${PROJECT_SOURCE_DIR}/thirdparties/cereal/include)
target_link_libraries(bench_fft benchmark::benchmark tfhe++)

add_executable(bench_batchboot_blockbinary batchboot_blockbinary.cpp)
target_include_directories(
  bench_batchboot_blockbinary PUBLIC ${PROJECT_SOURCE_DIR}/thirdparties/cereal/include)
target_link_libraries(bench_batchboot_blockbinary benchmark::benchmark tfhe++)

set(BATCHBOOT_BENCHMARK_SLOTS 8 CACHE STRING
    "Packed slot count for bench_batchboot_blockbinary (power of two <= 2048)")
set(BATCHBOOT_BENCHMARK_REPETITIONS 7 CACHE STRING
    "Timed repetitions for bench_batchboot_blockbinary")
target_compile_definitions(bench_batchboot_blockbinary PRIVATE
  BATCHBOOT_BENCHMARK_SLOTS=${BATCHBOOT_BENCHMARK_SLOTS}
  BATCHBOOT_BENCHMARK_REPETITIONS=${BATCHBOOT_BENCHMARK_REPETITIONS})
