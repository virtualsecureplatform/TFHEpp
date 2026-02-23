cmake_minimum_required(VERSION 3.16)
project(tfhe++ CXX C)
set(CMAKE_CXX_STANDARD 20)
option(DEBUG "Debug mode" OFF)
if(DEBUG)
  set(CMAKE_CXX_FLAGS
      "-march=native -g -Wall -Wextra -pedantic -Wno-sign-compare")
else()
  set(CMAKE_CXX_FLAGS
      "-march=native -O3 -g -funroll-loops -Wall -Wextra -pedantic -Wno-sign-compare"
  )
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

option(USE_80BIT_SECURITY "Use 80bit security parameter(faster)" OFF)
option(USE_CGGI19 "Use the parameter set proposed in CGGI19" OFF)
option(USE_CONCRETE "Use the parameter set proposed in CONCRETE" OFF)
option(USE_TFHE_RS "Use the parameter set proposed in TFHE-rs" OFF)
option(USE_TERNARY "Use ternary secret keys" OFF)
option(USE_COMPRESS "Use compressed ciphertexts" OFF)
option(USE_SUBSET_KEY "Use subset key" OFF)

option(USE_TERNARY_CMUX "Use ternary cmux" OFF)
option(USE_KEY_BUNDLE "Use key bundle algorithm" OFF)

option(USE_BLAKE3 "Use blake3 as CSPRNG" ON)
option(USE_RANDEN "Use randen as CSPRNG" OFF)

# Detect AVX-512 support at configure time
if(NOT DEFINED AVX512_SUPPORTED)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
    include(CheckCXXSourceRuns)
    set(CMAKE_REQUIRED_FLAGS "-mavx512f")
    check_cxx_source_runs("
      #include <immintrin.h>
      int main(int argc, char**) {
        double data[8] __attribute__((aligned(64)));
        for (int i = 0; i < 8; i++) data[i] = static_cast<double>(argc + i);
        __m512d a = _mm512_load_pd(data);
        __m512d b = _mm512_add_pd(a, a);
        _mm512_store_pd(data, b);
        return static_cast<int>(data[0]) - 2 * argc;
      }
    " AVX512_SUPPORTED)
    unset(CMAKE_REQUIRED_FLAGS)
  else()
    set(AVX512_SUPPORTED FALSE)
  endif()
endif()

if(AVX512_SUPPORTED AND EXISTS "/proc/cpuinfo")
  file(READ "/proc/cpuinfo" TFHEPP_CPUINFO)
  if(NOT TFHEPP_CPUINFO MATCHES "avx512f")
    message(STATUS "AVX-512 flags not found in /proc/cpuinfo; disabling AVX-512")
    set(AVX512_SUPPORTED FALSE)
  endif()
endif()

if(AVX512_SUPPORTED)
  message(STATUS "AVX-512 support detected, USE_AVX512 defaults to ON")
  set(USE_AVX512_DEFAULT ON)
else()
  message(STATUS "AVX-512 support not detected, USE_AVX512 defaults to OFF")
  set(USE_AVX512_DEFAULT OFF)
endif()

option(USE_AVX512 "Use AVX512 ver. SPQLIOS" ${USE_AVX512_DEFAULT})
option(USE_CONCRETE_FFT "Use concrete-fft" OFF)
option(USE_MKL "Use Intel MKL" OFF)
option(USE_FFTW3 "Use FFTW3" OFF)
option(USE_SPQLIOX_AARCH64 "Use spqliox_aarch64" OFF)
option(USE_SPQLIOS_ARITHMETIC "Use spqlios-arithmetic backend" OFF)
option(USE_HEXL "Use Intel HEXL" OFF)

option(ENABLE_TEST "Build tests" OFF)
option(ENABLE_BENCHMARK "Build benchmarks" OFF)
option(ENABLE_TUTORIAL "Build tutorial" OFF)
option(USE_PERF "Use Google Profiler" OFF)

option(ENABLE_SHARED "Build as shared libraries" OFF)

set(TFHEpp_DEFINITIONS
    ""
    PARENT_SCOPE)

if(USE_BLAKE3)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_BLAKE3"
      PARENT_SCOPE)
  add_compile_definitions(USE_BLAKE3)
  add_subdirectory(thirdparties/BLAKE3/c)
elseif(USE_RANDEN)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_RANDEN"
      PARENT_SCOPE)
  add_compile_definitions(USE_RANDEN)
endif()

if(USE_80BIT_SECURITY)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_80BIT_SECURITY"
      PARENT_SCOPE)
  add_compile_definitions(USE_80BIT_SECURITY)
elseif(USE_COMPRESS)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_COMPRESS"
      PARENT_SCOPE)
  add_compile_definitions(USE_COMPRESS)
elseif(USE_CGGI19)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_CGGI19"
      PARENT_SCOPE)
  add_compile_definitions(USE_CGGI19)
elseif(USE_CONCRETE)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_CONCRETE"
      PARENT_SCOPE)
  add_compile_definitions(USE_CONCRETE)
elseif(USE_TFHE_RS)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_TFHE_RS"
      PARENT_SCOPE)
  add_compile_definitions(USE_TFHE_RS)
elseif(USE_TERNARY_CMUX)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_TERNARY_CMUX;USE_TERNARY"
      PARENT_SCOPE)
  add_compile_definitions(USE_TERNARY)
  add_compile_definitions(USE_TERNARY_CMUX)
elseif(USE_TERNARY)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_TERNARY"
      PARENT_SCOPE)
  add_compile_definitions(USE_TERNARY)
endif()

if(NOT USE_TERNARY)
  if(USE_KEY_BUNDLE)
    set(TFHEpp_DEFINITIONS
        "${TFHEpp_DEFINITIONS};USE_KEY_BUNDLE"
        PARENT_SCOPE)
    add_compile_definitions(USE_KEY_BUNDLE)
  endif()
endif()

# if(USE_AVX512) string(APPEND CMAKE_CXX_FLAGS " -mprefer-vector-width=512")
# endif()

if(USE_FFTW3)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_FFTW3"
      PARENT_SCOPE)
  add_compile_definitions(USE_FFTW3)
  add_subdirectory(thirdparties/fftw)
elseif(USE_MKL)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_MKL;USE_INTERLEAVED_FORMAT"
      PARENT_SCOPE)
  add_compile_definitions(USE_MKL)
  add_compile_definitions(USE_INTERLEAVED_FORMAT)
  find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
  include_directories(${MKLROOT}/include)
  add_subdirectory(thirdparties/mkl)
elseif(USE_SPQLIOX_AARCH64)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_SPQLIOX_AARCH64"
      PARENT_SCOPE)
  add_compile_definitions(USE_SPQLIOX_AARCH64)
  add_subdirectory(thirdparties/spqliox_aarch64)
  # Check if the platform is macOS and the architecture is ARM64
  if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    include_directories(/opt/homebrew/include)
  endif()
elseif(USE_CONCRETE_FFT)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_CONCRETE_FFT;USE_INTERLEAVED_FORMAT"
      PARENT_SCOPE)
  add_compile_definitions(USE_CONCRETE_FFT)
  add_compile_definitions(USE_INTERLEAVED_FORMAT)
  add_subdirectory(thirdparties/concrete-fft)
elseif(USE_SPQLIOS_ARITHMETIC)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_SPQLIOS_ARITHMETIC"
      PARENT_SCOPE)
  add_compile_definitions(USE_SPQLIOS_ARITHMETIC)
  # Disable spqlios-arithmetic tests to avoid gtest dependency
  set(ENABLE_TESTING
      OFF
      CACHE BOOL "Disable spqlios-arithmetic tests" FORCE)
  add_subdirectory(thirdparties/spqlios-arithmetic)
  # Restore ENABLE_TEST for the parent project
  set(ENABLE_TESTING ${ENABLE_TEST})
else()
  if(USE_AVX512)
    set(TFHEpp_DEFINITIONS
        "${TFHEpp_DEFINITIONS};USE_AVX512"
        PARENT_SCOPE)
    add_compile_definitions(USE_AVX512)
  endif()
  add_subdirectory(thirdparties/spqlios)
endif()

if(USE_HEXL)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_HEXL"
      PARENT_SCOPE)
  add_compile_definitions(USE_HEXL)
  add_subdirectory(thirdparties/hexl)
  # set(CMAKE_CXX_FLAGS "-march=native -O3 -g -funroll-loops -Wall -Wextra
  # -pedantic -Wno-sign-compare -mprefer-vector-width=512" )
endif()

# For OpenMP
find_package(OpenMP)
if(OpenMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(USE_RANDEN)
  add_subdirectory(thirdparties/randen)
endif()
add_subdirectory(src)

if(ENABLE_TEST)
  add_subdirectory(test)
  add_subdirectory(thirdparties/AES)
endif()
if(ENABLE_BENCHMARK)
  add_subdirectory(benchmark)
endif()
if(ENABLE_TUTORIAL)
  add_subdirectory(tutorial)
endif()

install(TARGETS tfhe++ LIBRARY DESTINATION lib)
