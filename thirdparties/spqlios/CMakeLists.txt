# This is the nayuki fft processor for the tfhe library
enable_language(ASM)

if(USE_AVX512)
  set(SRCS_FMA spqlios-fft-avx512.s spqlios-ifft-avx512.s
               fft_processor_spqlios.cpp spqlios-fft-impl-avx512.cpp)
else()
  set(SRCS_FMA spqlios-fft-fma.s spqlios-ifft-fma.s fft_processor_spqlios.cpp
               spqlios-fft-impl.cpp)
endif()

set(SPQLIOS_HEADERS fft_processor_spqlios.h x86.h)

if(ENABLE_SHARED)
  add_library(spqlios SHARED ${SRCS_FMA} ${SPQLIOS_HEADERS})
else()
  add_library(spqlios STATIC ${SRCS_FMA} ${SPQLIOS_HEADERS})
endif()

target_include_directories(spqlios PUBLIC ${PROJECT_SOURCE_DIR}/include)
