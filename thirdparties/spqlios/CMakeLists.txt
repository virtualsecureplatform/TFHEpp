# This is the nayuki fft processor for the tfhe library
enable_language(ASM)

if(USE_AVX512)
  # The dedicated AVX-512 SPQLIOS FFT kernel is currently numerically
  # unstable on the TFHEpp test suite. Keep USE_AVX512 enabled for the
  # higher-level FD arithmetic intrinsics, but route the FFT backend through
  # the stable FMA implementation until the AVX-512 kernel is fixed.
  set(SRCS_FMA spqlios-fft-fma.s spqlios-ifft-fma.s fft_processor_spqlios.cpp
               spqlios-fft-impl.cpp)
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
