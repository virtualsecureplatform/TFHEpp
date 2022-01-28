set(SRCS_FFTWPROC fft_processor_fftw.cpp)

set(FFTWPROC_HEADERS fft_processor_fftw.h)

add_library(fftwproc STATIC ${SRCS_FFTWPROC} ${FFTWPROC_HEADERS})

if(USE_MKL)
  target_include_directories(fftwproc PUBLIC ${PROJECT_SOURCE_DIR}/include)
else()
  target_include_directories(fftwproc PUBLIC ${PROJECT_SOURCE_DIR}/include)
endif()
