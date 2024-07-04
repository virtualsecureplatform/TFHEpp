set(SRCS_MKLPROC fft_processor_mkl.cpp)

set(MKLPROC_HEADERS fft_processor_mkl.hpp)

add_library(mklproc STATIC ${SRCS_MKLPROC} ${MKLPROC_HEADERS})

target_include_directories(mklproc PUBLIC ${PROJECT_SOURCE_DIR}/include)
