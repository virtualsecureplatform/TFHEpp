set(SPQLIOX_AARCH64_SRCS fft_processor_spqliox_aarch64.cpp)
set(SPQLIOX_AARCH64_HEADERS fft_processor_spqliox_aarch64.h
                            spqliox_aarch64_impl.h)

add_library(spqliox_aarch64 STATIC ${SPQLIOX_AARCH64_SRCS}
                                   ${SPQLIOX_AARCH64_HEADERS})

target_include_directories(spqliox_aarch64 PUBLIC ${PROJECT_SOURCE_DIR}/include)
