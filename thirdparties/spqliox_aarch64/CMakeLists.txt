set(SPQLIOX_AARCH64_SRCS fft_processor_spqliox_aarch64.cpp)
set(SPQLIOX_AARCH64_HEADERS fft_processor_spqliox_aarch64.h
                            spqliox_aarch64_impl.h)

add_library(spqliox_aarch64 STATIC ${SPQLIOX_AARCH64_SRCS}
                                   ${SPQLIOX_AARCH64_HEADERS})

# Check if the platform is macOS and the architecture is ARM64
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
  # Add the include directory provided by Homebrew
  target_include_directories(
    spqliox_aarch64 PUBLIC ${PROJECT_SOURCE_DIR}/include /opt/homebrew/include)
else()
  target_include_directories(spqliox_aarch64
                             PUBLIC ${PROJECT_SOURCE_DIR}/include)
endif()
