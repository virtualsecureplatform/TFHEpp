file(GLOB tfhe++_sources "*.cpp")

add_library(tfhe++ ${tfhe++_sources})
target_include_directories(
  tfhe++
  PUBLIC ${PROJECT_SOURCE_DIR}/include
         ${PROJECT_SOURCE_DIR}/thirdparties/spqlios
         ${PROJECT_SOURCE_DIR}/thirdparties/fftw
         ${PROJECT_SOURCE_DIR}/thirdparties/randen
         ${PROJECT_SOURCE_DIR}/thirdparties/cereal/include)

if(USE_FFTW3)
  target_link_libraries(tfhe++ INTERFACE fftwproc randen fftw3)
else()
  target_link_libraries(tfhe++ INTERFACE spqlios randen)
endif()
