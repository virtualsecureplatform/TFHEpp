file(GLOB_RECURSE tfhe++_sources "*.cpp")

add_library(tfhe++ ${tfhe++_sources})
target_include_directories(
  tfhe++
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/thirdparties/spqlios
    ${PROJECT_SOURCE_DIR}/thirdparties/fftw
    ${PROJECT_SOURCE_DIR}/thirdparties/spqliox_aarch64
    ${PROJECT_SOURCE_DIR}/thirdparties/spqliox_aarch64/xbyak_aarch64/xbyak_aarch64
    ${PROJECT_SOURCE_DIR}/thirdparties/hexl/hexl/hexl/include
    ${PROJECT_SOURCE_DIR}/thirdparties/randen
    ${PROJECT_SOURCE_DIR}/thirdparties/cereal/include)

if(USE_RANDEN)
  target_link_libraries(tfhe++ INTERFACE randen)
endif()

if(USE_FFTW3)
  target_link_libraries(tfhe++ INTERFACE fftwproc)
  if(USE_MKL)
    target_include_directories(tfhe++ PUBLIC ${MKLROOT}/include/fftw)
    target_link_libraries(
      tfhe++
      INTERFACE "-Wl,--start-group"
                $ENV{MKLROOT}/interfaces/fftw3xc/libfftw3xc_gnu.a
                $ENV{MKLROOT}/lib/intel64/libmkl_cdft_core.a
                $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a
                $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a
                $ENV{MKLROOT}/lib/intel64/libmkl_core.a
                $ENV{MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a
                "-Wl,--end-group"
                ${CMAKE_DL_LIBS}
                pthread)
  else()
    target_link_libraries(tfhe++ INTERFACE fftw3)
  endif()
elseif(USE_SPQLIOX_AARCH64)
  target_link_libraries(tfhe++ INTERFACE spqliox_aarch64 xbyak_aarch64)
elseif(USE_HEXL)
  target_link_libraries(tfhe++ INTERFACE spqlios HEXL::hexl)
else()
  target_link_libraries(tfhe++ INTERFACE spqlios)
endif()
