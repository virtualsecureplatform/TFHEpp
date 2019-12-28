file(GLOB tfhe++_sources "*.cpp")

find_library(IPP_CORE libippcore.a PATHS /opt/intel/ipp/lib/intel64)
find_library(IPP_VM libippvm.a PATHS /opt/intel/ipp/lib/intel64)
find_library(IPP_S libipps.a PATHS /opt/intel/ipp/lib/intel64)

add_library(tfhe++
  SHARED
  ${tfhe++_sources}
  )
target_include_directories(tfhe++
  PUBLIC ${PROJECT_SOURCE_DIR}/include 
  /opt/intel/ipp/include
  )
target_link_libraries(tfhe++ ${IPP_CORE} ${IPP_VM} ${IPP_S})