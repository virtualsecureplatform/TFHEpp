file(GLOB tfhe++_sources "*.cpp")


add_library(tfhe++
  SHARED
  ${tfhe++_sources}
  )
target_include_directories(tfhe++
  PUBLIC ${PROJECT_SOURCE_DIR}/include 
  ${PROJECT_SOURCE_DIR}/thirdparties/spqlios
  )
target_link_libraries(tfhe++ INTERFACE spqlios)