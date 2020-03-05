file(GLOB tfhe++_sources "*.cpp")

if(USE_80BIT_SECURITY)
    add_definitions(-DENABLE_80BIT_SECURITY)
endif()

add_library(tfhe++
  ${tfhe++_sources}
  )
target_include_directories(tfhe++
  PUBLIC ${PROJECT_SOURCE_DIR}/include 
  ${PROJECT_SOURCE_DIR}/thirdparties/spqlios
  ${PROJECT_SOURCE_DIR}/thirdparties/randen
  )
target_link_libraries(tfhe++ INTERFACE spqlios randen)
