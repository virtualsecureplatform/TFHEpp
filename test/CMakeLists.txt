file(GLOB test_sources RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "*.cpp")

file(COPY test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)
file(COPY show_info.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)
file(COPY benchmark_docker.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)

foreach(test_source ${test_sources})
    string( REPLACE ".cpp" "" test_name ${test_source} )
    add_executable(${test_name} ${test_source})
    if(USE_PERF)
        add_compile_definitions(USE_PERF)
        target_link_libraries(${test_name} tfhe++  profiler)
    else()
        target_link_libraries(${test_name} tfhe++)
    endif()
endforeach(test_source ${test_sources})

# For GMP used in INTorus test
find_library(gmp NAMES gmp gmpxx)
if(gmp)
    add_compile_definitions(ENABLE_LSHTEST)
    target_link_libraries(INTorus gmp)
endif()
