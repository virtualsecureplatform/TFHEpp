if(USE_80BIT_SECURITY)
    add_compile_definitions(USE_80BIT_SECURITY)
elseif(USE_CGGI19)
    add_compile_definitions(USE_CGGI19)
endif()

file(GLOB test_sources RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "*.cpp")

foreach(test_source ${test_sources})
    string( REPLACE ".cpp" "" test_name ${test_source} )
    add_executable(${test_name} ${test_source})
    target_link_libraries(${test_name} tfhe++  profiler)
endforeach(test_source ${test_sources})

# For GMP used in INTorus test
find_library(gmp NAMES gmp gmpxx)
if(gmp)
    add_compile_definitions(ENABLE_LSHTEST)
    target_link_libraries(INTorus gmp)
endif()