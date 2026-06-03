file(GLOB_RECURSE test_sources RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "*.cpp")

file(COPY test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)
file(COPY show_info.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)
file(COPY benchmark_docker.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)

set(test_directories common tfhe bfv clpx ckks)

foreach(test_source ${test_sources})
    get_filename_component(test_name ${test_source} NAME_WE)
    get_filename_component(test_dir ${test_source} DIRECTORY)
    add_executable(${test_name} ${test_source})
    if(test_dir STREQUAL "")
        set(test_output_dir "${CMAKE_CURRENT_BINARY_DIR}")
    else()
        set(test_output_dir "${CMAKE_CURRENT_BINARY_DIR}/${test_dir}")
    endif()
    set_target_properties(${test_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${test_output_dir}"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG "${test_output_dir}"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE "${test_output_dir}"
        RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${test_output_dir}"
        RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL "${test_output_dir}")
    if(USE_PERF)
        add_compile_definitions(USE_PERF)
        target_link_libraries(${test_name} tfhe++ AES profiler)
    else()
        target_link_libraries(${test_name} tfhe++ AES)
    endif()
    if(test_name STREQUAL "ascon" OR
       test_name STREQUAL "ascon_aes_decrypt_benchmark")
        target_sources(${test_name} PRIVATE
            ${PROJECT_SOURCE_DIR}/thirdparties/ascon-c/crypto_aead/asconaead128/ref/aead.c
            ${PROJECT_SOURCE_DIR}/thirdparties/ascon-c/crypto_hash/asconxof128/ref/hash.c)
        target_include_directories(${test_name} PRIVATE
            ${PROJECT_SOURCE_DIR}/thirdparties/ascon-c/tests)
    endif()
    if(test_dir IN_LIST test_directories)
        list(APPEND test_${test_dir}_targets ${test_name})
    endif()
endforeach(test_source ${test_sources})

foreach(test_dir ${test_directories})
    set(test_targets_var "test_${test_dir}_targets")
    set(test_targets ${${test_targets_var}})
    if(test_targets)
        add_custom_target(test_${test_dir} DEPENDS ${test_targets})
    endif()
endforeach()

# For GMP used in INTorus test
find_library(gmp NAMES gmp gmpxx)
if(gmp)
    add_compile_definitions(ENABLE_LSHTEST)
    target_link_libraries(INTorus gmp)
endif()
