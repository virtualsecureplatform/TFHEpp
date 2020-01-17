SET(test_sources
    circuitbootstrapping.cpp
    cmux.cpp
    externalproduct.cpp
    gate.cpp
    gatebootstrapping.cpp
    nand.cpp
    polymul.cpp
    rom.cpp
    tlweenc.cpp
    trlweenc.cpp
    write.cpp)


foreach(test_source ${test_sources})
    string( REPLACE ".cpp" "" test_name ${test_source} )
    add_executable(${test_name} ${test_source})
    target_link_libraries(${test_name} tfhe++)
endforeach(test_source ${test_sources})