set(test_all_sources
    "catch.hpp"
    "catch_main.cc"
    "scan_number.h"
    "test_dtoa.cc"
    "test_strtod.cc"
    )

add_executable(test_all ${test_all_sources})

target_include_directories(
    test_all
    PUBLIC
        "${CMAKE_SOURCE_DIR}/ext/google_benchmark/"
        "${CMAKE_SOURCE_DIR}/ext/google_double_conversion/"
        "${CMAKE_SOURCE_DIR}/src/"
    )

target_link_libraries(
    test_all
    INTERFACE
        ${DN_INTERFACE}
    PRIVATE
        drachennest
        google_double_conversion
    )
