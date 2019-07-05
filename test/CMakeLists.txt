set(test_dtoa_sources
    "catch.hpp"
    "catch_main.cc"
    "scan_number.h"
    "test_dtoa.cc"
    )

add_executable(test_dtoa ${test_dtoa_sources})

target_include_directories(
    test_dtoa
    PUBLIC
        "${CMAKE_SOURCE_DIR}/ext/google_benchmark/"
        "${CMAKE_SOURCE_DIR}/ext/google_double_conversion/"
        "${CMAKE_SOURCE_DIR}/src/"
    )

target_link_libraries(
    test_dtoa
    INTERFACE
        ${DN_INTERFACE}
    PRIVATE
        drachennest
        google_double_conversion
    )
