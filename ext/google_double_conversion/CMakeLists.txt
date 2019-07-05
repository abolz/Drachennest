file(GLOB_RECURSE google_double_conversion_sources "*.h" "*.cc")

add_library(google_double_conversion STATIC ${google_double_conversion_sources})

target_include_directories(
    google_double_conversion
    PUBLIC
        "${CMAKE_SOURCE_DIR}/ext/google_double_conversion/"
    )

target_link_libraries(
    google_double_conversion
    INTERFACE
        ${DN_INTERFACE}
    )
