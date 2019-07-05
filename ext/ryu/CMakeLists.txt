file(GLOB_RECURSE ryu_sources "*.h" "*.c")

add_library(ryu STATIC ${ryu_sources})

target_include_directories(
    ryu
    PUBLIC
        "${CMAKE_SOURCE_DIR}/ext/ryu"
    )

target_link_libraries(
    ryu
    INTERFACE
        ${DN_PROJECT_NAME}
    )
