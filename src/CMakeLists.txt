file(GLOB drachennest_sources "*.h" "*.cc" "ryu/*.h" "ryu/*.cc")

add_library(drachennest STATIC ${drachennest_sources})

target_link_libraries(
    drachennest
    INTERFACE
        ${DN_INTERFACE}
    )
