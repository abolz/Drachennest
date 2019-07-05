file(GLOB_RECURSE google_benchmark_sources "*.h" "*.cc")

add_library(google_benchmark STATIC ${google_benchmark_sources})

target_include_directories(
    google_benchmark
    PUBLIC
        "${CMAKE_SOURCE_DIR}/ext/google_benchmark/include"
    )

target_link_libraries(
    google_benchmark 
    INTERFACE 
        ${DN_INTERFACE}
    PRIVATE
        shlwapi
    )
