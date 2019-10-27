set(bench_dtoa_sources "bench_dtoa.cc")

add_executable(bench_dtoa ${bench_dtoa_sources})

target_include_directories(
    bench_dtoa
    PUBLIC
        "${CMAKE_SOURCE_DIR}/ext/"
        "${CMAKE_SOURCE_DIR}/src/"
    )

target_link_libraries(
    bench_dtoa
    INTERFACE
        ${DN_INTERFACE}
    PRIVATE
        google_benchmark
        google_double_conversion
        ryu
    )
