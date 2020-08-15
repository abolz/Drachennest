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
        drachennest
        google_benchmark
        google_double_conversion
        ryu
    )

set(bench_strtod_sources "bench_strtod.cc")

add_executable(bench_strtod ${bench_strtod_sources})

target_include_directories(
    bench_strtod
    PUBLIC
        "${CMAKE_SOURCE_DIR}/ext/"
        "${CMAKE_SOURCE_DIR}/src/"
    )

target_link_libraries(
    bench_strtod
    INTERFACE
        ${DN_INTERFACE}
    PRIVATE
        drachennest
        google_benchmark
        google_double_conversion
        ryu
    )
