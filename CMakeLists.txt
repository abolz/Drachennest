cmake_minimum_required(VERSION 3.9)

#-------------------------------------------------------------------------------
# Drachennest
# Algorithms to convert binary to decimal floating-point
#-------------------------------------------------------------------------------

project(Drachennest
    VERSION 0.0.0.1
    DESCRIPTION "Algorithms to convert binary to decimal floating-point"
    LANGUAGES C CXX
    )

set(DN_INTERFACE ${PROJECT_NAME})

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

add_library(${DN_INTERFACE} INTERFACE)

target_compile_features(
    ${DN_INTERFACE}
    INTERFACE
        cxx_std_17
    )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

add_subdirectory(bench)
add_subdirectory(ext)
add_subdirectory(src)
add_subdirectory(test)
