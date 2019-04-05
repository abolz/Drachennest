local build_dir = "build/" .. _ACTION

newoption {
    trigger = "cxxflags",
    description = "Additional build options",
}

newoption {
    trigger = "linkflags",
    description = "Additional linker options",
}

--------------------------------------------------------------------------------
workspace "Grisu"
    configurations { "release", "debug" }
    architecture "x64"

    location    (build_dir)
    objdir      (build_dir .. "/obj")

    warnings "Extra"

    -- exceptionhandling "Off"
    -- rtti "Off"

    flags {
        "StaticRuntime",
    }

    configuration { "debug" }
        defines { "_DEBUG" }
        symbols "On"

    configuration { "release" }
        defines { "NDEBUG" }
        symbols "On" -- for profiling...
        optimize "On"
            -- On ==> -O2
            -- Full ==> -O3

    configuration { "gmake" }
        buildoptions {
            "-std=c++14",
            "-march=native",
            "-Wformat",
            -- "-Wsign-compare",
            -- "-Wsign-conversion",
            -- "-pedantic",
            -- "-fno-omit-frame-pointer",
            -- "-ftime-report",
        }

    configuration { "gmake", "debug", "linux" }
        buildoptions {
            -- "-fno-omit-frame-pointer",
            -- "-fsanitize=undefined",
            -- "-fsanitize=address",
            -- "-fsanitize=memory",
            -- "-fsanitize-memory-track-origins",
        }
        linkoptions {
            -- "-fsanitize=undefined",
            -- "-fsanitize=address",
            -- "-fsanitize=memory",
        }

    configuration { "vs*" }
        buildoptions {
            "/utf-8",
            -- "/std:c++latest",
            -- "/EHsc",
            -- "/arch:AVX2",
            -- "/GR-",
        }
        defines {
            -- "_CRT_SECURE_NO_WARNINGS=1",
            -- "_SCL_SECURE_NO_WARNINGS=1",
            -- "_HAS_EXCEPTIONS=0",
        }

    configuration { "windows" }
        characterset "MBCS"

    if _OPTIONS["cxxflags"] then
        configuration {}
            buildoptions {
                _OPTIONS["cxxflags"],
            }
    else
        configuration { "gmake" }
            -- buildoptions {
            --     "-std=c++17",
            -- }
        configuration { "vs*" }
            -- buildoptions {
            --     "/std:c++latest",
            --     -- "/arch:AVX2",
            -- }
    end

    if _OPTIONS["linkflags"] then
        configuration {}
            linkoptions {
                _OPTIONS["linkflags"],
            }
    end

--------------------------------------------------------------------------------
group "Libs"

project "benchmark"
    language "C++"
    kind "StaticLib"
    files {
        "ext/benchmark/*.cc",
        "ext/benchmark/*.h",
    }
    includedirs {
        "ext/",
    }

project "double-conversion"
    language "C++"
    kind "StaticLib"
    files {
        "ext/double-conversion/*.cc",
        "ext/double-conversion/*.h",
    }
    includedirs {
        "ext/",
    }

--------------------------------------------------------------------------------
group "Tests"

project "test"
    language "C++"
    kind "ConsoleApp"
    files {
        -- "lib/**.cc",
        -- "lib/**.h",
        "lib/grisu2.cc",
        "lib/grisu2.h",
        "src/**.cc",
        "src/**.h",
        "test/catch_main.cc",
        "test/catch.hpp",
        "test/test_dtoa.cc",
    }
    includedirs {
        "ext/",
        "lib/",
    }
    links {
        "double-conversion",
    }
    configuration { "gmake" }
        buildoptions {
            "-Wsign-compare",
            "-Wsign-conversion",
        }

project "bench_dtoa"
    language "C++"
    kind "ConsoleApp"
    files {
        "bench/bench_dtoa.cc",
        "lib/**.cc",
        "lib/**.h",
        "src/**.cc",
        "src/**.h",
    }
    includedirs {
        "ext/",
        "lib/",
    }
    links {
        "benchmark",
        "double-conversion",
    }
    configuration { "windows" }
        links { "shlwapi" }
