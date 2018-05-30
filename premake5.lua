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
    platforms { "x64", "x86" }

    filter { "platforms:x64" }
        architecture "x86_64"

    filter { "platforms:x86" }
        architecture "x86"

    filter {}

    location    (build_dir)
    objdir      (build_dir .. "/obj")

    warnings "Extra"

    -- exceptionhandling "Off"
    -- rtti "Off"

    flags {
        "StaticRuntime",
    }

    configuration { "debug", "x64" }
        targetdir (build_dir .. "/bin/debug_x64")
    configuration { "release", "x64" }
        targetdir (build_dir .. "/bin/release_x64")
    configuration { "debug", "x86" }
        targetdir (build_dir .. "/bin/debug_x86")
    configuration { "release", "x86" }
        targetdir (build_dir .. "/bin/release_x86")

    configuration { "debug" }
        defines { "_DEBUG" }
        symbols "On"

    configuration { "release" }
        defines { "NDEBUG" }
        symbols "On" -- for profiling...
        optimize "Full"
            -- On ==> -O2
            -- Full ==> -O3

    configuration { "gmake" }
        buildoptions {
            "-std=c++11",
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
            buildoptions {
                "-std=c++14",
            }
    end

    if _OPTIONS["linkflags"] then
        configuration {}
            linkoptions {
                _OPTIONS["linkflags"],
            }
    end

--------------------------------------------------------------------------------
group "Libs"

--project "benchmark"
--    language "C++"
--    kind "StaticLib"
--    files {
--        "ext/benchmark/*.cc",
--        "ext/benchmark/*.h",
--    }
--    includedirs {
--        "ext/",
--    }

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
        "src/**.h",
        "src/**.cc",
        "test/catch.hpp",
        "test/catch_main.cc",
        "test/test_dtoa.cc",
        "test/test_strtod.cc",
    }
    includedirs {
        "ext/",
    }
    links {
        "double-conversion",
    }

--project "bench_dtoa"
--    language "C++"
--    kind "ConsoleApp"
--    files {
--        "src/**.h",
--        "src/**.cc",
--        "test/bench_dtoa.cc",
--    }
--    includedirs {
--        "ext/",
--    }
--    links {
--        "benchmark",
--    }
--    configuration { "windows" }
--        links { "shlwapi" }
--
--project "bench_strtod"
--    language "C++"
--    kind "ConsoleApp"
--    files {
--        "src/**.h",
--        "src/**.cc",
--        "test/bench_strtod.cc",
--    }
--    includedirs {
--        "ext/",
--    }
--    links {
--        "benchmark",
--    }
--    configuration { "windows" }
--        links { "shlwapi" }
