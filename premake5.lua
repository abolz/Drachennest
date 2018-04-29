local build_dir = "build/" .. _ACTION

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

    configuration { "debug" }
        targetdir (build_dir .. "/bin/debug")

    configuration { "release" }
        targetdir (build_dir .. "/bin/release")

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

--------------------------------------------------------------------------------
group "Libs"

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
        "test/test.cc",
    }
    includedirs {
        "ext/",
    }
    links {
        "double-conversion",
    }
