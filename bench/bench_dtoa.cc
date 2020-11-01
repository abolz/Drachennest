#include "benchmark/benchmark.h"

#include "ryu_32.h"
#include "ryu_64.h"

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <chrono>
#include <random>
#include <vector>

#include <math.h>

#define BENCH_RYU()             0
#define BENCH_STD_PRINTF()      0
#define BENCH_STD_CHARCONV()    0
#define BENCH_SCHUBFACH()       1
#define BENCH_GRISU2()          0
#define BENCH_GRISU2B()         0
#define BENCH_GRISU3()          0
#define BENCH_DRAGONBOX()       0

#define BENCH_SINGLE()          0
#define BENCH_DOUBLE()          1

#define BENCH_TO_DECIMAL()      0

//==================================================================================================
//
//==================================================================================================

#if BENCH_RYU()
struct D2S
{
    static char const* Name() { return "ryu"; }
    char* operator()(char* buf, int /*buflen*/, float f) const { return ryu::Ftoa(buf, f); }
    char* operator()(char* buf, int /*buflen*/, double f) const { return ryu::Dtoa(buf, f); }
#if BENCH_TO_DECIMAL()
    static ryu::FloatingDecimal64 ToDec(double value) { return ryu::ToDecimal64(value); }
#endif
};
#endif

#if BENCH_STD_PRINTF()
struct D2S
{
    static char const* Name() { return "std::printf"; }
    char* operator()(char* buf, int buflen, float f) const { return buf + std::snprintf(buf, static_cast<size_t>(buflen), "%.9g", f); }
    char* operator()(char* buf, int buflen, double f) const { return buf + std::snprintf(buf, static_cast<size_t>(buflen), "%.17g", f); }
};
#endif

#if BENCH_STD_CHARCONV()
#include <charconv>
struct D2S
{
#if 0
    static char const* Name() { return "std::charconv::general"; }
    char* operator()(char* buf, int buflen, float f) const { return std::to_chars(buf, buf + buflen, f, std::chars_format::general).ptr; }
    char* operator()(char* buf, int buflen, double f) const { return std::to_chars(buf, buf + buflen, f, std::chars_format::general).ptr; }
#else
    static char const* Name() { return "std::charconv"; }
    char* operator()(char* buf, int buflen, float f) const { return std::to_chars(buf, buf + buflen, f).ptr; }
    char* operator()(char* buf, int buflen, double f) const { return std::to_chars(buf, buf + buflen, f).ptr; }
#endif
};
#endif

#if BENCH_SCHUBFACH()
#include "schubfach_32.h"
#include "schubfach_64.h"
struct D2S
{
    static char const* Name() { return "schubfach"; }
    char* operator()(char* buf, int /*buflen*/, float f) const { return schubfach::Ftoa(buf, f); }
    char* operator()(char* buf, int /*buflen*/, double f) const { return schubfach::Dtoa(buf, f); }
};
#endif

#if BENCH_GRISU2()
#include "grisu2.h"
struct D2S
{
    static char const* Name() { return "grisu2"; }
    char* operator()(char* buf, int /*buflen*/, double f) const { return grisu2::Dtoa(buf, f); }
};
#endif

#if BENCH_GRISU2B()
#include "grisu2b.h"
struct D2S
{
    static char const* Name() { return "grisu2b"; }
    char* operator()(char* buf, int /*buflen*/, double f) const { return grisu2b::Dtoa(buf, f); }
};
#endif

#if BENCH_GRISU3()
#include "grisu3.h"
struct D2S
{
    static char const* Name() { return "grisu3"; }
    char* operator()(char* buf, int /*buflen*/, double f) const { return grisu3::Dtoa(buf, f); }
};
#endif

#if BENCH_DRAGONBOX()
#include "dragonbox.h"
struct D2S
{
    static char const* Name() { return "dragonbox"; }
    char* operator()(char* buf, int /*buflen*/, double f) const { return dragonbox::Dtoa(buf, f); }
};
#endif

//==================================================================================================
//
//==================================================================================================

template <typename Target, typename Source>
static Target ReinterpretBits(Source const& source)
{
    static_assert(sizeof(Target) == sizeof(Source), "size mismatch");

    Target target;
    std::memcpy(&target, &source, sizeof(Source));
    return target;
}

class JenkinsRandom
{
    // A small noncryptographic PRNG
    // http://burtleburtle.net/bob/rand/smallprng.html

    uint32_t a;
    uint32_t b;
    uint32_t c;
    uint32_t d;

    static uint32_t Rotate(uint32_t value, int n) {
        return (value << n) | (value >> (32 - n));
    }

    uint32_t Gen() {
        const uint32_t e = a - Rotate(b, 27);
        a = b ^ Rotate(c, 17);
        b = c + d;
        c = d + e;
        d = e + a;
        return d;
    }

public:
    using result_type = uint32_t;

    static constexpr uint32_t min() { return 0; }
    static constexpr uint32_t max() { return UINT32_MAX; }

    explicit JenkinsRandom(uint32_t seed = 0) {
        a = 0xF1EA5EED;
        b = seed;
        c = seed;
        d = seed;
        for (int i = 0; i < 20; ++i) {
            static_cast<void>(Gen());
        }
    }

    uint32_t operator()() { return Gen(); }
};

static JenkinsRandom random;

//==================================================================================================
//
//==================================================================================================

static constexpr int BufSize = 64;
static constexpr int NumFloats = 1 << 13;

//template <typename Float>
//static inline void PrintFloat(Float v)
//{
//    char buf[32];
//    char* end = D2S_DoubleConversion{}(buf, 32, v);
//    end[0] = '\0';
//    printf("%s\n", buf);
//}

template <typename ...Args>
static inline char const* StrPrintf(char const* format, Args&&... args)
{
    char buf[1024];
    snprintf(buf, 1024, format, std::forward<Args>(args)...);
#ifdef _MSC_VER
    return _strdup(buf); // leak...
#else
    return strdup(buf); // leak...
#endif
}

#if BENCH_TO_DECIMAL()
template <typename D2S, typename Float>
static inline void BenchIt(benchmark::State& state, std::vector<Float> const& numbers)
{
    int index = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(D2S::ToDec(numbers[index]));
        index = (index + 1) & (NumFloats - 1);
    }
}
#else
template <typename D2S, typename Float>
static inline void BenchIt(benchmark::State& state, std::vector<Float> const& numbers)
{
    D2S d2s;
 
    int index = 0;

    uint64_t sum = 0;
    for (auto _ : state)
    {
        char buffer[BufSize];
        d2s(buffer, BufSize, numbers[index]);
        sum += static_cast<unsigned char>(buffer[0]);
        index = (index + 1) & (NumFloats - 1);
    }

    if (sum == UINT64_MAX)
        abort();
}
#endif

template <typename Float>
static inline void RegisterBenchmarks(char const* name, std::vector<Float> const& numbers)
{
    const char* float_name = sizeof(Float) == 4 ? "single" : "double";
    auto* bench = benchmark::RegisterBenchmark(StrPrintf("%s - %s   ", float_name, name), BenchIt<D2S, Float>, numbers);

    bench->ComputeStatistics("min", [](const std::vector<double>& v) -> double {
        return *(std::min_element(std::begin(v), std::end(v)));
    });
    bench->Repetitions(3);
    bench->ReportAggregatesOnly();
}

//----------------------------------------------------------------------------------------------------------
//
//----------------------------------------------------------------------------------------------------------
static inline void Register_RandomBits_double()
{
    std::vector<double> numbers(NumFloats);

    std::uniform_int_distribution<uint64_t> gen(1, 0x7FF0000000000000ull - 1);
    std::generate(numbers.begin(), numbers.end(), [&] { return ReinterpretBits<double>(gen(random)); });

    RegisterBenchmarks("Random-bits", numbers);
}

static inline void Register_RandomBits_single()
{
    std::vector<float> numbers(NumFloats);

    std::uniform_int_distribution<uint32_t> gen(1, 0x7F800000u - 1);
    std::generate(numbers.begin(), numbers.end(), [&] { return ReinterpretBits<float>(gen(random)); });

    RegisterBenchmarks("Random-bits", numbers);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------
template <typename Float>
static inline void Register_Uniform(Float low, Float high)
{
    std::vector<Float> numbers(NumFloats);

    std::uniform_real_distribution<Float> gen(low, high);
    std::generate(numbers.begin(), numbers.end(), [&] { return gen(random); });

    RegisterBenchmarks(StrPrintf("Uniform %.1g/%.1g", low, high), numbers);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------
static constexpr int64_t kPow10_i64[] = {
    1,
    10,
    100,
    1000,
    10000,
    100000,
    1000000,
    10000000,
    100000000,
    1000000000,
    10000000000,
    100000000000,
    1000000000000,
    10000000000000,
    100000000000000,
    1000000000000000,
    10000000000000000,
    100000000000000000,
    1000000000000000000,
};
static constexpr double kPow10_f64[] = {
    1.0e+00,
    1.0e+01,
    1.0e+02,
    1.0e+03,
    1.0e+04,
    1.0e+05,
    1.0e+06,
    1.0e+07,
    1.0e+08,
    1.0e+09,
    1.0e+10,
    1.0e+11,
    1.0e+12,
    1.0e+13,
    1.0e+14,
    1.0e+15,
    1.0e+16,
    1.0e+17,
    1.0e+18,
    1.0e+19,
    1.0e+20,
    1.0e+21,
    1.0e+22,
};

static inline void Register_Digits_double(const char* name, int digits, int e10)
{
    assert(digits >= 1);
    assert(digits <= 18);
    assert(e10 >= -22);
    assert(e10 <= 22);

    std::vector<double> numbers(NumFloats);

    std::uniform_int_distribution<int64_t> gen(kPow10_i64[digits - 1], kPow10_i64[digits] - 1);

    std::generate(numbers.begin(), numbers.end(), [&] {
        int64_t n = gen(random);
        if (n % 10 == 0)
            n |= 1;
        double v = static_cast<double>(n);
        if (e10 < 0)
            v /= kPow10_f64[-e10];
        else
            v *= kPow10_f64[e10];
        //PrintFloat(v);
        return v;
    });

    RegisterBenchmarks(name, numbers);
}

static constexpr int32_t kPow10_i32[] = {
    1,
    10,
    100,
    1000,
    10000,
    100000,
    1000000,
    10000000,
    100000000,
    1000000000,
};
static constexpr float kPow10_f32[] = {
    1.0e+00f,
    1.0e+01f,
    1.0e+02f,
    1.0e+03f,
    1.0e+04f,
    1.0e+05f,
    1.0e+06f,
    1.0e+07f,
    1.0e+08f,
    1.0e+09f,
    1.0e+10f,
};

static inline void Register_Digits_single(const char* name, int digits, int e10)
{
    assert(digits >= 1);
    assert(digits <= 9);
    assert(e10 >= -10);
    assert(e10 <= 10);

    std::vector<float> numbers(NumFloats);

    std::uniform_int_distribution<int32_t> gen(kPow10_i32[digits - 1], kPow10_i32[digits] - 1);

    std::generate(numbers.begin(), numbers.end(), [&] {
        int32_t n = gen(random);
        if (n % 10 == 0)
            n |= 1;
        float v = static_cast<float>(n);
        if (e10 < 0)
            v /= kPow10_f32[-e10];
        else
            v *= kPow10_f32[e10];
        //PrintFloat(v);
        return v;
    });

    RegisterBenchmarks(name, numbers);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

static inline std::vector<double> GenRandomDigitData_double(int digits, int count)
{
    std::uniform_real_distribution<double> gen(1, 2);

    std::vector<double> result;
    result.resize(count);

    for (int i = 0; i < count; ++i)
    {
        const double d = gen(random);
        const double rounded = ryu::Round10(d, -digits);
        result[i] = rounded;
    }

    return result;
}

static inline void Register_RandomDigits_double(const char* name, int digits)
{
    RegisterBenchmarks(name, GenRandomDigitData_double(digits, NumFloats));
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

static inline std::vector<float> GenRandomDigitData_float(int digits, int count)
{
    std::uniform_real_distribution<float> gen(1, 2);

    std::vector<float> result;
    result.resize(count);

    for (int i = 0; i < count; ++i)
    {
        const float d = gen(random);
        const float rounded = ryu::Round10(d, -digits);
        result[i] = rounded;
    }

    return result;
}

static inline void Register_RandomDigits_float(const char* name, int digits)
{
    RegisterBenchmarks(name, GenRandomDigitData_float(digits, NumFloats));
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

int main(int argc, char** argv)
{
#if defined(__clang__)
    printf("clang %d.%d\n", __clang_major__, __clang_minor__);
#elif defined(__GNUC__)
    printf("gcc %s\n", __VERSION__);
#elif defined(_MSC_VER)
    printf("msc %d\n", _MSC_FULL_VER);
#endif

    printf("Preparing benchmarks...\n");

#if BENCH_DOUBLE()
    Register_RandomBits_double();
    Register_RandomBits_double();
    Register_RandomBits_double();
    Register_Uniform(0.0, 1.0);
    Register_Uniform(0.0, 1.0e+308);
    Register_Uniform(1.0, 2.0);

#if 0
    for (int d = 1; d <= 18; ++d) {
        for (int e = -10; e <= 10; e += 1) {
            Register_Digits_double(StrPrintf("%2d,%3d", d, e), d, e);
        }
    }
#else
    for (int d = 0; d <= 16; ++d)
    {
        Register_RandomDigits_double(StrPrintf("1.%d-digits", d), d);
    }
#endif
#endif

#if BENCH_SINGLE()
    Register_RandomBits_single();
    Register_RandomBits_single();
    Register_RandomBits_single();
    Register_Uniform(0.0f, 1.0f);
    Register_Uniform(0.0f, 1.0e+38f);
    Register_Uniform(1.0f, 2.0f);

#if 0
    for (int d = 1; d <= 10; ++d) {
        for (int e = -10; e <= 10; e += 1) {
            Register_Digits_single(StrPrintf("%2d,%3d", d, e), d, e);
        }
    }
#else
    for (int d = 0; d <= 8; ++d)
    {
        Register_RandomDigits_float(StrPrintf("1.%d-digits", d), d);
    }
#endif
#endif

    printf("Benchmarking %s\n", D2S::Name());

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}
