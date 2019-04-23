#include "benchmark/benchmark.h"

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstring>

#include <algorithm>
#include <random>
#include <vector>
#include <chrono>

#include <math.h>

#include "charconv.h"
#include "double-conversion.h"
#include "grisu2.h"
#include "grisu3.h"
#include "ryu.h"

#if _MSC_VER >= 1920
#define HAS_CHARCONV 1
#endif

#define BENCH_SINGLE 1
#define BENCH_DOUBLE 1

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

struct Grisu2
{
    static char const* Name() { return "grisu2"; }
    char* operator()(char* buf, int buflen, float f) const { return grisu2_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) const { return grisu2_Dtoa(buf, buflen, f); }
};

struct Grisu3
{
    static char const* Name() { return "grisu3"; }
    char* operator()(char* buf, int buflen, float f) const { return grisu3_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) const { return grisu3_Dtoa(buf, buflen, f); }
};

struct Ryu
{
    static char const* Name() { return "ryu"; }
    char* operator()(char* buf, int buflen, float f) const { return ryu_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) const { return ryu_Dtoa(buf, buflen, f); }
};

#if HAS_CHARCONV
struct Charconv
{
    static char const* Name() { return "charconv"; }
    char* operator()(char* buf, int buflen, float f) const { return charconv_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) const { return charconv_Dtoa(buf, buflen, f); }
};
#endif

struct DoubleConversion
{
    static char const* Name() { return "double-conversion"; }
    char* operator()(char* buf, int buflen, float f) const { return double_conversion_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) const { return double_conversion_Dtoa(buf, buflen, f); }
};

//==================================================================================================
//
//==================================================================================================

static constexpr int NumFloats = 1 << 12;

template <typename ...Args>
static inline char const* SPrintf(char const* format, Args&&... args)
{
    char buf[1024];
    snprintf(buf, 1024, format, std::forward<Args>(args)...);
#ifdef _MSC_VER
    return _strdup(buf); // leak...
#else
    return strdup(buf); // leak...
#endif
}

template <typename D2S, typename Float>
static inline void BenchIt(benchmark::State& state, std::vector<Float> const& numbers)
{
    D2S d2s;

    int index = 0;

    for (auto _ : state)
    {
        char buffer[32];
        benchmark::DoNotOptimize(d2s(buffer, 32, numbers[index]));
        index = (index + 1) & (NumFloats - 1);
    }
}

template <typename Float>
static inline void RegisterBenchmarks(char const* name, std::vector<Float> const& numbers)
{
    using D2S = Grisu2;
    //using D2S = Grisu3;
    //using D2S = Ryu;
    //using D2S = Charconv;
    //using D2S = DoubleConversion;

    //const char* float_name = sizeof(Float) == 4 ? "F32" : "F64";
    const char* float_name = sizeof(Float) == 4 ? "single" : "double";
    //const char* float_name = SPrintf("p=%d", std::numeric_limits<Float>::digits);

//  auto* bench = benchmark::RegisterBenchmark(SPrintf("%s - %s (%s)   ", float_name, name, D2S::Name()), BenchIt<D2S, Float>, numbers);
    auto* bench = benchmark::RegisterBenchmark(SPrintf("%s - %s   ", float_name, name), BenchIt<D2S, Float>, numbers);

    bench->ComputeStatistics("min", [](const std::vector<double>& v) -> double {
        return *(std::min_element(std::begin(v), std::end(v)));
    });
    bench->Repetitions(3);
    bench->ReportAggregatesOnly(true);
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

    RegisterBenchmarks(SPrintf("Uniform %.1g/%.1g", low, high), numbers);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------
static inline void Register_Digits_double(int digits)
{
    assert(digits >= 1);
    assert(digits <= 15);
    static constexpr int64_t kPow10[] = {
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
    };
    const double scale = static_cast<double>(kPow10[digits - 1]);

    std::vector<double> numbers(NumFloats);

    std::uniform_int_distribution<int64_t> gen(kPow10[digits - 1], kPow10[digits] - 1);
    std::generate(numbers.begin(), numbers.end(), [&] {
        return static_cast<double>(gen(random) | 1) / scale;
    });

    RegisterBenchmarks(SPrintf("1.%d-digits", digits - 1), numbers);
}

static inline void Register_Digits_single(int digits)
{
    assert(digits >= 1);
    assert(digits <= 7);
    static constexpr int32_t kPow10[] = {
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
    const float scale = static_cast<float>(kPow10[digits - 1]);

    std::vector<float> numbers(NumFloats);

    std::uniform_int_distribution<int32_t> gen(kPow10[digits - 1], kPow10[digits] - 1);
    std::generate(numbers.begin(), numbers.end(), [&] {
        return static_cast<float>(gen(random) | 1) / scale;
    });

    RegisterBenchmarks(SPrintf("1.%d-digits", digits - 1), numbers);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------
#if 0

static void BM_RandomDigit32(benchmark::State& state, int digits, float scale)
{
    assert(digits >= 1);
    assert(digits <= 7);

    constexpr int const NumFloats = 1 << 12;

    if (scale <= 0)
        scale = static_cast<float>(kPow10[digits - 1]);

    std::uniform_int_distribution<int32_t> gen(kPow10[digits - 1], kPow10[digits] - 1);
    std::vector<float> numbers(NumFloats);
    std::generate(numbers.begin(), numbers.end(), [&] {
        return static_cast<float>(gen(random) | 1) / scale;
    });

    int index = 0;
    for (auto _: state) {
        char buffer[32];
        benchmark::DoNotOptimize(Ftoa(buffer, 32, numbers[index]));
        index = (index + 1) & (NumFloats - 1);
    }
}

BENCHMARK_CAPTURE(BM_RandomDigit32, _1_1, 1, -1.0f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _1_2, 2, -1.0f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _1_3, 3, -1.0f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _1_4, 4, -1.0f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _1_5, 5, -1.0f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _1_6, 6, -1.0f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _1_7, 7, -1.0f);

BENCHMARK_CAPTURE(BM_RandomDigit32, _1_over_10e4, 1, 1e4f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _2_over_10e4, 2, 1e4f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _3_over_10e4, 3, 1e4f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _4_over_10e4, 4, 1e4f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _5_over_10e4, 5, 1e4f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _6_over_10e4, 6, 1e4f);
BENCHMARK_CAPTURE(BM_RandomDigit32, _7_over_10e4, 7, 1e4f);

#endif

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

int main(int argc, char** argv)
{
#if BENCH_DOUBLE
    Register_RandomBits_double();

    Register_Uniform(0.0, 1.0);
    Register_Uniform(0.0, 1.0e+308);
    for (int e = 10; e < 20; ++e) {
        Register_Uniform(std::pow(10.0, e), std::pow(10.0, e+1));
    }

    for (int d = 1; d <= 15; ++d) {
        Register_Digits_double(d);
    }
#endif

#if BENCH_SINGLE
    Register_RandomBits_single();

    Register_Uniform(0.0f, 1.0f);
    Register_Uniform(0.0f, 1.0e+38f);

    for (int d = 1; d <= 7; ++d) {
        Register_Digits_single(d);
    }
#endif

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}
