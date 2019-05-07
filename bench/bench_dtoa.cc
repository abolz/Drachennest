#include "benchmark/benchmark.h"

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <random>
#include <vector>
#include <chrono>

#include <math.h>

#include "grisu2.h"
#include "grisu3.h"
#include "ryu.h"
#include "charconv.h"
#include "double-conversion.h"

#if _MSC_VER >= 1920
#define HAS_CHARCONV 1
#endif

#define BENCH_GRISU2 1
//#define BENCH_GRISU3 1
//#define BENCH_RYU 1
//#define BENCH_CHARCONV 1
//#define BENCH_DOUBLE_CONVERSION 1
//#define BENCH_SPRINTF 1

#define BENCH_SINGLE 0
#define BENCH_DOUBLE 1
#define BENCH_TO_DECIMAL 1

constexpr int BufSize = 64;

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

struct D2S_Grisu2
{
    static char const* Name() { return "Grisu2"; }
    char* operator()(char* buf, int buflen, float f) const { return grisu2_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) const { return grisu2_Dtoa(buf, buflen, f); }
};

struct D2S_Grisu3
{
    static char const* Name() { return "Grisu3 (Dragon4)"; }
    char* operator()(char* buf, int buflen, float f) const { return grisu3_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) const { return grisu3_Dtoa(buf, buflen, f); }
};

struct D2S_Ryu
{
    static char const* Name() { return "Ryu"; }
    char* operator()(char* buf, int buflen, float f) const { return ryu_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) const { return ryu_Dtoa(buf, buflen, f); }
};

#if HAS_CHARCONV
struct D2S_Charconv
{
    static char const* Name() { return "charconv"; }
    char* operator()(char* buf, int buflen, float f) const { return charconv_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) const { return charconv_Dtoa(buf, buflen, f); }
};
#endif

struct D2S_DoubleConversion
{
    static char const* Name() { return "double-conversion"; }
    char* operator()(char* buf, int buflen, float f) const { return double_conversion_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) const { return double_conversion_Dtoa(buf, buflen, f); }
};

struct D2S_SPrintf
{
    static char const* Name() { return "sprintf"; }
    char* operator()(char* buf, int buflen, float f) const { return buf + std::snprintf(buf, static_cast<size_t>(buflen), "%.9g", f); }
    char* operator()(char* buf, int buflen, double f) const { return buf + std::snprintf(buf, static_cast<size_t>(buflen), "%.17g", f); }
};

//==================================================================================================
//
//==================================================================================================

static constexpr int NumFloats = 1 << 12;

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

#if BENCH_TO_DECIMAL
#if BENCH_GRISU2
inline uint32_t ToDecimal(int& exponent, float  value) { return grisu2_Ftoa(exponent, value); }
inline uint64_t ToDecimal(int& exponent, double value) { return grisu2_Dtoa(exponent, value); }
#endif
#if BENCH_GRISU3
inline uint32_t ToDecimal(int& exponent, float  value) { return grisu3_Ftoa(exponent, value); }
inline uint64_t ToDecimal(int& exponent, double value) { return grisu3_Dtoa(exponent, value); }
#endif
#if BENCH_RYU
inline uint32_t ToDecimal(int& exponent, float  value) { return ryu_Ftoa(exponent, value); }
inline uint64_t ToDecimal(int& exponent, double value) { return ryu_Dtoa(exponent, value); }
#endif

template <typename, typename Float>
static inline void BenchIt(benchmark::State& state, std::vector<Float> const& numbers)
{
    int index = 0;

    uint64_t sum = 0;
    for (auto _ : state)
    {
        int exponent;
        sum += ToDecimal(exponent, numbers[index]) & 0xFF;
        index = (index + 1) & (NumFloats - 1);
    }

    if (sum == UINT64_MAX)
        abort();
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

#if BENCH_GRISU2
using D2S = D2S_Grisu2;
#endif
#if BENCH_GRISU3
using D2S = D2S_Grisu3;
#endif
#if BENCH_RYU
using D2S = D2S_Ryu;
#endif
#if BENCH_CHARCONV
using D2S = D2S_Charconv;
#endif
#if BENCH_DOUBLE_CONVERSION
using D2S = D2S_DoubleConversion;
#endif
#if BENCH_SPRINTF
using D2S = D2S_SPrintf;
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

    RegisterBenchmarks(StrPrintf("Uniform %.1g/%.1g", low, high), numbers);
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
    //const double scale = static_cast<double>(kPow10[digits - 1]);
    //const double scale = 1e8;

    std::vector<double> numbers(NumFloats);

    std::uniform_int_distribution<int64_t> gen(kPow10[digits - 1], kPow10[digits] - 1);
    std::uniform_int_distribution<int> gen_scale(0, digits);

    std::generate(numbers.begin(), numbers.end(), [&]
    {
        const double scale = static_cast<double>( kPow10[gen_scale(random)] );
        return static_cast<double>(gen(random) | 1) / scale;
    });

//  RegisterBenchmarks(StrPrintf("1.%d-digits", digits - 1), numbers);
    RegisterBenchmarks(StrPrintf("%d-digits", digits), numbers);
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
    //const float scale = static_cast<float>(kPow10[digits - 1]);
    //const float scale = 1e4f;

    std::vector<float> numbers(NumFloats);

    std::uniform_int_distribution<int32_t> gen(kPow10[digits - 1], kPow10[digits] - 1);
    std::uniform_int_distribution<int> gen_scale(0, digits);

    std::generate(numbers.begin(), numbers.end(), [&] 
    {
        const float scale = static_cast<float>( kPow10[gen_scale(random)] );
        return static_cast<float>(gen(random) | 1) / scale;
    });

    RegisterBenchmarks(StrPrintf("%d-digits", digits), numbers);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------
static inline void Register_Ints_double(int digits)
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

    std::vector<double> numbers(NumFloats);

    std::uniform_int_distribution<int64_t> gen(kPow10[digits - 1], kPow10[digits] - 1);
    std::generate(numbers.begin(), numbers.end(), [&] {
        return static_cast<double>(gen(random) | 1);
    });

    RegisterBenchmarks(StrPrintf("%d-digit int", digits), numbers);
}

static inline void Register_Ints_single(int digits)
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

    std::vector<float> numbers(NumFloats);

    std::uniform_int_distribution<int32_t> gen(kPow10[digits - 1], kPow10[digits] - 1);
    std::generate(numbers.begin(), numbers.end(), [&] {
        return static_cast<float>(gen(random) | 1);
    });

    RegisterBenchmarks(StrPrintf("%d-digit int", digits), numbers);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    printf("Benchmarking %s\n", D2S::Name());

#if BENCH_DOUBLE
    Register_RandomBits_double();
    Register_Uniform(0.0, 1.0);
    Register_Uniform(0.0, 1.0e+308);

    for (int e = 10; e < 20; ++e) {
        Register_Uniform(std::pow(10.0, e), std::pow(10.0, e+1));
    }

    //for (int d = 1; d <= 15; ++d) {
    //    Register_Ints_double(d);
    //}

    for (int d = 1; d <= 15; ++d) {
        Register_Digits_double(d);
    }
#endif

#if BENCH_SINGLE
    Register_RandomBits_single();
    Register_Uniform(0.0f, 1.0f);
    Register_Uniform(0.0f, 1.0e+38f);

    //for (int d = 1; d <= 7; ++d) {
    //    Register_Ints_single(d);
    //}

    for (int d = 1; d <= 7; ++d) {
        Register_Digits_single(d);
    }
#endif

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}
