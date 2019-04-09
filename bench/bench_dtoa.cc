#include "benchmark/benchmark.h"
#include "double-conversion/double-conversion.h"

#include <cfloat>
#include <climits>
#include <cmath>
#include <math.h>
#include <cstring>
#include <algorithm>
#include <charconv>
#include <random>
#include <vector>

#include "floaxie.h"
//#include "fmt.h"
#include "grisu2.h"
#include "milo.h"

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

static inline char* Dtoa(char* buffer, int buffer_length, double value)
{
    return grisu2_Dtoa(buffer, buffer + buffer_length, value);
    //return milo_Dtoa(buffer, buffer + buffer_length, value);
    //return floaxie_Dtoa(buffer, buffer + buffer_length, value);
    //return fmt_Dtoa(buffer, buffer + buffer_length, value);
    //return std::to_chars(buffer, buffer + buffer_length, value, std::chars_format::general).ptr;
}

static inline char* Ftoa(char* buffer, int buffer_length, float value)
{
    return grisu2_Ftoa(buffer, buffer + buffer_length, value);
    //return fmt_Ftoa(buffer, buffer + buffer_length, value);
    //return floaxie_Ftoa(buffer, buffer + buffer_length, value);
    //return std::to_chars(buffer, buffer + buffer_length, value, std::chars_format::general).ptr;
}

static inline char* Float32ToChars(char* buf, int buflen, float f)
{
    auto const& f2s = double_conversion::DoubleToStringConverter::EcmaScriptConverter();
    double_conversion::StringBuilder builder(buf, buflen);
    f2s.ToShortestSingle(f, &builder);
    int const length = builder.position();
    builder.Finalize();
    return buf + length;
}

static inline char* Float64ToChars(char* buf, int buflen, double f)
{
    auto const& f2s = double_conversion::DoubleToStringConverter::EcmaScriptConverter();
    double_conversion::StringBuilder builder(buf, buflen);
    f2s.ToShortest(f, &builder);
    int const length = builder.position();
    builder.Finalize();
    return buf + length;
}

//----------------------------------------------------------------------------------------------------------
//
//----------------------------------------------------------------------------------------------------------
#if 1

template <typename Float, typename Uint>
static void BM_RandomBits(benchmark::State& state, Uint lower, Uint upper) {
    constexpr int const NumFloats = 1 << 20;

    std::uniform_int_distribution<Uint> gen(lower, upper);
    std::vector<Float> numbers(NumFloats);
    std::generate(numbers.begin(), numbers.end(), [&] { return ReinterpretBits<Float>(gen(random)); });

    int index = 0;
    for (auto _ : state) {
        char buffer[64];
        benchmark::DoNotOptimize(Dtoa(buffer, 64, numbers[index]));
        index = (index + 1) & (NumFloats - 1);
    }
}

static inline void BM_RandomBits_double(benchmark::State& state) {
    BM_RandomBits<double, uint64_t>(state, 1, 0x7FF0000000000000 - 1);
}

BENCHMARK(BM_RandomBits_double);
BENCHMARK(BM_RandomBits_double);
BENCHMARK(BM_RandomBits_double);

#endif

//----------------------------------------------------------------------------------------------------------
//
//----------------------------------------------------------------------------------------------------------
#if 0

template <typename Float>
static void BM_Const(benchmark::State& state, Float value) {
    for (auto _ : state) {
        char buffer[32];
        benchmark::DoNotOptimize(value);
        benchmark::DoNotOptimize(Dtoa(buffer, 32, value));
    }
}

#if 0
BENCHMARK_CAPTURE(BM_Const, double_1, 1.0);
BENCHMARK_CAPTURE(BM_Const, double_12, 12.0);
BENCHMARK_CAPTURE(BM_Const, double_123, 123.0);
BENCHMARK_CAPTURE(BM_Const, double_1234, 1234.0);
BENCHMARK_CAPTURE(BM_Const, double_12345, 12345.0);
BENCHMARK_CAPTURE(BM_Const, double_123456, 123456.0);
BENCHMARK_CAPTURE(BM_Const, double_1234567, 1234567.0);
BENCHMARK_CAPTURE(BM_Const, double_12345678, 12345678.0);
BENCHMARK_CAPTURE(BM_Const, double_123456789, 123456789.0);
BENCHMARK_CAPTURE(BM_Const, double_1234567890, 1234567890.0);
BENCHMARK_CAPTURE(BM_Const, double_12345678901, 12345678901.0);
BENCHMARK_CAPTURE(BM_Const, double_123456789012, 123456789012.0);
BENCHMARK_CAPTURE(BM_Const, double_1234567890123, 1234567890123.0);
BENCHMARK_CAPTURE(BM_Const, double_12345678901234, 12345678901234.0);
BENCHMARK_CAPTURE(BM_Const, double_123456789012345, 123456789012345.0);
#endif
#if 1
BENCHMARK_CAPTURE(BM_Const, double_1_2, 1.2);
BENCHMARK_CAPTURE(BM_Const, double_1_23, 1.23);
BENCHMARK_CAPTURE(BM_Const, double_1_234, 1.234);
BENCHMARK_CAPTURE(BM_Const, double_1_2345, 1.2345);
BENCHMARK_CAPTURE(BM_Const, double_1_23456, 1.23456);
BENCHMARK_CAPTURE(BM_Const, double_1_234567, 1.234567);
BENCHMARK_CAPTURE(BM_Const, double_1_2345678, 1.2345678);
BENCHMARK_CAPTURE(BM_Const, double_1_23456789, 1.23456789);
BENCHMARK_CAPTURE(BM_Const, double_1_234567895, 1.234567895); // 1.234567890 would be trimmed
BENCHMARK_CAPTURE(BM_Const, double_1_2345678901, 1.2345678901);
BENCHMARK_CAPTURE(BM_Const, double_1_23456789012, 1.23456789012);
BENCHMARK_CAPTURE(BM_Const, double_1_234567890123, 1.234567890123);
BENCHMARK_CAPTURE(BM_Const, double_1_2345678901234, 1.2345678901234);
BENCHMARK_CAPTURE(BM_Const, double_1_23456789012345, 1.23456789012345);
#endif

#endif

//----------------------------------------------------------------------------------------------------------
//
//----------------------------------------------------------------------------------------------------------
static const int64_t kPow10[] = {
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

#if 0

static void BM_Ints(benchmark::State& state, int digits) {
    constexpr int NumFloats = 1 << 12;

    std::uniform_int_distribution<int64_t> gen(kPow10[digits - 1], kPow10[digits] - 1);
    std::vector<double> numbers(NumFloats);
    std::generate(numbers.begin(), numbers.end(), [&] { return static_cast<double>(gen(random)); });

    int index = 0;
    for (auto _ : state) {
        char buffer[32];
        benchmark::DoNotOptimize(Dtoa(buffer, 32, numbers[index]));
        index = (index + 1) & (NumFloats - 1);
    }
}

BENCHMARK_CAPTURE(BM_Ints, _1,   1);
BENCHMARK_CAPTURE(BM_Ints, _2,   2);
BENCHMARK_CAPTURE(BM_Ints, _3,   3);
BENCHMARK_CAPTURE(BM_Ints, _4,   4);
BENCHMARK_CAPTURE(BM_Ints, _5,   5);
BENCHMARK_CAPTURE(BM_Ints, _6,   6);
BENCHMARK_CAPTURE(BM_Ints, _7,   7);
BENCHMARK_CAPTURE(BM_Ints, _8,   8);
BENCHMARK_CAPTURE(BM_Ints, _9,   9);
BENCHMARK_CAPTURE(BM_Ints, _10, 10);
BENCHMARK_CAPTURE(BM_Ints, _11, 11);
BENCHMARK_CAPTURE(BM_Ints, _12, 12);
BENCHMARK_CAPTURE(BM_Ints, _13, 13);
BENCHMARK_CAPTURE(BM_Ints, _14, 14);
BENCHMARK_CAPTURE(BM_Ints, _15, 15);

#endif

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------
#if 0

static void BM_Uniform(benchmark::State& state, double lower, double upper) {
    constexpr int const NumFloats = 1 << 12;

    const double Lo = std::pow(10.0, static_cast<double>(lower));
    const double Hi = std::pow(10.0, static_cast<double>(upper));

    std::uniform_real_distribution<double> gen(Lo, Hi);
    std::vector<double> numbers(NumFloats);
    std::generate(numbers.begin(), numbers.end(), [&] { return gen(random); });

    int index = 0;
    for (auto _ : state) {
        char buffer[64];
        benchmark::DoNotOptimize(Dtoa(buffer, 64, numbers[index]));
        index = (index + 1) & (NumFloats - 1);
    }
}

BENCHMARK_CAPTURE(BM_Uniform, _30_25, -30, -25);
BENCHMARK_CAPTURE(BM_Uniform, _25_20, -25, -20);
BENCHMARK_CAPTURE(BM_Uniform, _20_15, -20, -15);
BENCHMARK_CAPTURE(BM_Uniform, _15_20, -15, -20);
BENCHMARK_CAPTURE(BM_Uniform, _10_5,  -10,  -5);
BENCHMARK_CAPTURE(BM_Uniform, _5_1,    -5,  -1);
BENCHMARK_CAPTURE(BM_Uniform, _1_1,    -1,   1);
BENCHMARK_CAPTURE(BM_Uniform, _1_5,     1,   5);
BENCHMARK_CAPTURE(BM_Uniform, _5_10,    5,  10);
BENCHMARK_CAPTURE(BM_Uniform, _10_15,  10,  15);
BENCHMARK_CAPTURE(BM_Uniform, _15_20,  15,  20);
BENCHMARK_CAPTURE(BM_Uniform, _20_25,  20,  25);
BENCHMARK_CAPTURE(BM_Uniform, _25_30,  25,  30);
BENCHMARK_CAPTURE(BM_Uniform, _30_35,  30,  35);
BENCHMARK_CAPTURE(BM_Uniform, _35_40,  35,  40);
BENCHMARK_CAPTURE(BM_Uniform, _40_45,  40,  45);
BENCHMARK_CAPTURE(BM_Uniform, _45_50,  45,  50);

#endif

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------
#if 1

static void BM_RandomDigit(benchmark::State& state, int digits, double scale)
{
    assert(digits >= 1);
    assert(digits <= 15);

    constexpr int const NumFloats = 1 << 12;

    if (scale <= 0)
        scale = static_cast<double>(kPow10[digits - 1]);

    std::uniform_int_distribution<int64_t> gen(kPow10[digits - 1], kPow10[digits] - 1);
    std::vector<double> numbers(NumFloats);
    std::generate(numbers.begin(), numbers.end(), [&] {
        return static_cast<double>(gen(random) | 1) / scale;
    });

    int index = 0;
    for (auto _: state) {
        char buffer[32];
        benchmark::DoNotOptimize(Dtoa(buffer, 32, numbers[index]));
        index = (index + 1) & (NumFloats - 1);
    }
}

BENCHMARK_CAPTURE(BM_RandomDigit, _1_1,  1,  -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_2,  2,  -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_3,  3,  -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_4,  4,  -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_5,  5,  -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_6,  6,  -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_7,  7,  -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_8,  8,  -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_9,  9,  -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_10, 10, -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_11, 11, -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_12, 12, -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_13, 13, -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_14, 14, -1.0);
BENCHMARK_CAPTURE(BM_RandomDigit, _1_15, 15, -1.0);

BENCHMARK_CAPTURE(BM_RandomDigit,  _1_over_10e8,  1, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit,  _2_over_10e8,  2, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit,  _3_over_10e8,  3, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit,  _4_over_10e8,  4, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit,  _5_over_10e8,  5, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit,  _6_over_10e8,  6, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit,  _7_over_10e8,  7, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit,  _8_over_10e8,  8, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit,  _9_over_10e8,  9, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit, _10_over_10e8, 10, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit, _11_over_10e8, 11, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit, _12_over_10e8, 12, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit, _13_over_10e8, 13, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit, _14_over_10e8, 14, 1e8);
BENCHMARK_CAPTURE(BM_RandomDigit, _15_over_10e8, 15, 1e8);

#endif

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------
#if 1

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

BENCHMARK_MAIN()
